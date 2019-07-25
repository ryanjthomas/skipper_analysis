import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib import colors

def mad(data):
  return np.median(np.abs(data))

def plot_2d(data,cmap="spectral",xlim=None, ylim=None, title=None, units="ADU"):
  data=np.array(data)
  fig=plt.figure()
  mesh=plt.pcolormesh(data-np.min(data),cmap=plt.get_cmap(cmap), norm=colors.LogNorm())
    
  ax=mesh.axes
  ax.set_aspect("equal")
  if xlim is not None:
    ax.set_xlim(right=xlim)
  if ylim is not None:
    ax.set_ylim(top=ylim)
  cbar=plt.colorbar()
  cbar.set_label(units)
  plt.xlabel("X"); plt.ylabel("Y");
  plt.title(title)
  plt.show(False)
  return fig, mesh

class SkipperImage:
  def __init__(self, fname, extension=0):
    self.hdu=fits.open(fname)[extension]
    self.header=self.hdu.header.copy()
    self.ndcms=self.header.get("NDCMS")
    self.nrows=self.header.get("NAXIS2")
    self.ncols=self.header.get("NAXIS1")
    self.data=self.hdu.data
    #Total exposure time in seconds
    self.exptot=self.header.get("MEXP")/1000+self.header.get("MREAD")/1000

    self.image_means=None
    self.image_rmses=None
    self.baseline=None
    self.charge_mask=None
    #Default parameters
    self.pre_skips=0
    self.post_skips=0
    self.invert=False
    self.nskips=self.ndcms
    self.baseline_subtracted=False
    
  def set_params(self,nskips=None, pre_skips=None, post_skips=None, invert=None):
    '''
    Simple handler function to change parameters for handling the skipper image. Returns number of changed parameters.
    '''
    changed=0
    
    if nskips is not None and nskips > 0 and nskips != self.nskips:
      self.nskips=nskips
      changed+=1
      
    if pre_skips is not None and pre_skips >= 0 and pre_skips != self.pre_skips:
      self.pre_skips=pre_skips
      changed+=1
      
    if post_skips is not None and post_skips >= 0 and post_skips != self.post_skips:
      self.post_skips=post_skips
      changed+=1

    if invert is not None and type(invert) is bool and invert != self.invert:
      self.invert=invert
      changed+=1

    total_skips=self.pre_skips+self.post_skips
    if total_skips >= self.ndcms:
      print("Warning, trying to skip more skips than exist, reducing skipped skips...")
      #Reduced both proportional to the original setting
      self.pre_skips=int((self.ndcms-1)*self.pre_skips/(total_skips))
      self.post_skips=int((self.ndcms-1)*self.post_skips/(total_skips))
      print("New pre_skips is: " + str(self.pre_skips) + ". New post_skips is: " +str(self.post_skips))

    return changed
    
  def combine_skips(self,force=False,*args, **kwargs):
    '''
    Creates an image from the raw skipper data, skipping the first "pre_skips" and last "post_skip" skips.
    '''
    new_params=self.set_params(*args, **kwargs)

    if new_params==0 and self.image_means is not None and not force:
      return False

    means=np.zeros((self.nrows,self.ncols/self.ndcms))
    rmses=np.zeros((self.nrows,self.ncols/self.ndcms))

    #Actually make the image
    for i in range(self.ncols/self.ndcms):
      xslice=np.s_[i*self.ndcms+self.pre_skips:(i+1)*self.ndcms-self.post_skips]
      means[:,i]=np.mean(self.data[:,xslice],axis=1)
      rmses[:,i]=np.std(self.data[:,xslice],axis=1)

    if self.invert:
      means=-1*means
    #Store the result
    self.image_means=means
    self.image_rmses=rmses
    self.baseline_subtracted=False

    return True

  def write_combined_fits(self,fname, *args, **kwargs):
    '''
    Write the combined skips to a .fits file.
    '''

    self.set_params(*args, **kwargs)
    
    header=self.header.copy()
    header["NAXIS1"]=self.ncols/self.ndcms

    self.combine_skips()

    mean_hdu=fits.PrimaryHDU(self.image_means)
    mean_hdu.header=header
    rms_hdu=fits.ImageHDU(self.image_rmses, header, "RMSES")

    hdul=fits.HDUList([mean_hdu, rms_hdu])
    hdul.writeto(fname,overwrite=False)
    
  def compute_baseline(self, *args, **kwargs):
    new_image=self.combine_skips(*args, **kwargs)

    if self.baseline is not None and new_image is False:
      return False
    self.baseline=np.zeros(self.image_means.shape)

    self.baseline_y=np.median(self.image_means, axis=1)
    self.baseline[:,:]+=self.baseline_y[:,np.newaxis]
    self.baseline_x=np.median(self.image_means-self.baseline, axis=0)
    self.baseline[:,:]+=self.baseline_x[np.newaxis,:]
    return True
  
  def subtract_baseline(self, *args, **kwargs):
    new_image=self.combine_skips(*args, **kwargs)
    self.compute_baseline()
    self.image_means-=self.baseline
    self.baseline_subtracted=True
    return True

  def draw_image(self, cmap="spectral", *args,**kwargs):
    self.combine_skips(*args, **kwargs)
    plot_2d(self.image_means)

  def compute_charge_mask(self, *args, **kwargs):
    self.compute_baseline(*args, **kwargs)
    
    if self.baseline_subtracted:
      pass
      
