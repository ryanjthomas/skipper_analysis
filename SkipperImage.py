import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib import colors

from scipy.optimize import curve_fit

def mad(data):
  return np.median(np.abs(data))


def plot_2d(data,cmap="rainbow",xlim=None, ylim=None, title=None, units="ADU", vmin=-50, vmax=None):
  data=np.array(data)
  fig=plt.figure()
  # mesh=plt.pcolormesh(data-np.min(data),cmap=plt.get_cmap(cmap), norm=colors.LogNorm())
  mesh=plt.pcolormesh(data,cmap=plt.get_cmap(cmap), norm=colors.SymLogNorm(10, vmin=vmin, vmax=vmax))
    
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

#For fitting
def gauss(x, A,mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def fit_func(func, data, nbins=5000, plot_fit=False):
  '''
  Simple wrapper function to perform binned fit of data to func.
  '''
  hist, bin_edges=np.histogram(data, density=True, bins=nbins)
  bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
  coeff, var_matrix = curve_fit(func, bin_centres, hist)
  
  if plot_fit:
    hist_fit = func(bin_centres, *coeff)
    plt.figure()
    plt.plot(bin_centres, hist, label="Data")
    plt.plot(bin_centres, hist_fit, label="Fit")
    plt.show(False)

  return coeff, var_matrix

class SkipperImage:
  def __init__(self, fname, extension=0):
    self.hdu=fits.open(fname)[extension]
    self.header=self.hdu.header.copy()
    self.ndcms=self.header.get("NDCMS")
    self.nrows=self.header.get("NAXIS2")
    self.ncols=self.header.get("NAXIS1")
    self.ncols_phys=int(self.ncols/self.ndcms)
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

    #Derived values
    self.nskips=self.ndcms
    self.baseline_subtracted=False
    self.noise=None
    self.sigma_coeffs=None
    self.sigma_errs=None
    self.gen_charge=None
    self.image_gen_charge=None
    
  def set_params(self,nskips=None, pre_skips=None, post_skips=None, invert=None, *args, **kwars):
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

    means=np.zeros((self.nrows,self.ncols_phys))
    rmses=np.zeros((self.nrows,self.ncols_phys))

    #Actually make the image
    for i in range(self.ncols_phys):
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

    header=self.header.copy()
    header["NAXIS1"]=self.ncols/self.ndcms

    self.combine_skips(*args,**kwargs)

    mean_hdu=fits.PrimaryHDU(self.image_means)
    mean_hdu.header=header
    rms_hdu=fits.ImageHDU(self.image_rmses, header, "RMSES")

    hdul=fits.HDUList([mean_hdu, rms_hdu])
    hdul.writeto(fname,clobber=False)
    
  def compute_baseline(self, *args, **kwargs):
    '''
    Computes the baseline for the CCD image, as the median value along each row and column. Also forms the baseline per pixel as the sum of the baselines in x and y.
    '''
    new_image=self.combine_skips(*args, **kwargs)

    if self.baseline is not None and not new_image:
      return False
    self.baseline=np.zeros(self.image_means.shape)

    self.baseline_x=np.median(self.image_means-self.baseline, axis=0)
    self.baseline[:,:]+=self.baseline_x[np.newaxis,:]
    self.baseline_y=np.median(self.image_means-self.baseline, axis=1)
    self.baseline[:,:]+=self.baseline_y[:,np.newaxis]
    return True
  
  def subtract_baseline(self, *args, **kwargs):
    '''
    Subtracts baseline from the image. Computes the baseline and skipper means, if not already done.
    '''
    new_image=self.combine_skips(*args, **kwargs)
    if self.baseline_subtracted and not new_image:
      return False
    self.compute_baseline(*args, **kwargs)
    self.image_means-=self.baseline
    self.baseline_subtracted=True
    return True

  def draw_image(self, cmap="spectral", *args,**kwargs):
    self.combine_skips(*args, **kwargs)
    plot_2d(self.image_means)

  def compute_charge_mask(self, *args, **kwargs):
    '''
    Not implemented
    '''
    self.compute_baseline(*args, **kwargs)
    
    if self.baseline_subtracted:
      pass
      
  
  def compute_noise(self,nbins=5000, plot_fit=False,*args, **kwargs):
    '''
    Performs a binned fit to the data to estimate the pixel noise of the CCD image
    '''

    #self.combine_skips(*args, **kwargs)
    #Need to subtract the baseline for the fit to work properly
    #TODO: fix that so this doesn't need to be done (maybe split left/right sides of image? Estimate mu?)
    self.subtract_baseline(*args, **kwargs)
    coeff, var_matrix=fit_func(gauss,self.image_means, nbins=5000, plot_fit=plot_fit)

    self.noise=coeff[2]
    self.sigma_coeffs=coeff
    self.sigma_errs=var_matrix

    return coeff, var_matrix

  def compute_gen_charge(self, *args, **kwargs):
    '''
    Estimates the amount of spurious charge generated in each pixel from the skipping procedure
    '''
    self.set_params(*args, **kwargs)
    self.image_gen_charge=np.zeros((self.nrows, self.ncols_phys))
    #TODO: reimplement this using slices to eliminate the loop
    for i in range(self.ncols_phys):
      x1=i*self.ndcms+self.pre_skips
      x2=(i+1)*self.ndcms-self.post_skips-1
      self.image_gen_charge[:,i]=self.data[:,x2]-self.data[:,x1]
    self.gen_charge=np.median(self.image_gen_charge)
    return self.gen_charge

  def compute_charge_loss(*args, **kwargs):
    self.set_params(*args, **kwargs)
      
  def compute_statistics(self, *args, **kwargs):
    self.combine_skips(*args, **kwargs)
    self.compute_gen_charge(*args, **kwargs)
    self.compute_charge_loss(*args, **kwargs)
    
