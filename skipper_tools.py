from numpy import np
from astropy.io import fits

class SkipperImage:
  def __init__(self, fname, extension=0):
    self.hdu=fits.open(fname)[extension]
    self.header=hdu.header.copy()
    self.ndcms=self.header.get("NDCMS")
    self.nrows=self.header.get("NAXIS2")
    self.ncols=self.header.get("NAXIS1")
    self.data=self.hdu.data
    #Total exposure time in seconds
    self.exptot=self.header.get("MEXP")/1000+self.header.get("MREAD")/1000

    self.image_means=None
    self.image_rmses=None
    #Default parameters
    self.pre_skips=0
    self.post_skips=0
    self.invert=False
    self.nskips=self.ndcms

  def set_params(nskips=None, pre_skips=None, post_skips=None, invert=None):
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
      self.pre_skips=int(self.ndcms*self.pre_skips/(total_skips))-1
      self.post_skips=int(self.ndcms*self.post_skips(total_skips))-1
      print("New pre_skips is: " + str(self.pre_skips) + ". New post_skips is: " +str(self.post_skips))

    return changed
    
  def combine_skips(*args, **kwargs)
    '''
    Creates an image from the raw skipper data, skipping the first "pre_skips" and last "post_skip" skips.
    '''
    self.set_params(*args, **kwargs)

    if self.pre_skips=pre_skips and self.post_skips=post_skips and self.image_means is not None:
      return

    #Actually make the image
    for i in range(self.ncols/self.ndcms):
      xslice=np.s_[i*self.ndcms+pre_skips:(i+1)*ndcms-post_skips]
      means[:,i]=np.mean(data[:,xslice],axis=1)
      rmses[:,i]=np.std(data[:,xslice],axis=1)

    if invert:
      means=-1*means
    #Store the result
    self.image_means=means
    self.image_rmses=rmses
    self.pre_skips=pre_skips
    self.post_skips=post_skips
    self.invert=invert
    return

  def write_combined_fits(fname, *args, **kwargs)
    '''
    Write the combined skips to a .fits file.
    '''

    self.set_params(*args, **kwargs)
    
    header=self.header.copy()
    header["NAXIS1"]=self.ncols/self.ndcms

    self.combine_skips(pre_skips, post_skips, invert)

    mean_hdu=fits.PrimaryHDU(self.means)
    mean_hdu.header=header
    rms_hdu=fits.ImageHDU(self.rmses, header, "RMSES")

    hdul=fits.HDUList([mean_hdu, rms_hdu])
    hdul.writeto(fname,overwrite=False)
    
    
