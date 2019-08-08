import numpy as np
from astropy.io import fits

try:
  import matplotlib.pyplot as plt
  from matplotlib import colors
  mpl=True
except:
  print("Warning: matplotlib not found, plotting will not work...")
  mpl=False

try:
  from pyds9 import DS9
  pyds9=True
except:
  pyds9=False
  
fitter=None
fitters=[]  
#Had issues getting iminuit to work, curve_fit works fine for now
try:
  import iminuit
  import probfit
  fitter="iminuit"
  fitters.append("iminuit")
  print("Using iminuit as fitter")
except:
  pass
try:
  from scipy.optimize import curve_fit
  if fitter is None:
    print("Iminuit not found, using scipy as fitter...")
    fitter="scipy"
    fitters.append("scipy")
except:
    print("Warning, no fitter found, cannot fit to noise...")


def mad(data):
  return np.median(np.abs(data-np.median(data)))


def plot_2d(data,cmap="rainbow",xlim=None, ylim=None, title=None, units="ADU", vmin=-50, vmax=None):
  if not mpl:
    return None, None
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
def gauss(x, norm,mean, sigma):
  return norm*np.exp(-(x-mean)**2/(2.*sigma**2))

#For fitting
def gauss_norm(x, mean, sigma):
  return 1./(sigma*(2*np.pi)**.5)*np.exp(-(x-mean)**2/(2.*sigma**2))

  
def fit_gauss(data, nbins=None, plot_fit=False, fitter=fitter, fit_range=None, rms_range_factor=8, *args, **kwargs):
  '''
  Wrapper function to perform binned fit of data to func. If given, will fit in the range "fit_range" (should be a tuple/array if [min,max]). Otherwise fits in the range median-rms_range_factor*MAD.
  '''
  if len(fitters) < 1:
    print("Error, no valid fitters loaded, exiting...")
    return None, None
  #Guess for initial parameters and to get a reasonable fitting range
  mean_guess=np.median(data)
  rms_guess=mad(data)
  if fit_range is None:
    fit_range=[mean_guess-rms_guess*rms_range_factor, mean_guess+rms_guess*rms_range_factor]
  fit_data=data[(data > fit_range[0]) & (data < fit_range[1])]

  #Get a reasonable number of bins
  if nbins is None:
    nbins=int((fit_range[1]-fit_range[0])/2)
    print("Bins for fit: " + str(nbins))
    if nbins < 40:
      print("Warning, number of bins is fairly small, fit may not be good...")

  if fitter not in fitters:
    print("Error, selected fitter is not loaded, loading another fitter")
    fitter=fitters[0]
      
  if fitter=="scipy":
    hist, bin_edges=np.histogram(fit_data, density=True, bins=nbins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist)
  
    if plot_fit and mpl:
      hist_fit = gauss(bin_centres, *coeff)
      fig=plt.figure()
      plt.plot(bin_centres, hist, label="Data")
      plt.plot(bin_centres, hist_fit, label="Fit")
      plt.show(False)
    #To keep compatibility between scipy fitter (which requires norm parameter A) and iminuit (which does not),
    #don't return the normalization coeff.
    coeff=coeff[1:]
    var_matrix=var_matrix[1:,1:]
    if plot_fit:
      return coeff, var_matrix, fig
    else:
      return coeff, var_matrix

  elif fitter=="iminuit":
    #ext_gauss=probfit.Extended(probfit.gaussian)
    #likelihood=probfit.BinnedLH(ext_gauss, fit_data,bins=nbins, extended=True)#, extended_bound=(-1000,1000))
    likelihood=probfit.BinnedLH(probfit.gaussian, fit_data,bins=nbins)#, extended_bound=(-1000,1000))
    #likelihood=probfit.UnbinnedLH(probfit.gaussian, fit_data, extended=False)
    minuit = iminuit.Minuit(likelihood,mean=mean_guess,sigma=rms_guess, print_level=1)
    minuit.migrad()
    #      minuit.minos()
    if plot_fit:
      fig=plt.figure()
      likelihood.draw(minuit)
    #Probably a better way to handle this
    coeff=[]
    coeff.append(minuit.values[0])
    coeff.append(minuit.values[1])
    var_matrix=[]
    var_matrix.append(minuit.matrix()[0])
    var_matrix.append(minuit.matrix()[1])
#    var_matrix.append(minuit.matrix()[2])
    if plot_fit:
      return coeff, var_matrix, fig
    else:
      return coeff, var_matrix
  else:
    print("Error, fitter not valid, exiting without doing fit...")
    return None, None
    
class SkipperImage:
  def _reset_derived(self):
    self.image_means=None
    self.image_rmses=None
    self.baseline=None
    self.charge_mask=None
    
    self.baseline_subtracted=False
    self.noise=None
    self.noise_coeffs=None
    self.noise_errs=None
    self.gen_charge=None
    self.image_gen_charge=None
    return True
    
  def __init__(self, fname, extension=0):
    self.hdu=fits.open(fname)[extension]
    self.header=self.hdu.header.copy()
    self.ndcms=self.header.get("NDCMS")
    self.nrows=self.header.get("NAXIS2")
    self.ncols=self.header.get("NAXIS1")
    self.ncols_phys=int(self.ncols/self.ndcms)
    self.npix=self.nrows*self.ncols_phys
    self.data=self.hdu.data
    #Total exposure time in seconds
    self.exptot=self.header.get("MEXP")/1000+self.header.get("MREAD")/1000
    self.nskips=self.ndcms
    #Default parameters
    self.pre_skips=0
    self.post_skips=0
    self.invert=False

    #Derived values    
    self._reset_derived()

    self.figures=[]
    self.ds9=None
    
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


    if changed>0:
      #All our old computed values become useless, so clear them out
      self._reset_derived()
      
    return changed
    
  def combine_skips(self,recompute=False,*args, **kwargs):
    '''
    Creates an image from the raw skipper data, skipping the first "pre_skips" and last "post_skip" skips.
    '''
    new_params=self.set_params(*args, **kwargs)

    if new_params==0 and self.image_means is not None and not recompute:
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
    new_image=self.compute_baseline(*args, **kwargs)
    if self.baseline_subtracted and not new_image:
      return False

    self.image_means-=self.baseline
    self.baseline_subtracted=True
    return True

  def draw_image(self, cmap="spectral", *args,**kwargs):
    '''
    Function to draw the image after combining skips.
    '''
    
    if not mpl:
      print("Error, matplotlib not found, cannot draw image...")
      return False

    self.combine_skips(*args, **kwargs)
    if pyds9:
      if self.ds9 is None:
        self.ds9=DS9()
      print("Drawing to DS9")
      self.ds9.set_np2arr(self.image_means)
      self.ds9.set("scale histequ")
      return True
    else:
      print("Drawing using matplotlib")
      fig, mesh=plot_2d(self.image_means)
      self.figures.append(fig)
      return True
    
  def compute_charge_mask(self, *args, **kwargs):
    '''
    Not implemented
    '''
    self.compute_baseline(*args, **kwargs)
    
    if self.baseline_subtracted:
      pass
      
  
  def compute_noise(self,nbins=None, plot_fit=False,*args, **kwargs):
    '''
    Performs a binned fit to the data to estimate the pixel noise of the CCD image
    '''

    #self.combine_skips(*args, **kwargs)
    #Need to subtract the baseline for the fit to work properly
    #TODO: fix that so this doesn't need to be done (maybe split left/right sides of image? Estimate mu?)
    self.subtract_baseline(*args, **kwargs)
    plot_fit=plot_fit and mpl

    fit_results=fit_gauss(self.image_means.flatten(), nbins=nbins, plot_fit=plot_fit, *args, **kwargs)
    coeff=fit_results[0]
    var_matrix=fit_results[1]
    if plot_fit:
      self.figures.append(fit_results[2])
    self.noise=coeff[-1]
    self.noise_coeffs=coeff
    self.noise_errs=var_matrix

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
    
