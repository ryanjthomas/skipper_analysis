#!/usr/bin/python

import numpy as np
from astropy.io import fits

import sys


if __name__=="__main__":
  invert=True
  skipskips=0
  if len(sys.argv) < 2:
    sys.exit(1)
  fname=sys.argv[1]
  hdu=fits.open(fname)
  header=hdu[0].header.copy()
  nskips=header.get("NDCMS")
  nrows=header.get("NAXIS2")
  ncols=header.get("NAXIS1")
  data=hdu[0].data
  means=np.zeros((nrows,ncols/nskips))
  rmses=np.zeros((nrows,ncols/nskips))

  if skipskips >=nskips:
    print("Warning: skipskips >= nskips, reducing skipskips")
    skipskips=nskips-1
  for i in range(ncols/nskips):
    xslice=np.s_[i*nskips+skipskips:(i+1)*nskips]
    means[:,i]=np.mean(data[:,i*nskips+skipskips:(i+1)*nskips],axis=1)
    rmses[:,i]=np.std(data[:,i*nskips+skipskips:(i+1)*nskips],axis=1)
  header["NAXIS1"]= ncols/nskips
  if invert:
    means=-1*means
  mean_hdu=fits.PrimaryHDU(means)
  mean_hdu.header=header
  rms_hdu=fits.ImageHDU(rmses, header, "RMSES")
  hdul=fits.HDUList([mean_hdu, rms_hdu])
  hdul.writeto("test.fits", overwrite=True)
  
