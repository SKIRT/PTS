#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_images

# -----------------------------------------------------------------

import astropy.io.fits as pyfits
import numpy as np
from scipy import interpolate
from scipy import integrate


def main():

    obsPath = "SKIRTInput/"
    modPath = "SKIRTOutput/iteration5_J14/"
    outPath = "modelChecks/iteration5_J14/h67/"
    
    # Total circumstellar dust derived from the MAPPINGS III SED of the ionizing stars
    TotalCircumDust = 2143279.799
    
    # bands = ['FUV','NUV','u','g','r','i','z','W1','W2','W3','W4','3.6','4.5','5.8','8','24','70','100','160','250','350','500']
    # simIDs = [ 0,    2,    3,  4,  5,  6,  7,  8,   10,  21,  24,  9,    10,   12,   15, 25,  26,  27,   28,   29,   30,   31]
    
    #bands  = ['FUV','NUV','r','W1','W2','W3','W4','3.6','4.5','5.8','8','24','350']
    #simIDs   = np.array([2,3,5,41,50,78,96,43,51,57,66,98,110])
    
    bands  = ['FUV','NUV','u','g','r','i','z','W1','W2','W3','W4','3.6','4.5','5.8','8','24','70','100','160','250','350','500']

    #datacube = modPath+"M31_212_i77.5_total.fits"
    #datacube = modPath+"M31_212full_i77.5_total.fits"
    #datacube = modPath+"M31_212_h50_i77.5_total.fits"
    datacube = modPath+"M31_212_h67_i77.5_total.fits"

    input       = np.loadtxt(datacube.replace('total.fits','sed.dat'))
    wavelengths = input[:,0]
    
    cube = pyfits.open(datacube)
    suffix = "all"
    #suffix = "Shear_"
    
    for i in range(0,len(bands)):
    
        hdulist = pyfits.open(obsPath+"new"+bands[i]+"MJySr.fits")
        obsIm =  hdulist[0].data[0:,0:]
        obsHdr = hdulist[0].header


        filter = "Files/Filters/transmission_"+bands[i]+".dat"
        modImage = convolveFilter2D(cube[0].data[0:,0:,0:],wavelengths,filter)

        obsIm[obsIm == 0] = np.nan
        #res = np.absolute((obsIm - mod)/mod)
        res = (obsIm - modImage)/obsIm
       
        # Set S/N cutoff for residuals
        SNR = 2.

        hdulist = pyfits.open(obsPath+"new"+bands[i]+"_error.fits")
        errIm =  hdulist[0].data[0:,0:]
        mask = errIm > 1./SNR
        res[mask] = np.nan

        if bands[i] == 'FUV':
            datacube = datacube.replace('total','transparent')
            transCube = pyfits.open(datacube)
            transHdr = transCube[0].header

            filter = "Files/Filters/transmission_FUV.dat"
            transImage = convolveFilter2D(transCube[0].data[0:,0:,0:],wavelengths,filter)

            A = -2.5*(np.log10(modImage) -np.log10(transImage))
            hdu = pyfits.PrimaryHDU(A,transHdr)
            hdu.writeto(outPath+"diffuseAtt"+suffix+bands[i]+".fits",clobber=True)

            # Load geometry for ionizing stars
            hdulist = pyfits.open(obsPath+"ionizingStars.fits")
            ionStarsMap = hdulist[0].data[0:,0:]

            # normalize geometry to total circumstellar dust
            ionStarsMap = ionStarsMap / np.sum(ionStarsMap) * TotalCircumDust
            print "Normalization of ionizing stars Geometry", np.sum(ionStarsMap)
            
            # write out circumstellar dust map
            hdu = pyfits.PrimaryHDU(ionStarsMap,hdulist[0].header)
            hdu.writeto(outPath+"CircumMdust.fits",clobber=True)
            
            #create dust mass surface density map assuming 137.008346269 pc as pixel scale
            SigmaCircumDust = ionStarsMap / 0.137**2
            
            # Convert to Attenuation following Kreckel et al 2013, equation (4)
            # And a Av to AFUV conversion factor of 8.17604E-01 / 2.23403E-01 from http://www.mpia.de/homes/brent/SBmodels/attenuation_data.txt
            circumAttFUV = 8.17604e-01 / 2.23403e-01 * 0.67* SigmaCircumDust * 1.e-5
            
            hdu = pyfits.PrimaryHDU(circumAttFUV,transHdr)
            hdu.writeto(outPath+"circumAtt"+suffix+bands[i]+".fits",clobber=True)
        
            # add this map to the FUV attenuation map from diffuse dust
            A = A+circumAttFUV
            
            hdu = pyfits.PrimaryHDU(A,transHdr)
            hdu.writeto(outPath+"Att"+suffix+bands[i]+".fits",clobber=True)
            
        
            # Compute residuals
            hdulist = pyfits.open(obsPath+"M31_A_FUV.fits")
            hdu = pyfits.PrimaryHDU( (hdulist[0].data[0:,0:] - A)/A ,hdulist[0].header)
            hdu.writeto(outPath+"res"+suffix+"AFUV.fits",clobber=True)
        
            hdulist = pyfits.open(obsPath+"M31_obs_A_FUV.fits")
            hdu = pyfits.PrimaryHDU( (hdulist[0].data[0:,0:] - A)/A ,hdulist[0].header)
            hdu.writeto(outPath+"res"+suffix+"obsAFUV.fits",clobber=True)

        
        
        
        hdu = pyfits.PrimaryHDU(res,obsHdr)
        hdu.writeto(outPath+"res"+suffix+bands[i]+".fits",clobber=True)
        
        hdu = pyfits.PrimaryHDU(modImage,obsHdr)
        hdu.writeto(outPath+"mod"+suffix+bands[i]+".fits",clobber=True)


def convolveFilter2D(cube,wl,filter):
    
    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])
    
    input = np.loadtxt(filter)
    response_wl     = input[:,0] * 1.e-4
    response_trans  = input[:,1]
    
    minwl = response_wl[0]
    maxwl = response_wl[len(response_wl)-1]
    intwl = np.copy(wl)
    intwl[intwl > maxwl] = -1
    intwl[intwl < minwl] = -1
    wlrange = intwl > 0
    intwl = intwl[wlrange]
    transmission = np.zeros(len(wl))
    
    interpfunc = interpolate.interp1d(response_wl,response_trans, kind='linear')
    transmission[wlrange] = interpfunc(intwl)
    tot_trans = integrate.simps(transmission,wl)
    
    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = integrate.simps(transmission*sed,wl)
    
    
    return slice/tot_trans

    #plt.plot(wl,transmission,'bo')
    #plt.plot(wl,transmission,'k-')
    #plt.plot(response_wl,response_trans,'r+')
    #plt.xscale('log')
    #plt.show()


if __name__ == '__main__':
    main()