#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_fir_colours

# -----------------------------------------------------------------

import astropy.io.fits as pyfits
import numpy as np

def main():

    obsPath = "SKIRTInput/"
    outPath = "modelChecks/iteration5_J14/"
    SNR = 2.
    
    # Check observed FIR colours
    print 'Making observed FIR colour maps.'
    
    
    ####### 70 micron #######
    hdulist = pyfits.open(obsPath+"new70MJySr.fits")
    ima70 = hdulist[0].data
    hdr70 = hdulist[0].header
    # Set S/N cut
    hdulist = pyfits.open(obsPath+"new70_error.fits")
    err70 = hdulist[0].data
    ima70[err70 > 1./SNR] = np.nan


    ####### 100 micron #######
    hdulist = pyfits.open(obsPath+"new100MJySr.fits")
    ima100 = hdulist[0].data
    hdr100 = hdulist[0].header

    hdulist = pyfits.open(obsPath+"new100_error.fits")
    err100 = hdulist[0].data
    ima100[err100 > 1./SNR] = np.nan


    ####### 160 micron #######
    hdulist = pyfits.open(obsPath+"new160MJySr.fits")
    ima160 = hdulist[0].data
    hdr160 = hdulist[0].header

    hdulist = pyfits.open(obsPath+"new160_error.fits")
    err160 = hdulist[0].data
    ima160[err160 > 1./SNR] = np.nan


    ####### 250 micron #######
    hdulist = pyfits.open(obsPath+"new250MJySr.fits")
    ima250 = hdulist[0].data
    hdr250 = hdulist[0].header

    hdulist = pyfits.open(obsPath+"new250_error.fits")
    err250 = hdulist[0].data
    ima250[err250 > 1./SNR] = np.nan

    ####### 350 micron #######
    hdulist = pyfits.open(obsPath+"new350MJySr.fits")
    ima350 = hdulist[0].data
    hdr350 = hdulist[0].header
    
    hdulist = pyfits.open(obsPath+"new350_error.fits")
    err350 = hdulist[0].data
    ima350[err350 > 1./SNR] = np.nan


    ####### 500 micron #######
    hdulist = pyfits.open(obsPath+"new500MJySr.fits")
    ima500 = hdulist[0].data
    hdr500 = hdulist[0].header
    
    hdulist = pyfits.open(obsPath+"new500_error.fits")
    err500 = hdulist[0].data
    ima500[err500 > 1./SNR] = np.nan


    ####### 70/100 colour #######
    col70_100 = np.log10(ima70/ima100)
    hdu = pyfits.PrimaryHDU(col70_100,hdr70)
    hdu.writeto(obsPath+"obs70_100.fits",clobber=True)

    ####### 100/160 colour #######
    col100_160 = np.log10(ima100/ima160)
    hdu = pyfits.PrimaryHDU(col100_160,hdr100)
    hdu.writeto(obsPath+"obs100_160.fits",clobber=True)

    ####### 160/250 colour #######
    col160_250 = np.log10(ima160/ima250)
    hdu = pyfits.PrimaryHDU(col160_250,hdr160)
    hdu.writeto(obsPath+"obs160_250.fits",clobber=True)

    ####### 250/350 colour #######
    col250_350 = np.log10(ima250/ima350)
    hdu = pyfits.PrimaryHDU(col250_350,hdr250)
    hdu.writeto(obsPath+"obs250_350.fits",clobber=True)
    
    ####### 350/500 colour #######
    col350_500 = np.log10(ima350/ima500)
    hdu = pyfits.PrimaryHDU(col350_500,hdr350)
    hdu.writeto(obsPath+"obs350_500.fits",clobber=True)


    # Check model FIR colours
    print 'Making model FIR colour maps.'
    
    hdulist = pyfits.open(outPath+"modall70.fits")
    im70 =  hdulist[0].data[0:,0:]
    hdr70 = hdulist[0].header
    
    hdulist = pyfits.open(outPath+"modall100.fits")
    im100 =  hdulist[0].data[0:,0:]
    hdr100 = hdulist[0].header
    
    hdulist = pyfits.open(outPath+"modall160.fits")
    im160 =  hdulist[0].data[0:,0:]
    hdr160 = hdulist[0].header
    
    hdulist = pyfits.open(outPath+"modall250.fits")
    im250 =  hdulist[0].data[0:,0:]
    hdr250 = hdulist[0].header
    
    hdulist = pyfits.open(outPath+"modall350.fits")
    im350 =  hdulist[0].data[0:,0:]
    hdr350 = hdulist[0].header

    hdulist = pyfits.open(outPath+"modall500.fits")
    im500 =  hdulist[0].data[0:,0:]
    hdr500 = hdulist[0].header


    mod70_100  = np.log10(im70/im100)
    mask = np.isnan(col70_100)
    mod70_100[mask] = np.nan

    mod100_160 = np.log10(im100/im160)
    mask = np.isnan(col100_160)
    mod100_160[mask] = np.nan

    mod160_250 = np.log10(im160/im250)
    mask = np.isnan(col160_250)
    mod160_250[mask] = np.nan

    mod250_350 = np.log10(im250/im350)
    mask = np.isnan(col250_350)
    mod250_350[mask] = np.nan

    mod350_500 = np.log10(im350/im500)
    mask = np.isnan(col350_500)
    mod350_500[mask] = np.nan


    hdu = pyfits.PrimaryHDU(mod70_100,hdr70)
    hdu.writeto(outPath+"modall70_100.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(mod100_160,hdr100)
    hdu.writeto(outPath+"modall100_160.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(mod160_250,hdr160)
    hdu.writeto(outPath+"modall160_250.fits",clobber=True)
    
    hdu = pyfits.PrimaryHDU(mod250_350,hdr250)
    hdu.writeto(outPath+"modall250_350.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(mod350_500,hdr350)
    hdu.writeto(outPath+"modall350_500.fits",clobber=True)


    res160_250 = (col160_250 - mod160_250) / col160_250
    res250_350 = (col250_350 - mod250_350) / col250_350

    print "Median 160/250 residual: "+str(np.nanmedian(np.abs(res160_250)))
    print "Standard deviation 160/250 residual: "+str(np.nanstd(res160_250))

    print "Median 250/350 residual: "+str(np.nanmedian(np.abs(res250_350)))
    print "Standard deviation 250/350 residual: "+str(np.nanstd(res250_350))


def cropImage(im, hdr, region):
    xmin, ymin, xmax, ymax = readBox(region)
    cim = im[ymin-1:ymax-1,xmin-1:xmax-1]
    hdr['CRPIX1'] = hdr['CRPIX1'] - xmin + 1
    hdr['CRPIX2'] = hdr['CRPIX2'] - ymin + 1
    return cim, hdr

def readBox(region):
    with open(region) as f:
        for _ in xrange(3):
            line = f.readline()
        line = f.readline().replace('(',',').replace(')',',').split(',')
    
    xc = float(line[1])
    yc = float(line[2])
    xsize = float(line[3])
    ysize = float(line[4])

    return int(round(xc-xsize/2.)), int(round(yc-ysize/2.)), int(round(xc+xsize/2.)), int(round(yc+ysize/2.))


if __name__ == '__main__':
    main()