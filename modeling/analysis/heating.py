#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating Contains the DustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .component import AnalysisComponent

# -----------------------------------------------------------------

class DustHeatingAnalyser(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DustHeatingAnalyser, self).__init__(config)

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new HeatingAnalyser instance
        analyser = cls()

        # Set the modeling path
        analyser.config.path = arguments.path

        # Return the new instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        outpath = "modelChecks/iteration5_J14/"
        inpath = "SKIRTOutput/iteration5_J14/"
        inSED = "M31_212full_i77.5_sed.dat"

        Lsun = 3.846e26  # Watts

        # load SEDs
        input = np.loadtxt(inpath + inSED)
        wavelengths = input[:, 0]

        # Load the widths of the wavelength bins. Crucial for integration!
        delta_wls = np.loadtxt("SKIRTOutput/iteration5_J14/M31_reference_wavelengths512.dat", usecols=(1,))

        # only consider wavelengths longwards of 10 micron. For speed and memory
        startwl = next(wl[0] for wl in enumerate(wavelengths) if wl[1] > 10.)
        coldstartwl = next(wl[0] for wl in enumerate(wavelengths) if wl[1] > 100.)

        # produce wavelength ranges for all, warm and cold dust.
        dustwls = wavelengths[startwl:]
        warmwls = wavelengths[startwl:coldstartwl - 1]
        coldwls = wavelengths[coldstartwl:]
        delta_wls = delta_wls[startwl:]

        # Compute global heating fracions

        flux_all = input[startwl:, 1]

        input = np.loadtxt(inpath + inSED.replace('_i', '_old_i'))
        flux_old = input[startwl:, 1]

        input = np.loadtxt(inpath + inSED.replace('_i', '_young_i'))
        flux_young = input[startwl:, 1]

        Fold = 100. * (0.5 * flux_old + 0.5 * (flux_all - flux_young)) / flux_all
        Fyoung = 100. * (0.5 * flux_young + 0.5 * (flux_all - flux_old)) / flux_all

        Fold_alternative1 = 100. * flux_old / (flux_old + flux_young)
        Fyoung_alternative1 = 100. * flux_young / (flux_old + flux_young)

        Fold_alternative2 = 100. * flux_old / flux_all
        Fyoung_alternative2 = 100. * flux_young / flux_all

        Fold_alternative3 = 100. * np.sqrt(flux_old * (flux_all - flux_young)) / flux_all
        Fyoung_alternative3 = 100. * np.sqrt(flux_young * (flux_all - flux_old)) / flux_all

        Fold_alternative4 = 100. * (0.5 * flux_old + 0.5 * (flux_all - flux_young)) / (flux_old + flux_young)
        Fyoung_alternative4 = 100. * (0.5 * flux_young + 0.5 * (flux_all - flux_old)) / (flux_old + flux_young)

        JyToLsun = 1.e-26 * 4 * np.pi * (0.785e6 * 3.086e+16) ** 2 * 3.e14 / (dustwls ** 2) / Lsun  # in Lsun/micron

        totFlux = 0
        totFlux_young = 0
        totFlux_old = 0
        for i in range(len(flux_all) - 1):
            totFlux += delta_wls[i] * flux_all[i] * JyToLsun[i]
            totFlux_young += delta_wls[i] * flux_young[i] * JyToLsun[i]
            totFlux_old += delta_wls[i] * flux_old[i] * JyToLsun[i]

        print('Total heating from old stars: ', totFlux_old / totFlux)
        print('Total heating from young stars: ', totFlux_young / totFlux)

        plt.figure(figsize=(7, 5))
        plt.ylabel('$F^\prime_{\lambda,\mathrm{unev.}}$ [$\%$]', fontsize=20)
        plt.xlabel('$\lambda/\mu\mathrm{m}$', fontsize=20)
        plt.xlim(10., 1.e3)
        plt.ylim(0., 60.)
        plt.xscale('log')
        plt.tick_params(labelsize=20)
        # plt.subplots_adjust(bottom=0.1)
        # plt.plot(dustwls,Fold, 'r-', label="old stars")
        plt.plot(dustwls, Fyoung, 'k-', label="Young SPs")
        # plt.plot(dustwls,Fyoung_alternative1, 'r-', label="alt 1")
        # plt.plot(dustwls,Fyoung_alternative2, 'g-', label="alt 2")
        # plt.plot(dustwls,Fyoung_alternative3, 'c-', label="alt 3")
        plt.plot(dustwls, Fyoung_alternative4, 'k-', label="alt 4")
        plt.fill_between(dustwls, Fyoung, Fyoung_alternative4, color='grey', alpha='0.5')

        plt.tight_layout()
        # plt.legend(loc='upper left',numpoints=1,markerscale=1.5,fontsize=14)
        plt.savefig(outpath + inSED.replace('sed.dat', 'heating.pdf'), format='pdf')

        # plt.show()
        plt.close()

        # Make heating map
        inCube = inSED.replace('sed.dat', 'total.fits')
        makeHeatMap(inpath, outpath, inCube, startwl, dustwls, delta_wls)
        # makeWarmHeatMap(inpath,outpath,inCube,startwl,coldstartwl-1, warmwls)
        # makeColdHeatMap(inpath,outpath,inCube,coldstartwl, coldwls)

# -----------------------------------------------------------------

def makeHeatMap(inpath,outpath,inCube,startwl,dustwls, delta_wls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:,0:,0:]


    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    hdu = pyfits.PrimaryHDU( Fold,hdr_all)
    hdu.writeto(outpath+"heatingFold.fits",clobber=True)

    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all
    hdu = pyfits.PrimaryHDU( Fyoung,hdr_all)
    hdu.writeto(outpath+"heatingFyoung.fits",clobber=True)


    pixelTot = integratePixelSEDs(cube_all,     dustwls, delta_wls)
    pixelOld = integratePixelSEDs(cube_old,     dustwls, delta_wls)
    pixelYoung = integratePixelSEDs(cube_young, dustwls, delta_wls)

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    hdu = pyfits.PrimaryHDU(pixelTot,hdr_wcs)
    hdu.writeto(outpath+"Ldust_tot.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelOld,hdr_wcs)
    hdu.writeto(outpath+"Ldust_old.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelYoung,hdr_wcs)
    hdu.writeto(outpath+"Ldust_young.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelOld/pixelTot,hdr_wcs)
    hdu.writeto(outpath+"heatingTotOld.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelYoung/pixelTot,hdr_wcs)
    hdu.writeto(outpath+"heatingTotYoung.fits",clobber=True)


# OLD AND INCORRECT?
#    tot_all = 100.* (dustwls[len(dustwls)-1] - dustwls[0])
#
#    tot_old   = integrateHeating(Fold,dustwls) / tot_all
#    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
#    hdu.writeto(outpath+"heatingTotOld.fits",clobber=True)
#
#    tot_young = integrateHeating(Fyoung,dustwls) / tot_all
#    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
#    hdu.writeto(outpath+"heatingTotYoung.fits",clobber=True)

def makeWarmHeatMap(inpath,outpath,inCube,startwl,stopwl,warmwls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:stopwl,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:stopwl,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:stopwl,0:,0:]

    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all

    tot_all = 100.* (warmwls[len(warmwls)-1] - warmwls[0])

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    tot_old   = integrateHeating(Fold,warmwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
    hdu.writeto(outpath+"heatingTotWarmOld.fits",clobber=True)

    tot_young = integrateHeating(Fyoung,warmwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
    hdu.writeto(outpath+"heatingTotWarmYoung.fits",clobber=True)

def makeColdHeatMap(inpath,outpath,inCube,startwl,coldwls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:,0:,0:]


    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all

    tot_all = 100.* (coldwls[len(coldwls)-1] - coldwls[0])

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    tot_old   = integrateHeating(Fold,coldwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
    hdu.writeto(outpath+"heatingTotColdOld.fits",clobber=True)

    tot_young = integrateHeating(Fyoung,coldwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
    hdu.writeto(outpath+"heatingTotColdYoung.fits",clobber=True)


def integratePixelSEDs(cube, wls, dwls):


    Lsun = 3.846e26 # Watts
    MjySr_to_LsunMicron = 1.e6 * (36./206264.806247)**2 * 1.e-26 * 4*np.pi*(0.785e6*3.086e+16)**2 * 3.e14/(wls**2) / Lsun

    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])

    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = np.sum(sed * MjySr_to_LsunMicron * dwls)

    return slice


def integrateHeating(cube,dustwls):

    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])

    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = integrate.simps(sed,dustwls)

    return slice

# -----------------------------------------------------------------
