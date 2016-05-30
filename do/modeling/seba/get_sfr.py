#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.get_sfr

# -----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import glob
import os
from scipy import interpolate
from scipy import integrate

def main():

    Lsun = 3.846e26 # Watts
    D = 0.785e6*3.086e+16 # Distance in meter
    
    # Wavelengths
    wls, delta_wls = np.loadtxt("SKIRTOutput/iteration5_J14/M31_reference_wavelengths512.dat", usecols=(0,1), unpack=True)
    
    # skirt best fit SED
    input    = np.loadtxt("SKIRTOutput/iteration5_J14/M31_212full_i77.5_sed.dat")
    allFlux = input[:,1] * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(wls*1.e-6) / Lsun
    
    # skirt best fit stellar SED. BEWARE OF THE MAPPING CONTRIBUTION IN THE FIR
    input       = np.loadtxt("modelChecks/iteration5_J14/nodust/M31_212full_stellar_i77.5_sed.dat")
    stellarFlux = input[:,1] * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(wls*1.e-6) / Lsun

    # Integrate total infrared
    # 3-1100 micron as defined in Kennicutt & Evans 2012, table 1
    idx = wls > 3.
    TIR = np.sum(allFlux[idx]/wls[idx] * delta_wls[idx])

    # Residuals in all bands
    bands  = ['FUV','NUV', 'u','g','r','i','z', 'W1','3.6','4.5','W2', '5.8','8',
              'W3', 'IRAS1', 'W4', '24','IRAS2', 'IRAS3', '70', '100', 'IRAS4', '160', '250','350','Planck350','500', 'Planck550', 'Planck850']
    
    modAllFlux = np.array([])
    modStarFlux = np.array([])
    for band in bands:
        filter = "Files/Filters/transmission_"+band+".dat"
        modAllFlux = np.append(modAllFlux,convolveFilter(allFlux,wls,filter))
        modStarFlux = np.append(modStarFlux,convolveFilter(stellarFlux,wls,filter))

    #for band, flux in zip(bands, modAllFlux):
        #print band+': %e Lsun' %flux


    #for band, flux in zip(bands, modStarFlux):
        #print band+': %e Lsun' %flux

    print 'SFR from TIR:       %f Msun/yr' % (10**(np.log10(1.e7 * Lsun * TIR) - 43.41))

    j =  bands.index('FUV')
    print 'SFR from FUV:       %f Msun/yr' % (10**(np.log10(1.e7 * Lsun * modStarFlux[j]) - 43.35))

    j =  bands.index('NUV')
    print 'SFR from NUV:       %f Msun/yr' % (10**(np.log10(1.e7 * Lsun * modStarFlux[j]) - 43.17))

    j =  bands.index('24')
    print 'SFR from 24 micron: %f Msun/yr' % (10**(np.log10(1.e7 * Lsun * modAllFlux[j]) - 42.69))

    j =  bands.index('70')
    print 'SFR from 70 micron: %f Msun/yr' % (10**(np.log10(1.e7 * Lsun * modAllFlux[j]) - 43.23))


def convolveFilter(flux,wl,filter):
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
    transmission = np.zeros(len(flux))
    
    interpfunc = interpolate.interp1d(response_wl,response_trans, kind='linear')
    transmission[wlrange] = interpfunc(intwl)
    tot_trans = integrate.simps(transmission,wl)
    tot_flux  = integrate.simps(transmission*flux,wl)
    
    
    return tot_flux/tot_trans


if __name__ == '__main__':
    main()