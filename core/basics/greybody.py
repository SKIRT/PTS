#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.greybody Modified black body radiation function and fitting
#
# The class in this module implements the modified black body function representing the idealized radiation emitted by
# a body of dust at long wavelengths, and allows fitting the dust temperature and mass to a number of flux data points.

# -----------------------------------------------------------------

import numpy as np
from functools import partial
from scipy.optimize import curve_fit

# -----------------------------------------------------------------

# universal constants and units
c = 2.99792458e8        # m/s
h = 6.62606957e-34      # J s
k = 1.3806488e-23       # J/K
pc = 3.08567758e16      # m
Msun = 1.9891e30        # kg

# typical opacities at 350 micron
kappa350_Cortese = 0.192   # m2/kg (Cortese)
kappa350_Zubko = 0.330     # m2/kg (Zubko)

## This static function returns the black body emissivity \f$B_\nu\f$ (W/m2/Hz) for the specified wavelength or
# wavelengths \em wave (micron), and for the given dust temperature \em T (K).
def Bnu(wave, T):
    nu = c / (wave*1e-6)                                    # Hz
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    return Bnu

## This static function returns the grey body flux (Jy) for the specified wavelength or wavelengths \em wave (micron),
# for the given observer distance \em D (pc), power-law exponent \em beta, opacity \em kappa at 350 micron (m2/kg),
# dust temperature \em T (K) and dust mass \em M (Msun).
def greybody(D, beta, kappa350, wave, T, M):
    nu = c / (wave*1e-6)                                    # Hz
    nu350 = c / 350e-6                                      # Hz
    kappa = kappa350 * (nu/nu350)**beta                     # m2/kg
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    flux = M*Msun * kappa * Bnu / (D*pc)**2                 # W/m2/Hz
    return flux * 1e26                                      # Jy

# ----------------------------------------------------------------------

## This class represents a grey body with predefined material constants and distance.
class GreyBody:
    ## This constructor creates a grey body instance with bound values for the observer distance \em D (pc),
    # power-law exponent \em beta, and opacity \em kappa at 350 micron (m2/kg).
    def __init__(self, D, beta, kappa350):
        self._bound_greybody = partial(greybody, D, beta, kappa350)

    # This function serves to make the GreyBody instance callable. It returns the grey body flux (Jy) for the
    # specified wavelength or wavelengths \em wave (micron), dust temperature \em T (K) and dust mass \em M (Msun),
    # using the values for distance, power-law exponent, and opacity specified in the constructor of the object.
    def __call__(self, wave, T, M):
        return self._bound_greybody(wave, T, M)

    # This function returns the temperature \em T (K) and the dust mass \em M (Msun) for the best grey body fit with
    # the given flux data data points (wavelengths in micron and fluxes in Jy) and uncertainties (relative weights),
    # using the values for distance, power-law exponent, and opacity specified in the constructor of the object.
    # Optionally, the caller can specify initial values for the temperature (K) and/or dust mass (Msun).
    def fit(self, lambdav, fluxv, sigmav, Tinit=20, Minit=1e7):
        if np.all(np.asarray(fluxv)>0):
            popt, pcov = curve_fit(self._bound_greybody, lambdav, fluxv, p0=(Tinit,Minit), sigma=sigmav,
                                    absolute_sigma=False, maxfev=5000)
            return popt
        else:
            return (0,0)

# ----------------------------------------------------------------------
