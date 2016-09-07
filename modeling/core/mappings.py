#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.mappings Contains the Mappings class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.sed import IntrinsicSED
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

# Determine the path to the Mappings SED directory
mappings_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "Mappings")

# -----------------------------------------------------------------

class Mappings(object):
    
    """
    This class ...
    """

    def __init__(self, metallicity, compactness, pressure, covering_factor, sfr=1.):

        """
        The constructor ...
        :param metallicity:
        :param compactness:
        :param pressure:
        :param covering_factor:
        """

        # Load the data files
        lambdav, _Zrelv, _logCv, _logpv, j0_dict, j1_dict = self.load_data()

        # Convert the input parameters to the parameters that are assumed in MAPPINGS III.
        # - the metallicity is converted from an absolute value Z to a value Zrel relative to the
        #   sun, where Zsun = 0.0122 as in Asplund et al. (2005). The same value is used in the MAPPINGS III
        #   models of Groves et al. (2008).
        # - the pressure is converted from the actual pressure in SI units (i.e. in Pa = N/m^2) to
        #   log(p/k), with k Boltzmann's constant, and in units of K/cm^3

        # Relative metallicity
        rel_metallicity = metallicity / 0.0122

        # Log(pressure)
        #pressure_in_Kcm3 = (pressure / boltzman).to("K/cm3").value # in SKIRT, the pressure is in program units (=SI units)
        pressure_in_Kcm3 = pressure.to("K/cm3").value
        log_p = np.log10(pressure_in_Kcm3)

        # Log(compactness)
        log_c = compactness  # compactness is already logarithmic

        # Covering factor
        fPDR = covering_factor

        # Ensure that the value of the parameters is within the boundaries of the parameter space
        rel_metallicity = max(rel_metallicity, 0.05)
        rel_metallicity = min(rel_metallicity, 2.0 - 1e-8)
        log_c = max(log_c, 4.0)
        log_c = min(log_c, 6.5 - 1e-8)
        log_p = max(log_p, 4.0)
        log_p = min(log_p, 8.0 - 1e-8)

        # Find the appropriate SED from interpolating in the library
        #i_ = find_nearest(np.array(_Zrelv), rel_metallicity)
        i = locate_clip(_Zrelv, rel_metallicity)
        #if i_ != i: print(i_, i)
        hZrel = (rel_metallicity-_Zrelv[i])/(_Zrelv[i+1]-_Zrelv[i])

        #j_ = find_nearest(np.array(_logCv), log_c)
        j = locate_clip(_logCv, log_c)
        #if j_ != j: print(j_, j)
        hlogC = (log_c - _logCv[j])/(_logCv[j+1]-_logCv[j])

        #k_ = find_nearest(np.array(_logpv), log_p)
        k = locate_clip(_logpv, log_p)
        #if k_ != k: print(k_, k)
        hlogp = (log_p - _logpv[k])/(_logpv[k+1]-_logpv[k])

        j0LLLv = j0_dict[i  , j  , k  ]
        j0RLLv = j0_dict[i+1, j  , k  ]
        j0LRLv = j0_dict[i  , j+1, k  ]
        j0RRLv = j0_dict[i+1, j+1, k  ]
        j0LLRv = j0_dict[i  , j  , k+1]
        j0RLRv = j0_dict[i+1, j  , k+1]
        j0LRRv = j0_dict[i  , j+1, k+1]
        j0RRRv = j0_dict[i+1, j+1, k+1]
        j1LLLv = j1_dict[i  , j  , k  ]
        j1RLLv = j1_dict[i+1, j  , k  ]
        j1LRLv = j1_dict[i  , j+1, k  ]
        j1RRLv = j1_dict[i+1, j+1, k  ]
        j1LLRv = j1_dict[i  , j  , k+1]
        j1RLRv = j1_dict[i+1, j  , k+1]
        j1LRRv = j1_dict[i  , j+1, k+1]
        j1RRRv = j1_dict[i+1, j+1, k+1]

        # Initialize jv
        jv = np.zeros(len(lambdav))

        # Fill jv
        for k in range(len(lambdav)):

            j0 = (1.0 - hZrel) * (1.0 - hlogC) * (1.0 - hlogp) * j0LLLv[k] \
             + hZrel * (1.0 - hlogC) * (1.0 - hlogp) * j0RLLv[k] \
             + (1.0 - hZrel) * hlogC * (1.0 - hlogp) * j0LRLv[k] \
             + hZrel * hlogC * (1.0 - hlogp) * j0RRLv[k] \
             + (1.0 - hZrel) * (1.0 - hlogC) * hlogp * j0LLRv[k] \
             + hZrel * (1.0 - hlogC) * hlogp * j0RLRv[k] \
             + (1.0 - hZrel) * hlogC * hlogp * j0LRRv[k] \
             + hZrel * hlogC * hlogp * j0RRRv[k]

            j1 = (1.0 - hZrel) * (1.0 - hlogC) * (1.0 - hlogp) * j1LLLv[k] \
             + hZrel * (1.0 - hlogC) * (1.0 - hlogp) * j1RLLv[k] \
             + (1.0 - hZrel) * hlogC * (1.0 - hlogp) * j1LRLv[k] \
             + hZrel * hlogC * (1.0 - hlogp) * j1RRLv[k] \
             + (1.0 - hZrel) * (1.0 - hlogC) * hlogp * j1LLRv[k] \
             + hZrel * (1.0 - hlogC) * hlogp * j1RLRv[k] \
             + (1.0 - hZrel) * hlogC * hlogp * j1LRRv[k] \
             + hZrel * hlogC * hlogp * j1RRRv[k]

            jk = (1.0 - fPDR) * j0 + fPDR * j1
            jv[k] = jk

        wavelength_column = lambdav
        luminosity_column = jv * sfr

        # Create the SED
        self.sed = IntrinsicSED()
        self.sed.table = tables.new([wavelength_column, luminosity_column], ["Wavelength", "Luminosity"])
        self.sed.table["Wavelength"].unit = Unit("micron")
        self.sed.table["Luminosity"].unit = Unit("W/micron")

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # number of items in the library read by the constructor
        Nlambda = 1800
        NZrel = 5
        NlogC = 6
        Nlogp = 5

        #_lambdav.resize(Nlambda);
        #_lambdav = [None] * Nlambda
        _lambdav = None

        _Zrelv = [0.05, 0.20, 0.40, 1., 2.]
        Zrelnamev = ["Z005", "Z020", "Z040", "Z100", "Z200"]

        # vector<QString> logCnamev(NlogC);
        _logCv = [4., 4.5, 5., 5.5, 6., 6.5]
        logCnamev = ["C40", "C45", "C50", "C55", "C60", "C65"]

        # vector<QString> logpnamev(Nlogp);
        _logpv = [4., 5., 6., 7., 8.]
        logpnamev = ["p4", "p5", "p6", "p7", "p8"]

        # The j0 and j1 values
        j0_dict = dict()
        j1_dict = dict()

        # Read in the emissivity vectors
        for i in range(NZrel):
            for j in range(NlogC):
                for k in range(Nlogp):

                    filename = "Mappings_" + Zrelnamev[i] + "_" + logCnamev[j] + "_" + logpnamev[k] + ".dat"
                    path = fs.join(mappings_path, filename)

                    # Check whether the file exists
                    if not fs.is_file(path): raise IOError("The file '" + path + "' does not exist")

                    wavelengths, j0, j1 = np.loadtxt(path, unpack=True)
                    Zrel = _Zrelv[i]
                    logC = _logCv[j]
                    logp = _logpv[k]
                    #j0_dict[(Zrel, logC, logp)] = j0
                    #j1_dict[(Zrel, logC, logp)] = j1

                    j0_dict[(i, j, k)] = j0
                    j1_dict[(i, j, k)] = j1

                    if _lambdav is not None: assert len(_lambdav) == len(wavelengths)
                    else: _lambdav = wavelengths

        return _lambdav, _Zrelv, _logCv, _logpv, j0_dict, j1_dict

    # -----------------------------------------------------------------

    def luminosity_at(self, wavelength, unit="W/micron"):

        """
        This function ...
        :param wavelength:
        :param unit:
        :return:
        """

        return self.sed.luminosity_at(wavelength).to(unit)

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, fltr, unit="W/micron"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        #luminosity = filter.integrate(self.sed["Wavelength"], self.sed["Luminosity"])
        luminosity = fltr.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.luminosities(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * Unit("W/micron")

        # Return the luminosity in the desired unit
        return luminosity.to(unit)

# -----------------------------------------------------------------

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

# -----------------------------------------------------------------

def locate_clip(array, value):

    """
    This function ...
    :param array: actually a list
    :param value:
    :return:
    """

    #n = array.size
    n = len(array)
    if value < array[0]: return 0
    return locate_basic_impl(array, value, n-1)

# -----------------------------------------------------------------

def locate_basic_impl(xv, x, n):

    """
    This function ...
    :param xv:
    :param x:
    :param n:
    :return:
    """

    jl = -1
    ju = n

    while ju - jl > 1:

        jm = (ju + jl) >> 1
        if x < xv[jm]: ju = jm
        else: jl = jm

    return jl

# -----------------------------------------------------------------
