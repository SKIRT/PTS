#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.seds Contains several pre-defined SEDs.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from .sed import SED

# -----------------------------------------------------------------

seds_path = fs.join(introspection.pts_dat_dir("modeling"), "seds")

# -----------------------------------------------------------------

def load_example_mappings_sed():

    """
    This function ...
    :return:
    """

    # Determine the path to the SED file
    sed_path = fs.join(seds_path, "mapsed.dat")

    # Get the data
    wavelength_column, flux_column = np.loadtxt(sed_path, usecols=(0, 1), unpack=True)

    # Create an SED instance
    #sed = SED.initialize("W/m2")

    sed = SED.from_arrays(wavelength_column, flux_column, "micron", "W/m2", density=True)

    # Set the columns
    #sed.table = tables.new([wavelength_column, flux_column], ["Wavelength", "Flux"])
    #sed.table["Wavelength"].unit = "micron"
    #sed.table["Flux"].unit = "W/m2" # = lambda * F_Lambda !

    # Return the SED
    return sed

# -----------------------------------------------------------------

def load_example_bruzualcharlot_sed():

    """
    This function ...
    :return:
    """

    # Determine the path to the SED file
    sed_path = fs.join(seds_path, "bcsed.dat")

    # Get the data
    wavelength_column, flux_column = np.loadtxt(sed_path, usecols=(0, 1), unpack=True)

    # Create the SED instance
    #sed = SED()

    sed = SED.from_arrays(wavelength_column, flux_column, "micron", "W/m2", density=True)

    # Set the columns
    #sed.table = tables.new([wavelength_column, flux_column], ["Wavelength", "Flux"])
    #sed.table["Wavelength"].unit = "micron"
    #sed.table["Flux"].unit = "W/m2" # = lambda * F_lambda !

    # Return the SED
    return sed

# -----------------------------------------------------------------

def load_example_zubko_sed():

    """
    This function ...
    :return:
    """

    # Determine the path to the SED file
    sed_path = fs.join(seds_path, "zubkosed.dat")

    # Get the data
    wavelength_column, flux_column = np.loadtxt(sed_path, usecols=(0, 2), unpack=True)

    # Create the SED instance
    #sed = SED()

    sed = SED.from_arrays(wavelength_column, flux_column, "micron", "W/m2", density=True)

    # Set the columns
    #sed.table = tables.new([wavelength_column, flux_column], ["Wavelength", "Flux"])
    #sed.table["Wavelength"].unit = "micron"
    #sed.table["Flux"].unit = "W/m2" # = lambda * F_lambda

    # Return the SED
    return sed

# -----------------------------------------------------------------

def load_example_themis_sed():

    """
    This function ...
    :return:
    """

    raise NotImplementedError("Not yet implemented")

    # Determine the path to the SED file
    # ...

# -----------------------------------------------------------------
