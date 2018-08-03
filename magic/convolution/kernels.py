#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.convolution.kernels Contains functions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.filter.filter import parse_filter
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# The path to the PTS kernels directory
if not fs.is_directory(introspection.pts_ext_dir): fs.create_directory(introspection.pts_ext_dir)
kernels_path = fs.join(introspection.pts_ext_dir, "kernels")
if not fs.is_directory(kernels_path): fs.create_directory(kernels_path)

# -----------------------------------------------------------------

# GALEX - SPIRE:
# Common-Resolution Convolution Kernels for Space- and Ground-Based Telescopes (G. Aniano et. al)

# Planck:
# DustPedia: Multiwavelength Photometry and Imagery of 875 Nearby Galaxies in 42 Ultraviolet–Microwave Bands (Clark et al.)
# - wavelength (micron)
# Planck 350 278
# Planck 550 290
# Planck 850 296
# Planck 1380 301
# Planck 2100 438
# Planck 3000 581
# Planck 4260 799
# Planck 6810 1630
# Planck 10600 1940

# FWHMS for different bands
fwhms = {"GALEX FUV": 4.48 * u("arcsec"),
         "GALEX NUV": 5.05 * u("arcsec"),
         "UVOT U": 2.3 * u("arcsec"), # *
         "UVOT B": 2.0 * u("arcsec"), # *
         "UVOT V": 1.6 * u("arcsec"), # *
         "UVOT UVW1": 1.7 * u("arcsec"), # *
         "UVOT UVM2": 2.0 * u("arcsec"), # *
         "UVOT UVW2": 2.3 * u("arcsec"), # *
         "Halpha": 2.0 * u("arcsec"),
         "IRAC I1": 1.90 * u("arcsec"),
         "IRAC I2": 1.81 * u("arcsec"),
         "IRAC I3": 2.11 * u("arcsec"),
         "IRAC I4": 2.82 * u("arcsec"),
         "WISE W1": 5.79 * u("arcsec"),
         "WISE W2": 6.37 * u("arcsec"),
         "WISE W3": 6.60 * u("arcsec"),
         "WISE W4": 11.89 * u("arcsec"),
         "MIPS 24mu": 6.43 * u("arcsec"),
         "MIPS 70mu": 18.74 * u("arcsec"),
         "MIPS 160mu": 38.78 * u("arcsec"),
         "Pacs blue": 5.67 * u("arcsec"),
         "Pacs green": 7.04 * u("arcsec"),
         "Pacs red": 11.18 * u("arcsec"),
         "SPIRE PSW": 18.15 * u("arcsec"),
         "SPIRE PMW": 24.88 * u("arcsec"),
         "SPIRE PLW": 36.09 * u("arcsec"),
         "Planck 350": 278. * u("arcsec"),
         "Planck 550": 290. * u("arcsec"),
         "Planck 850": 296. * u("arcsec"),
         "Planck 1380": 301. * u("arcsec"),
         "Planck 2100": 438. * u("arcsec"),
         "Planck 3000": 581. * u("arcsec"),
         "Planck 4260": 799. * u("arcsec"),
         "Planck 6810": 1630. * u("arcsec"),
         "Planck 10600": 1940. * u("arcsec"),
         "IRAS 12mu": 270 * u("arcsec"),
         "IRAS 25mu": 276 * u("arcsec"),
         "IRAS 60mu": 282 * u("arcsec"),
         "IRAS 100mu": 300 * u("arcsec")}

# -----------------------------------------------------------------

# REFERENCES:
# * from https://www.mssl.ucl.ac.uk/www_astro/uvot/uvot_instrument/filterwheel/filterwheel.html

# -----------------------------------------------------------------

variable_fwhms = ["SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z", "2MASS H", "2MASS J", "2MASS Ks", "UKIDSS J",
                  "UKIDSS H", "UKIDSS K", "Halpha", "Ha"]

# -----------------------------------------------------------------

def get_variable_filters():

    """
    This function ...
    :return:
    """

    return [parse_filter(filter_name) for filter_name in variable_fwhms]

# -----------------------------------------------------------------

def get_filters_and_fwhms_lists():

    """
    This function ...
    :return:
    """

    filters = []
    fwhm_list = []

    # Loop over the filters
    for filter_name, fwhm in fwhms.iteritems():

        # Parse into filter
        fltr = parse_filter(filter_name)

        # add to lists
        filters.append(fltr)
        fwhm_list.append(fwhm)

    # Return the lists
    return filters, fwhm_list

# -----------------------------------------------------------------

def get_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # Determine the aniano name for the PSF
    #psf_name = aniano_names[str(fltr)]
        
    #if has_average_variable_fwhm(fltr):return(get_average_variable_fwhm(fltr))

    if has_variable_fwhm(fltr): raise ValueError("The specified filter has a variable FWHM")

    # Get the FWHM of the PSF
    #if "Gauss" in psf_name or "Moffet" in psf_name:
    #    fwhm = float(psf_name.split("_")[1]) * parse_unit("arcsec")
    #elif psf_name in fwhms:
    #    fwhm = fwhms[psf_name]
    #else: fwhm = None

    # Return the FWHM
    #return fwhm

    filters, fwhm_list = get_filters_and_fwhms_lists()
    #print([str(filt) for filt in filters])
    #print(str(fltr))
    #print([str(filter) for filter in filters])
    #print(fltr)
    index = filters.index(fltr)
    return fwhm_list[index]

# -----------------------------------------------------------------

def has_variable_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    return fltr in get_variable_filters()

# -----------------------------------------------------------------

def has_average_variable_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...dustpedia.core.properties import has_fwhm as has_dustpedia_fwhm
    return has_dustpedia_fwhm(fltr)

# -----------------------------------------------------------------

def get_average_variable_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...dustpedia.core.properties import get_fwhm as get_dustpedia_fwhm
    return get_dustpedia_fwhm(fltr)

# -----------------------------------------------------------------

class Kernels(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @abstractmethod
    def get_kernel(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def get_psf(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        pass

# -----------------------------------------------------------------
