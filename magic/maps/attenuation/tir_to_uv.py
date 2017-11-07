#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.dust.tir_to_uv Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....magic.core.list import NamedFrameList
from ....core.tools.stringify import tostr
from ....core.basics.log import log
from ....core.units.unit import PhotometricUnit
from ...core.image import Image
from ...core.frame import Frame

# -----------------------------------------------------------------

def make_tir_to_uv(tir, fuv, **kwargs):

    """
    This function ...
    :param tir: 
    :param fuv:
    :param kwargs:
    :return: 
    """

    # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

    # Debugging
    log.debug("Creating a TIR to UV map ...")

    # Get the distance
    if "distance" in kwargs: distance = kwargs.pop("distance")
    else:
        distance = tir.distance
        if distance is None: distance = fuv.distance
        if distance is None: raise ValueError("Distance could not be determined")

    ## FUV IN W/M2

    #print("FUV:", fuv.unit)
    #print("TIR:", tir.unit)

    # Create frame list
    #frames = NamedFrameList(fuv=fuv.copy(), tir=tir.copy())

    # FIRST CONVERT THEM BOTH (UNITS ARE IN FACT DIFFERENT, ONE IS DENSITY, OTHER IS NOT)
    fuv = fuv.copy()
    fuv.convert_to("W/m2", density=True, density_strict=True, brightness=False, brightness_strict=True, distance=distance) # here it is a neutral density!

    # UNIT CONVERSION OF TIR MAP
    if isinstance(tir, Image): tir_frame = tir.primary
    elif isinstance(tir, Frame): tir_frame = tir.copy()
    else: raise ValueError("Something went wrong")

    tir_map_data = tir_frame.data.astype('float64') # Necessary for extreme conversion factors
    #factor = tir.convert_to("W/m2", density=False, density_strict=True, **kwargs) # here it is bolometric!
    factor = tir_frame.unit.conversion_factor("W/m2", density=False, density_strict=True, brightness=False, brightness_strict=True, distance=distance, pixelscale=tir.pixelscale) # HERE IT IS BOLOMETRIC, AND NO FLUX!!!!
    log.debug("Conversion factor for the TIR map from " + tostr(tir.unit, add_physical_type=True) + " to W/m2 is " + str(factor))

    tir_frame._data = tir_map_data * factor
    tir_frame.unit = PhotometricUnit("W/m2", density=False, density_strict=True, brightness=False, brightness_strict=True)

    #frames.convert_to_same_unit("W/m2", density=True)

    ## TIR IN W/M2

    # Convert to same pixelscale and convolve to same resolution
    frames = NamedFrameList(fuv=fuv, tir=tir_frame)
    frames.convolve_and_rebin()
    
    # CALCULATE TIR TO FUV RATIO
    #print(frames.names)
    
    # The ratio of TIR and FUV
    tir_to_fuv = frames["tir"] / frames["fuv"]

    # Check unit and WCS
    tir_to_fuv.unit = None
    assert tir_to_fuv.wcs is not None

    # Check filter
    tir_to_fuv.filter = None

    # Set other properties than unit and wcs
    tir_to_fuv.pixelscale = frames.pixelscale
    tir_to_fuv.psf_filter = frames.psf_filter
    tir_to_fuv.fwhm = frames.fwhm
    tir_to_fuv.distance = frames.distance

    # Return a frame
    return tir_to_fuv

# -----------------------------------------------------------------
