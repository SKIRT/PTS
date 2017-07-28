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
from ....magic.core.frame import Frame
from ....magic.core.list import NamedFrameList

# -----------------------------------------------------------------

def make_tir_to_uv(tir, fuv, **kwargs):

    """
    This function ...
    :param tir: 
    :param fuv:
    :param kwargs:
    :return: 
    """

    # Conversions necessary? -> YES!

    ## Convert the FUV map from Lsun to W/m2
    #assert self.frames["GALEX FUV"].unit == "Lsun"
    ## Convert the TIR map from Lsun to W / m2
    #conversion_factor = 1.0
    # Conversion from Lsun to W
    #conversion_factor *= solar_luminosity.to("W").value
    # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
    #distance = self.galaxy_properties.distance
    #conversion_factor /= (4. * np.pi * distance ** 2).to("m2").value
    # FUV in W/M2
    #self.fuv_si = self.frames["GALEX FUV"] * conversion_factor
    #self.fuv_si.unit = "W/m2"

    ## Convert the TIR map from Lsun to W / m2
    #conversion_factor = 1.0
    # Conversion from Lsun to W
    #conversion_factor *= solar_luminosity.to("W").value
    # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
    #distance = self.galaxy_properties.distance
    #conversion_factor /= (4. * np.pi * distance ** 2).to("m2").value
    ## CONVERT AND SET NEW UNIT
    #self.tir_si = Frame(tir_map * conversion_factor)
    #self.tir_si.unit = "W/m2"

    # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

    ## FUV IN W/M2

    #print("FUV:", fuv.unit)
    #print("TIR:", tir.unit)

    # Create frame list
    #frames = NamedFrameList(fuv=fuv.copy(), tir=tir.copy())

    # FIRST CONVERT THEM BOTH (UNITS ARE IN FACT DIFFERENT, ONE IS DENSITY, OTHER IS NOT)
    fuv = fuv.copy()
    fuv.convert_to("W/m2", density=True, density_strict=True, **kwargs) # here it is a neutral density!
    
    tir = tir.copy()
    tir_mapData = tir.data.astype('float64') # Necessary for extreme conversion factors

    factor = tir.convert_to("W/m2", density=False, density_strict=True, **kwargs) # here it is bolometric!

    tir._data = tir_mapData * factor

    #frames.convert_to_same_unit("W/m2", density=True)

    ## TIR IN W/M2

    # Convert to same pixelscale and convolve to same resolution
    frames = NamedFrameList(fuv=fuv, tir=tir)
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
