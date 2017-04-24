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

# -----------------------------------------------------------------

def make_tir_to_uv(tir, fuv):

    """
    This function ...
    :param tir: 
    :param fuv: 
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



    ## TIR IN W/M2

    # CALCULATE TIR TO FUV RATIO

    # The ratio of TIR and FUV
    tir_to_fuv = tir_si / fuv_si
    #log_tir_to_fuv = Frame(np.log10(self.tir_to_fuv)

    # Return a frame
    return Frame(tir_to_fuv, wcs=uv_si.wcs)

# -----------------------------------------------------------------
