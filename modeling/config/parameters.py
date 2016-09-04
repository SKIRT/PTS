#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Choices for free parameters
choices = ["distance", "ionizing_scaleheight", "sfr_compactness", "fuv_young", "old_scaleheight", "position_angle",
           "dust_mass", "fuv_ionizing", "metallicity", "young_scaleheight", "sfr_covering", "dust_scaleheight",
           "i1_old", "sfr_pressure", "inclination"]

# Descriptions of the different parameters
descriptions = dict()
descriptions["distance"] = "distance of the galaxy"
descriptions["ionizing_scaleheight"] = "scale height of the ionizing stellar component"
descriptions["sfr_compactness"] = "compactness parameter of the star formation regions"
descriptions["fuv_young"] = "FUV luminosity of the young stellar component"
descriptions["old_scaleheight"] = "scale height of the old stellar disk component"
descriptions["position_angle"] = "position angle of the galaxy"
descriptions["dust_mass"] = "dust mass"
descriptions["fuv_ionizing"] = "FUV luminosity of the ionizing stellar component"
descriptions["metallicity"] = "metallicity"
descriptions["young_scaleheight"] = "scale height of the young stellar component"
descriptions["sfr_covering"] = "covering factor of the star formation regions"
descriptions["dust_scaleheight"] = "scale height of the dust component"
descriptions["i1_old"] = "I1 luminosity of the old stellar component"
descriptions["sfr_pressure"] = "pressure on the star formation regions"
descriptions["inclination"] = "inclination of the galactic plane"

# Types and ranges of the different parameters
types_and_ranges = dict()
types_and_ranges["distance"] = ("quantity", None)
types_and_ranges["ionizing_scaleheight"] = ("quantity", None)
types_and_ranges["sfr_compactness"] = ("quantity", None)
types_and_ranges["fuv_young"] = ("real", "0.0 W/micron>4.e16 W/micron")
types_and_ranges["old_scaleheight"] = ("quantity", None)
types_and_ranges["position_angle"] = ("angle", None)
types_and_ranges["dust_mass"] = ("quantity", "0.5e7 Msun>3.e7 Msun")
types_and_ranges["fuv_ionizing"] = ("real", "0.0 W/micron>5.e16 W/micron")
types_and_ranges["metallicity"] = ("real", None)
types_and_ranges["young_scaleheight"] = ("quantity", None)
types_and_ranges["sfr_covering"] = ("real", None)
types_and_ranges["dust_scaleheight"] = ("quantity", None)
types_and_ranges["i1_old"] = ("real", None)
types_and_ranges["sfr_pressure"] = ("quantity", None)
types_and_ranges["inclination"] = ("angle", None)

# Default units of the different parameters
units = dict()
units["distance"] = "Mpc"
units["ionizing_scaleheight"] = "pc"
# sfr_compactness: no unit
units["fuv_young"] = "W/micron"
units["old_scaleheight"] = "pc"
units["position_angle"] = "deg"
units["dust_mass"] = "Msun"
units["fuv_ionizing"] = "W/micron"
# metallicity: no unit
units["young_scaleheight"] = "pc"
# sfr_covering: no unit
units["dust_scaleheight"] = "pc"
units["i1_old"] = "W/micron"
units["sfr_pressure"] = "K/m3"
units["inclination"] = "deg"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add the required setting of the list of free parameters
definition.add_required("free_parameters", "string_list", "the parameters to be used as free parameters during the fitting", choices=choices, choice_descriptions=descriptions)

# -----------------------------------------------------------------
