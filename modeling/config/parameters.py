#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# GENERAL:

# -----------------------------------------------------------------

# Types of parameters
possible_parameter_types = ["dimless", "mass", "grainsize", "length", "angle", "posangle", "luminosity", "pressure"]

# -----------------------------------------------------------------

# Default units for different parameter types
default_units = dict()
default_units["dimless"] = ""
default_units["mass"] = "Msun"
default_units["grainsize"] = "micron"
default_units["length"] = "pc"
default_units["angle"] = "deg"
default_units["posangle"] = "deg"
default_units["luminosity"] = "W/micron"
default_units["pressure"] = "K/m3"

# -----------------------------------------------------------------

possible_parameter_types_descriptions = dict()
possible_parameter_types_descriptions["dimless"] = "dimensionless quantity (no unit)"
possible_parameter_types_descriptions["mass"] = "mass (default unit: " + default_units["mass"] + ")"
possible_parameter_types_descriptions["grainsize"] = "grain size (default unit: " + default_units["grainsize"] + ")"
possible_parameter_types_descriptions["length"] = "physical length (default unit: " + default_units["length"] + ")"
possible_parameter_types_descriptions["angle"] = "angle (default unit: " + default_units["angle"] + ")"
possible_parameter_types_descriptions["posangle"] = "position angle (default unit: " + default_units["posangle"] + ")"
possible_parameter_types_descriptions["luminosity"] = "(spectral) luminosity (default unit: " + default_units["luminosity"] + ")"
possible_parameter_types_descriptions["pressure"] = "pressure (default unit: " + default_units["pressure"] + ")"

# -----------------------------------------------------------------

parsing_types_for_parameter_types = dict()
parsing_types_for_parameter_types["dimless"] = "real"
parsing_types_for_parameter_types["mass"] = "quantity"
parsing_types_for_parameter_types["grainsize"] = "quantity"
parsing_types_for_parameter_types["length"] = "quantity"
parsing_types_for_parameter_types["angle"] = "angle"
parsing_types_for_parameter_types["posangle"] = "angle"
parsing_types_for_parameter_types["luminosity"] = "photometric_quantity"
parsing_types_for_parameter_types["pressure"] = "quantity"

# -----------------------------------------------------------------

# FOR GALAXY MODELING:

# -----------------------------------------------------------------

# Choices and descriptions of the different parameters
choices = dict()
choices["distance"] = "distance of the galaxy"
choices["ionizing_scaleheight"] = "scale height of the ionizing stellar component"
choices["sfr_compactness"] = "compactness parameter of the star formation regions"
choices["fuv_young"] = "FUV luminosity of the young stellar component"
choices["old_scaleheight"] = "scale height of the old stellar disk component"
choices["position_angle"] = "position angle of the galaxy"
choices["dust_mass"] = "dust mass"
choices["fuv_ionizing"] = "FUV luminosity of the ionizing stellar component"
choices["metallicity"] = "metallicity"
choices["young_scaleheight"] = "scale height of the young stellar component"
choices["sfr_covering"] = "covering factor of the star formation regions"
choices["dust_scaleheight"] = "scale height of the dust component"
choices["i1_old"] = "I1 luminosity of the old stellar component"
choices["sfr_pressure"] = "pressure on the star formation regions"
choices["inclination"] = "inclination of the galactic plane"

# -----------------------------------------------------------------

# Types and ranges of the different parameters
#types_and_ranges = dict()
#types_and_ranges["distance"] = ("quantity", None)
#types_and_ranges["ionizing_scaleheight"] = ("quantity", None)
#types_and_ranges["sfr_compactness"] = ("quantity", None)
#types_and_ranges["fuv_young"] = ("quantity", "0.0 W/micron>1e37 W/micron")
#types_and_ranges["old_scaleheight"] = ("quantity", None)
#types_and_ranges["position_angle"] = ("angle", None)
#types_and_ranges["dust_mass"] = ("quantity", "0.5e7 Msun>3.e7 Msun")
#types_and_ranges["fuv_ionizing"] = ("quantity", "0.0 W/micron>1e34 W/micron")
#types_and_ranges["metallicity"] = ("real", None)
#types_and_ranges["young_scaleheight"] = ("quantity", None)
#types_and_ranges["sfr_covering"] = ("real", None)
#types_and_ranges["dust_scaleheight"] = ("quantity", None)
#types_and_ranges["i1_old"] = ("real", None)
#types_and_ranges["sfr_pressure"] = ("quantity", None)
#types_and_ranges["inclination"] = ("angle", None)

# -----------------------------------------------------------------

types = dict()
types["distance"] = "length"
types["ionizing_scaleheight"] = "length"
types["sfr_compactness"] = "dimless"
types["fuv_young"] = "luminosity"
types["old_scaleheight"] = "length"
types["position_angle"] = "angle"
types["dust_mass"] = "mass"
types["fuv_ionizing"] = "luminosity"
types["metallicity"] = "dimless"
types["young_scaleheight"] = "length"
types["sfr_covering"] = "dimless"
types["dust_scaleheight"] = "length"
types["i1_old"] = "luminosity"
types["sfr_pressure"] = "pressure"
types["inclination"] = "angle"

# -----------------------------------------------------------------

default_ranges = dict()
default_ranges["fuv_young"] = "0.0 W/micron>1e37 W/micron"
default_ranges["dust_mass"] = "0.5e7 Msun>3.e7 Msun"
default_ranges["fuv_ionizing"] = "0.0 W/micron>1e34 W/micron"

# -----------------------------------------------------------------

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
definition.add_required("free_parameters", "string_list", "parameters to be used as free parameters during the fitting", choices=choices)

# -----------------------------------------------------------------
