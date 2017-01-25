#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.unit import parse_unit as u

# -----------------------------------------------------------------

# GENERAL:

# -----------------------------------------------------------------

# Types of parameters
possible_parameter_types = ["dimless", "mass", "grainsize", "length", "angle", "posangle", "luminosity",
                            "spectral luminosity density", "flux", "spectral flux density", "intensity",
                            "spectral intensity density", "surface brightness", "spectral surface brightness density",
                            "pressure"]

# -----------------------------------------------------------------

def is_photometric_density(parameter_type):

    """
    This function ...
    :param parameter_type:
    :return:
    """

    return parameter_type.startswith("spectral") and parameter_type.endswith("density")

# -----------------------------------------------------------------

def unit_parsing_type(parameter_type):

    """
    This function ...
    :param parameter_type:
    :return:
    """

    ptype = "photometric_density_unit" if is_photometric_density(parameter_type) else "unit"
    return ptype

# -----------------------------------------------------------------

# Default units for different parameter types
default_units = dict()
default_units["dimless"] = None
default_units["mass"] = u("Msun")
default_units["grainsize"] = u("micron")
default_units["length"] = u("pc")
default_units["angle"] = u("deg")
default_units["posangle"] = u("deg")
default_units["luminosity"] = u("Lsun") # bolometric luminosity
default_units["spectral luminosity density"] = u("W/micron", density=True) # spectral luminosity density
default_units["flux"] = u("W/m2", density=True) # bolometric flux
default_units["spectral flux density"] = u("Jy", density=True) # spectral flux density
default_units["intensity"] = u("W/sr") # bolometric intensity
default_units["spectral intensity density"] = u("W/sr/micron", density=True) # spectral intensity density
default_units["surface brightness"] = u("W/m2/sr") # bolometric surface brightness
default_units["spectral surface brightness density"] = u("W/m2/sr/micron", density=True) # spectral surface brightness density
default_units["pressure"] = u("K/m3")

# -----------------------------------------------------------------

possible_parameter_types_descriptions = dict()
possible_parameter_types_descriptions["dimless"] = "dimensionless quantity (no unit)"
possible_parameter_types_descriptions["mass"] = "mass (default unit: " + default_units["mass"] + ")"
possible_parameter_types_descriptions["grainsize"] = "grain size (default unit: " + default_units["grainsize"] + ")"
possible_parameter_types_descriptions["length"] = "physical length (default unit: " + default_units["length"] + ")"
possible_parameter_types_descriptions["angle"] = "angle (default unit: " + default_units["angle"] + ")"
possible_parameter_types_descriptions["posangle"] = "position angle (default unit: " + default_units["posangle"] + ")"
possible_parameter_types_descriptions["luminosity"] = "bolometric luminosity (default unit: " + default_units["luminosity"] + ")"
possible_parameter_types_descriptions["spectral luminosity density"] = "spectral luminosity (default unit: " + default_units["spectral luminosity density"] + ")"
possible_parameter_types_descriptions["flux"] = "bolometric flux (default unit: " + default_units["flux"] + ")"
possible_parameter_types_descriptions["spectral flux density"] = "spectral flux (default unit: " + default_units["spectral flux density"] + ")"
possible_parameter_types_descriptions["intensity"] = "bolometric intensity (default unit: " + default_units["intensity"] + ")"
possible_parameter_types_descriptions["spectral intensity density"] = "spectral intensity (default unit: " + default_units["spectral intensity density"] + ")"
possible_parameter_types_descriptions["surface brightness"] = "bolometric surface brightness (default unit: " + default_units["surface brightness"] + ")"
possible_parameter_types_descriptions["spectral surface brightness density"] = "spectral surface brightness (default unit: " + default_units["spectral surface brightness density"] + ")"
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
parsing_types_for_parameter_types["spectral luminosity density"] = "photometric_density_quantity"
parsing_types_for_parameter_types["flux"] = "photometric_quantity"
parsing_types_for_parameter_types["spectral flux density"] = "photometric_density_quantity"
parsing_types_for_parameter_types["intensity"] = "photometric_quantity"
parsing_types_for_parameter_types["spectral intensity density"] = "photometric_density_quantity"
parsing_types_for_parameter_types["surface brightness"] = "photometric_quantity"
parsing_types_for_parameter_types["spectral surface brightness density"] = "photometric_density_quantity"
parsing_types_for_parameter_types["pressure"] = "quantity"

# -----------------------------------------------------------------

# FOR GALAXY MODELING:

# -----------------------------------------------------------------

# Choices and descriptions of the different parameters
parameter_descriptions = dict()
parameter_descriptions["distance"] = "distance of the galaxy"
parameter_descriptions["ionizing_scaleheight"] = "scale height of the ionizing stellar component"
parameter_descriptions["sfr_compactness"] = "compactness parameter of the star formation regions"
parameter_descriptions["fuv_young"] = "FUV luminosity of the young stellar component"
parameter_descriptions["old_scaleheight"] = "scale height of the old stellar disk component"
parameter_descriptions["position_angle"] = "position angle of the galaxy"
parameter_descriptions["dust_mass"] = "dust mass"
parameter_descriptions["fuv_ionizing"] = "FUV luminosity of the ionizing stellar component"
parameter_descriptions["metallicity"] = "metallicity"
parameter_descriptions["young_scaleheight"] = "scale height of the young stellar component"
parameter_descriptions["sfr_covering"] = "covering factor of the star formation regions"
parameter_descriptions["dust_scaleheight"] = "scale height of the dust component"
parameter_descriptions["i1_old"] = "I1 luminosity of the old stellar component"
parameter_descriptions["sfr_pressure"] = "pressure on the star formation regions"
parameter_descriptions["inclination"] = "inclination of the galactic plane"

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
units["sfr_compactness"] = None
units["fuv_young"] = "W/micron"
units["old_scaleheight"] = "pc"
units["position_angle"] = "deg"
units["dust_mass"] = "Msun"
units["fuv_ionizing"] = "W/micron"
units["metallicity"] = None
units["young_scaleheight"] = "pc"
units["sfr_covering"] = None
units["dust_scaleheight"] = "pc"
units["i1_old"] = "W/micron"
units["sfr_pressure"] = "K/m3"
units["inclination"] = "deg"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add the required setting of the list of free parameters
definition.add_required("free_parameters", "string_list", "parameters to be used as free parameters during the fitting", choices=parameter_descriptions)

# -----------------------------------------------------------------
