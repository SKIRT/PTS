#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.units.parsing import parse_unit as u
from pts.core.units.stringify import represent_unit as ru

# -----------------------------------------------------------------

# GENERAL:

# -----------------------------------------------------------------

# Define parameter type names
dimless_name = "dimless"
mass_name = "mass"
grainsize_name = "grainsize"
length_name = "length"
angle_name = "angle"
posangle_name = "posangle"
luminosity_name = "luminosity"
spectral_luminosity_density_name = "spectral luminosity density"
flux_name = "flux"
spectral_flux_density_name = "spectral flux density"
intensity_name = "intensity"
spectral_intensity_density_name = "spectral intensity density"
surface_brightness_name = "surface brightness"
spectral_surface_brightness_density_name = "spectral surface brightness density"
pressure_name = "pressure"

# -----------------------------------------------------------------

# Types of parameters
possible_parameter_types = [dimless_name, mass_name, grainsize_name, length_name, angle_name, posangle_name, luminosity_name,
                            spectral_luminosity_density_name, flux_name, spectral_flux_density_name, intensity_name,
                            spectral_intensity_density_name, surface_brightness_name, spectral_surface_brightness_density_name,
                            pressure_name]

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

    if parsing_types_for_parameter_types[parameter_type] == "photometric_quantity":
        if is_photometric_density(parameter_type): ptype = "photometric_density_unit"
        else: ptype = "photometric_unit"
    else: ptype = "unit"

    return ptype

# -----------------------------------------------------------------

# Default units for different parameter types
default_units = dict()
default_units[dimless_name] = None
default_units[mass_name] = u("Msun")
default_units[grainsize_name] = u("micron")
default_units[length_name] = u("pc")
default_units[angle_name] = u("deg")
default_units[posangle_name] = u("deg")
default_units[luminosity_name] = u("Lsun") # bolometric luminosity
default_units[spectral_luminosity_density_name] = u("W/micron", density=True) # spectral luminosity density
default_units[flux_name] = u("W/m2", density=True) # bolometric flux
default_units[spectral_flux_density_name] = u("Jy", density=True) # spectral flux density
default_units[intensity_name] = u("W/sr") # bolometric intensity
default_units[spectral_intensity_density_name] = u("W/sr/micron", density=True) # spectral intensity density
default_units[surface_brightness_name] = u("W/m2/sr") # bolometric surface brightness
default_units[spectral_surface_brightness_density_name] = u("W/m2/sr/micron", density=True) # spectral surface brightness density
default_units[pressure_name] = u("K/m3")

# -----------------------------------------------------------------

possible_parameter_types_descriptions = dict()
possible_parameter_types_descriptions[dimless_name] = "dimensionless quantity (no unit)"
possible_parameter_types_descriptions[mass_name] = "mass (default unit: " + ru(default_units["mass"]) + ")"
possible_parameter_types_descriptions[grainsize_name] = "grain size (default unit: " + ru(default_units["grainsize"]) + ")"
possible_parameter_types_descriptions[length_name] = "physical length (default unit: " + ru(default_units["length"]) + ")"
possible_parameter_types_descriptions[angle_name] = "angle (default unit: " + ru(default_units["angle"]) + ")"
possible_parameter_types_descriptions[posangle_name] = "position angle (default unit: " + ru(default_units["posangle"]) + ")"
possible_parameter_types_descriptions[luminosity_name] = "bolometric luminosity (default unit: " + ru(default_units["luminosity"]) + ")"
possible_parameter_types_descriptions[spectral_luminosity_density_name] = "spectral luminosity (default unit: " + ru(default_units["spectral luminosity density"]) + ")"
possible_parameter_types_descriptions[flux_name] = "bolometric flux (default unit: " + ru(default_units["flux"]) + ")"
possible_parameter_types_descriptions[spectral_flux_density_name] = "spectral flux (default unit: " + ru(default_units["spectral flux density"]) + ")"
possible_parameter_types_descriptions[intensity_name] = "bolometric intensity (default unit: " + ru(default_units["intensity"]) + ")"
possible_parameter_types_descriptions[spectral_intensity_density_name] = "spectral intensity (default unit: " + ru(default_units["spectral intensity density"]) + ")"
possible_parameter_types_descriptions[surface_brightness_name] = "bolometric surface brightness (default unit: " + ru(default_units["surface brightness"]) + ")"
possible_parameter_types_descriptions[spectral_surface_brightness_density_name] = "spectral surface brightness (default unit: " + ru(default_units["spectral surface brightness density"]) + ")"
possible_parameter_types_descriptions[pressure_name] = "pressure (default unit: " + ru(default_units["pressure"]) + ")"

# -----------------------------------------------------------------

parsing_types_for_parameter_types = dict()
parsing_types_for_parameter_types[dimless_name] = "real"
parsing_types_for_parameter_types[mass_name] = "mass_quantity"
parsing_types_for_parameter_types[grainsize_name] = "length_quantity"
parsing_types_for_parameter_types[length_name] = "length_quantity"
parsing_types_for_parameter_types[angle_name] = "angle"
parsing_types_for_parameter_types[posangle_name] = "angle"
parsing_types_for_parameter_types[luminosity_name] = "photometric_quantity"
parsing_types_for_parameter_types[spectral_luminosity_density_name] = "photometric_density_quantity"
parsing_types_for_parameter_types[flux_name] = "photometric_quantity"
parsing_types_for_parameter_types[spectral_flux_density_name] = "photometric_density_quantity"
parsing_types_for_parameter_types[intensity_name] = "photometric_quantity"
parsing_types_for_parameter_types[spectral_intensity_density_name] = "photometric_density_quantity"
parsing_types_for_parameter_types[surface_brightness_name] = "photometric_quantity"
parsing_types_for_parameter_types[spectral_surface_brightness_density_name] = "photometric_density_quantity"
parsing_types_for_parameter_types[pressure_name] = "quantity"

# -----------------------------------------------------------------

# FOR GALAXY MODELING:

# -----------------------------------------------------------------

# Define names
distance_name = "distance"
ionizing_scaleheight_name = "ionizing_scaleheight"
sfr_compactness_name = "sfr_compactness"
fuv_young_name = "fuv_young"
old_scaleheight_name = "old_scaleheight"
position_angle_name = "position_angle"
dust_mass_name = "dust_mass"
fuv_ionizing_name = "fuv_ionizing"
metallicity_name = "metallicity"
young_scaleheight_name = "young_scaleheight"
sfr_covering_name = "sfr_covering"
dust_scaleheight_name = "dust_scaleheight"
i1_old_name = "i1_old"
sfr_pressure_name = "sfr_pressure"
inclination_name = "inclination"

# -----------------------------------------------------------------

# Set modeling parameter labels
modeling_parameter_labels = [distance_name, ionizing_scaleheight_name, sfr_compactness_name, fuv_young_name,
                             old_scaleheight_name, position_angle_name, dust_mass_name, fuv_ionizing_name,
                             metallicity_name, young_scaleheight_name, sfr_covering_name, dust_scaleheight_name,
                             i1_old_name, sfr_pressure_name, inclination_name]

# -----------------------------------------------------------------

# Choices and descriptions of the different parameters
parameter_descriptions = dict()
parameter_descriptions[distance_name] = "distance of the galaxy"
parameter_descriptions[ionizing_scaleheight_name] = "scale height of the ionizing stellar component"
parameter_descriptions[sfr_compactness_name] = "compactness parameter of the star formation regions"
parameter_descriptions[fuv_young_name] = "FUV luminosity of the young stellar component"
parameter_descriptions[old_scaleheight_name] = "scale height of the old stellar disk component"
parameter_descriptions[position_angle_name] = "position angle of the galaxy"
parameter_descriptions[dust_mass_name] = "dust mass"
parameter_descriptions[fuv_ionizing_name] = "FUV luminosity of the ionizing stellar component"
parameter_descriptions[metallicity_name] = "metallicity"
parameter_descriptions[young_scaleheight_name] = "scale height of the young stellar component"
parameter_descriptions[sfr_covering_name] = "covering factor of the star formation regions"
parameter_descriptions[dust_scaleheight_name] = "scale height of the dust component"
parameter_descriptions[i1_old_name] = "I1 luminosity of the old stellar component"
parameter_descriptions[sfr_pressure_name] = "pressure on the star formation regions"
parameter_descriptions[inclination_name] = "inclination of the galactic plane"

# -----------------------------------------------------------------

# Short descriptions
parameter_descriptions_short = dict()
parameter_descriptions_short[distance_name] = "distance"
parameter_descriptions_short[ionizing_scaleheight_name] = "ionizing stars scaleheight"
parameter_descriptions_short[sfr_compactness_name] = "SFR compactness"
parameter_descriptions_short[fuv_young_name] = "young stars FUV luminosity"
parameter_descriptions_short[old_scaleheight_name] = "old stars scaleheight"
parameter_descriptions_short[position_angle_name] = "position angle"
parameter_descriptions_short[dust_mass_name] = "dust mass"
parameter_descriptions_short[fuv_ionizing_name] = "ionizing stars FUV luminosity"
parameter_descriptions_short[metallicity_name] = "metallicity"
parameter_descriptions_short[young_scaleheight_name] = "young stars scaleheight"
parameter_descriptions_short[sfr_covering_name] = "SFR covering factor"
parameter_descriptions_short[dust_scaleheight_name] = "dust scaleheight"
parameter_descriptions_short[i1_old_name] = "old stars I1 luminosity"
parameter_descriptions_short[sfr_pressure_name] = "SFR pressure"
parameter_descriptions_short[inclination_name] = "inclination"

# -----------------------------------------------------------------

types = dict()
types[distance_name] = "length"
types[ionizing_scaleheight_name] = "length"
types[sfr_compactness_name] = "dimless"
types[fuv_young_name] = "luminosity"
types[old_scaleheight_name] = "length"
types[position_angle_name] = "angle"
types[dust_mass_name] = "mass"
types[fuv_ionizing_name] = "luminosity"
types[metallicity_name] = "dimless"
types[young_scaleheight_name] = "length"
types[sfr_covering_name] = "dimless"
types[dust_scaleheight_name] = "length"
types[i1_old_name] = "luminosity"
types[sfr_pressure_name] = "pressure"
types[inclination_name] = "angle"

# -----------------------------------------------------------------

# Define default ranges
default_ranges = dict()
default_ranges[fuv_young_name] = "1e34 W/micron>1e38 W/micron"
default_ranges[dust_mass_name] = "1e6 Msun>1e8 Msun"
default_ranges[fuv_ionizing_name] = "1e34 W/micron>1e38 W/micron"

# -----------------------------------------------------------------

# Default units of the different parameters
units = dict()
units[distance_name] = "Mpc"
units[ionizing_scaleheight_name] = "pc"
units[sfr_compactness_name] = None
units[fuv_young_name] = "W/micron"
units[old_scaleheight_name] = "pc"
units[position_angle_name] = "deg"
units[dust_mass_name] = "Msun"
units[fuv_ionizing_name] = "W/micron"
units[metallicity_name] = None
units[young_scaleheight_name] = "pc"
units[sfr_covering_name] = None
units[dust_scaleheight_name] = "pc"
units[i1_old_name] = "W/micron"
units[sfr_pressure_name] = "K/m3"
units[inclination_name] = "deg"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add the required setting of the list of free parameters
definition.add_required("free_parameters", "string_list", "parameters to be used as free parameters during the fitting", choices=parameter_descriptions)

# -----------------------------------------------------------------
