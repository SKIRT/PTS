#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_galaxies Plot the positions of the galaxies in the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.services.s4g import get_galaxy_names, has_galaxy
from pts.core.basics.map import Map
from pts.core.basics.log import log
from pts.core.basics.table import SmartTable

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Maximum number of galaxies
definition.add_optional("ngalaxies", "positive_integer", "max number of galaxies")

# Get configuration
config = parse_arguments("plot_galaxies", definition)

# -----------------------------------------------------------------

# Create the database instance
database = DustPediaDatabase()
#username, password = get_account()
#database.login(username, password)

# -----------------------------------------------------------------

sample = DustPediaSample()
galaxy_names = sample.get_names()

# -----------------------------------------------------------------

s4g_names = get_galaxy_names()

# -----------------------------------------------------------------

galaxies = []

# -----------------------------------------------------------------

# Loop over the names
for galaxy_name in galaxy_names:

    # Get info
    info = database.get_galaxy_info(galaxy_name)

    # Has S4G decomposition
    has_s4g = galaxy_name in s4g_names

    # Initialize the parameters
    sfr = None
    sfr_error = None
    dust_mass = None
    dust_mass_error = None
    dust_luminosity = None
    dust_luminosity_error = None
    dust_temperature = None
    dust_temperature_error = None
    stellar_mass = None
    stellar_mass_error = None
    stellar_luminosity = None
    stellar_luminosity_error = None

    # There are results from CIGALE fits
    if database.has_cigale_parameters(galaxy_name):

        # Get the parameters
        parameters = database.get_cigale_parameters(galaxy_name)

        # Dust parameters
        dust_mass = parameters["dust_mass"]
        dust_mass_error = parameters["dust_mass_error"]
        dust_luminosity = parameters["dust_luminosity"]
        dust_luminosity_error = parameters["dust_luminosity_error"]
        fuv_attenuation = parameters["fuv_attenuation"]
        fuv_attenuation_error = parameters["fuv_attenuation_error"]

        # Stellar parameters
        sfr = parameters["sfr"]
        sfr_error = parameters["sfr_error"]
        stellar_mass = parameters["stellar_mass"]
        stellar_mass_error = parameters["stellar_mass_error"]
        stellar_luminosity = parameters["stellar_luminosity"]
        stellar_luminosity_error = parameters["stellar_luminosity_error"]

    # There are results from black body fits
    elif database.has_dust_black_body_table_parameters(galaxy_name):

        # Get parameters
        parameters = database.get_dust_black_body_table_parameters(galaxy_name)

        # Dust parameters
        dust_mass = parameters["dust_mass"]
        dust_mass_error = parameters["dust_mass_error"]
        dust_luminosity = parameters["dust_luminosity"]
        dust_luminosity_error = parameters["dust_luminosity_error"]
        dust_temperature = parameters["dust_temperature"]
        dust_temperature_error = parameters["dust_temperature_error"]

    else: log.warning("No model parameters available for galaxy '" + galaxy_name + "'")

    # Check presence of data
    has_galex = database.has_galex(galaxy_name)
    has_sdss = database.has_sdss(galaxy_name)
    has_2mass = database.has_2mass(galaxy_name)
    has_spitzer = database.has_spitzer(galaxy_name)
    has_wise = database.has_wise(galaxy_name)
    has_pacs = database.has_pacs(galaxy_name)
    has_spire = database.has_spire(galaxy_name)
    has_planck = database.has_planck(galaxy_name)

    # Set galaxy properties
    properties = Map()

    # Set basic info
    properties.name = info.name
    properties.ra = info.position.ra
    properties.dec = info.position.dec
    properties.stage = info.stage
    properties.type = info.type
    properties.velocity = info.velocity
    properties.d25 = info.d25
    properties.inclination = info.inclination

    # Set model properties
    properties.sfr = sfr
    properties.sfr_error = sfr_error
    properties.dust_mass = dust_mass
    properties.dust_mass_error = dust_mass_error
    properties.dust_luminosity = dust_luminosity
    properties.dust_luminosity_error = dust_luminosity_error
    properties.dust_temperature = dust_temperature
    properties.dust_temperature_error = dust_temperature_error
    properties.stellar_mass = stellar_mass
    properties.stellar_mass_error = stellar_mass_error
    properties.stellar_luminosity = stellar_luminosity
    properties.stellar_luminosity_error = stellar_luminosity_error

    # Set flags
    properties.has_galex = has_galex
    properties.has_sdss = has_sdss
    properties.has_2mass = has_2mass
    properties.has_spitzer = has_spitzer
    properties.has_wise = has_wise
    properties.has_pacs = has_pacs
    properties.has_spire = has_spire
    properties.has_planck = has_planck

    # Add properties
    galaxies.append(properties)

    # Limit of number of galaxies
    if config.ngalaxies is not None and len(galaxies) == config.ngalaxies: break

# -----------------------------------------------------------------

# Create table
table = SmartTable.from_dictionaries(*galaxies, first="name", ignore_none=True)

# -----------------------------------------------------------------

# Determine the path
path = "galaxies.dat"

# Save the table
#table.saveto(path)

# -----------------------------------------------------------------

import numpy as np

#print(table["stellar_luminosity"].mask)

# for name in table.column_names:
#     mask = table.mask[name]
#     for index in range(len(table)):
#         masked = mask[index]
#         if masked: table[name][index] = None #np.nan #np.ma.core.MaskedConstant

#print(table)
#exit()

filename = "galaxies"
table.saveto_pts(filename)

#table.saveto_csv(filename)
#table.saveto_ecsv(filename)
#table.saveto_html(filename)
#table.saveto_latex(filename)
#table.saveto_votable(filename)
#table.saveto_ascii(filename)

# -----------------------------------------------------------------
