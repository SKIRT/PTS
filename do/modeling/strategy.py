#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.strategy Determine the modelling strategy for a certain galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.dustpedia.core.sample import DustPediaSample
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
definition.add_required("galaxy", "string", "galaxy name")

# Get configuration
setter = ArgumentConfigurationSetter("strategy")
config = setter.run(definition)

# -----------------------------------------------------------------

sample = DustPediaSample()

# Database
username, password = get_account()
database = DustPediaDatabase()
database.login(username, password)

# Get the DustPedia name
name = sample.get_name(config.galaxy)

# Get info
info = database.get_galaxy_info(name)

# FIRST IMAGES AND THAN THIS FAILS FOR SOME REASON!!
path = fs.join(introspection.pts_temp_dir, time.unique_name("cutouts_" + name) + ".png")
database.download_photometry_cutouts(name, path)

# Get image names
image_names = database.get_image_names(name, error_maps=False)

print("")
print(fmt.red + fmt.underlined + "General info:" + fmt.reset)
print("")

print(" - position: " + str(sample.get_position(name)))
print(" - D25: " + str(sample.get_d25(name)))
print(" - R25: " + str(sample.get_r25(name)))
print("")
print(" - stage: " + str(info["Hubble Stage"][0]))
print(" - V: " + str(info["V"][0]))
print(" - inclination: " + str(info["Inclination"][0]))

print("")
print(fmt.red + fmt.underlined + "Images:" + fmt.reset)
print("")

for name in image_names: print(" - " + name)
print("")

#urls = database.get_image_urls(name, error_maps=False)
#print(urls)

#database.reset(username, password)

# Open the cutouts file
fs.open_file(path)

# -----------------------------------------------------------------
