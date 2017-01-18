#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.model_info Show information about the radiative transfer model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import filesystem as fs
from pts.modeling.component.component import load_modeling_configuration, load_fitting_configuration, load_modeling_history
from pts.modeling.component.galaxy import load_preparation_statistics
from pts.magic.core.frame import Frame

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("model_info", "Show information about the radiative transfer model")
config = setter.run(definition)

# -----------------------------------------------------------------

# Get the modeling path
modeling_path = config.path

# -----------------------------------------------------------------

# Get modeling and fitting config
modeling_config = load_modeling_configuration(modeling_path)
fitting_config = load_fitting_configuration(modeling_path)

# Get modeling history
history = load_modeling_history(modeling_path)

# Get preparation statistics
statistics = load_preparation_statistics(modeling_path)

# Maps paths
maps_path = fs.join(modeling_path, "maps")
old_stars_path = fs.join(maps_path, "old_stars.fits")
young_stars_path = fs.join(maps_path, "young_stars.fits")
ionizing_stars_path = fs.join(maps_path, "ionizing_stars.fits")
dust_path = fs.join(maps_path, "dust.fits")

# Open the old stars map
old_stars = Frame.from_file(old_stars_path)

# -----------------------------------------------------------------

convolution_filter = statistics.convolution_filter
rebinning_filter = statistics.rebinning_filter

# -----------------------------------------------------------------

print("")
print("Galaxy name: " + modeling_config.name + " (" + modeling_config.ngc_name + ")")
print("Modeling method: " + modeling_config.method)
print("Model pixelscale: " + str(old_stars.average_pixelscale) + " (" + str(rebinning_filter) + ")")
print("Model resolution (FWHM): " + str(old_stars.fwhm) + " (" + str(convolution_filter) + ")")
print("Reference fluxes (for the SED fitting): ")
print("")
for filter_name in fitting_config.filters: print(" - " + filter_name)
print("")

# -----------------------------------------------------------------
