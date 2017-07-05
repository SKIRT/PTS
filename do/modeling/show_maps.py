#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_maps Show the maps

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.magic.plot.imagegrid import StandardImageGridPlotter
from pts.magic.core.frame import Frame
from pts.modeling.core.environment import map_sub_names
from pts.modeling.maps.component import get_maps_sub_name
from pts.modeling.component.component import load_modeling_history

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_required("sub_name", "string", "maps sub name for plotting", map_sub_names)
definition.add_positional_optional("method", "string", "method for plotting")

# Get configuration
config = parse_arguments("plot_maps", definition)

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Load the modeling history
history = load_modeling_history(modeling_path)

# -----------------------------------------------------------------

# Create plotter
plotter = StandardImageGridPlotter()

# Loop over the maps
maps = get_maps_sub_name(environment, history, config.sub_name, method=config.method)
for name in maps:

    # Add the map
    plotter.add_image(maps[name], name)

# -----------------------------------------------------------------

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
