#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plot_curves

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.scatter import Scatter2D
from pts.core.basics.plot import MPLFigure
from pts.core.tools import filesystem as fs
from pts.magic.tools.plotting import plot_scatters, plot_densities, plot_scatters_density

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("files", "filepath_list", "list of the data file paths")
definition.add_flag("xlog", "use log scale for x axis")
definition.add_flag("ylog", "use log scale for y axis")
definition.add_flag("loglog", "use log-log scale")
definition.add_optional("xlimits", "real_pair_or_quantity_pair", "x interval")
definition.add_optional("ylimits", "real_pair_or_quantity_pair", "y interval")
definition.add_optional("path", "new_path", "path for the plot file")

# Special options
definition.add_flag("density", "plot the density of points")
definition.add_flag("points", "plot the points (as opposed to just the density)", True)
definition.add_optional("density_method", "string", "density method", default="hexbin", choices=["contours", "seaborn", "kde", "hexbin"])

# Labels
definition.add_optional("xlabel", "string", "x axis label")
definition.add_optional("ylabel", "string", "y axis label")

# Parse
config = parse_arguments("plot_scatters", definition)

# -----------------------------------------------------------------

# Set options
if config.loglog: config.xlog = config.ylog = True

# -----------------------------------------------------------------

# Load the curves
scatters = OrderedDict()
for filepath in config.files:
    name = fs.strip_extension(fs.name(filepath))
    scatters[name] = Scatter2D.from_file(filepath) #Scatter.from_data_file(filepath)

# -----------------------------------------------------------------

# Check
if not config.points and not config.density: raise ValueError("Must plot points and/or density")

# -----------------------------------------------------------------

# Scatter points (with or without density)
if config.points: output = plot_scatters(scatters, xlimits=config.xlimits, ylimits=config.ylimits, xlog=config.xlog, ylog=config.ylog, density=config.density, size=12, show=False)

# Density field
else: output = plot_densities(scatters, xlimits=config.xlimits, ylimits=config.ylimits, xlog=config.xlog, ylog=config.ylog, method=config.density_method, show=False)

# Set ticks in nice way
output.plot.set_xticks()
output.plot.set_yticks()

# Set axes labels
if config.xlabel is not None: output.plot.set_xlabel(config.xlabel)
if config.ylabel is not None: output.plot.set_ylabel(config.ylabel)

# Show
if config.path is not None: output.figure.saveto(config.path)
else: output.figure.show()

# -----------------------------------------------------------------
