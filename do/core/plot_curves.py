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

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.curve import Curve
from pts.core.basics.plot import MPLPlot
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("files", "filepath_list", "list of the data file paths")
definition.add_flag("loglog", "use log-log scale", False)

# Parse
config = parse_arguments("plot_curves", definition)

# -----------------------------------------------------------------

curves = dict()
for filepath in config.files:
    name = fs.strip_extension(fs.name(filepath))
    curves[name] = Curve.from_data_file(filepath)

# -----------------------------------------------------------------

plot = MPLPlot()
for name in curves: plot.add_curve(curves[name], name)

# -----------------------------------------------------------------

if config.loglog:
    plot.set_x_log_scale()
    plot.set_y_log_scale()

# -----------------------------------------------------------------

# Plot
plot.finish()

# -----------------------------------------------------------------
