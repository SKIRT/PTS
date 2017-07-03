#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.compare_sed 

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.data.sed import SED, ObservedSED
from pts.core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

# Set the log level
level = "DEBUG"

# Initialize the logger
log = logging.setup_log(level=level)
log.start("Starting compare_sed ...")

# -----------------------------------------------------------------

fit_path = fs.absolute_path("../../../../../")
modeling_path = fs.directory_of(fit_path)
galaxy_name = fs.name(modeling_path)

# -----------------------------------------------------------------

# Create the plotter
plotter = SEDPlotter()

# Add the SEDs
filename = galaxy_name + "_earth_sed.dat"
sed = SED.from_skirt(filename)

# Load the reference SED
reference_sed_path = fs.join(modeling_path, "sed.dat")
reference_sed = ObservedSED.from_file(reference_sed_path)

plotter.add_sed(sed, "simulation")
plotter.add_sed(reference_sed, "reference")

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
