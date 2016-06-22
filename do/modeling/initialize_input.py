#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.initialize_input Initialize the input directory for the fitting procedure.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.fitting.initialization import InputInitializer
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add optional arguments
config.add_section("wavelengths")
config.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
config.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.09)
config.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 2000)
config.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 100)
config.sections["wavelengths"].add_optional("min_zoom", float, "the minimum wavelength of the zoomed-in grid", 1)
config.sections["wavelengths"].add_optional("max_zoom", float, "the maximum wavelength of the zoomed-in grid", 30)
config.sections["wavelengths"].add_optional("npoints_zoom", int, "the number of wavelength points in the zoomed-in grid", 100)

config.add_optional("packages", float, "the number of photon packages per wavelength", 2e5)
config.add_flag("selfabsorption", "enable dust self-absorption")
config.add_optional("dust_grid", str, "the type of dust grid to use (bintree, octtree or cartesian)", "bintree")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("initialization") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting initialize_input ...")

# -----------------------------------------------------------------

# Create a InputInitializer object
initializer = InputInitializer(config.get_settings())

# Run the input initialization
initializer.run()

# -----------------------------------------------------------------
