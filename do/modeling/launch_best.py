#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.launch_best Launch a simulation for the best fitting model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.analysis.launch import BestModelLauncher
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

config.add_section("wavelengths")
config.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
config.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.05)
config.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 1000)
config.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 300)
config.sections["wavelengths"].add_optional("min_zoom", float, "the minium wavelength of the zoomed-in grid", 1)
config.sections["wavelengths"].add_optional("max_zoom", float, "the maximum wavelength of the zoomed-in grid", 30)
config.sections["wavelengths"].add_optional("npoints_zoom", int, "the number of wavelength points in the zoomed-in grid", 300)

config.add_optional("packages", float, "the number of photon packages per wavelength", 1e6)
config.add_optional("selfabsorption", bool, "whether self-absorption should be enabled", True)
config.add_optional("remote", str, "the remote host on which to launch the simulations", "nancy")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting launch_best ...")

# -----------------------------------------------------------------

# Create a BestModelLauncher object
launcher = BestModelLauncher(config.get_settings())

# Run the launcher
launcher.run()

# -----------------------------------------------------------------
