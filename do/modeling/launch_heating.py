#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.launch_heating Launch simulations for analysing the heating contributions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.analysis.heating.launch import DustHeatingContributionLauncher
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

config.add_section("wavelengths")
config.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
config.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.1)
config.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 10)
config.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 25)

config.add_optional("packages", float, "the number of photon packages per wavelength", 1e7)
config.add_optional("remotes", "string_list", "the list of remote hosts on which to launch the simulations", ["nancy"])

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting launch_heating ...")

# -----------------------------------------------------------------

# Create a DustHeatingContributionLauncher object
launcher = DustHeatingContributionLauncher(config.get_settings())

# Run the launcher
launcher.run()

# -----------------------------------------------------------------
