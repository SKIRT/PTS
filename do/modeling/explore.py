#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.explore Explore the parameter space for the radiative transfer modeling.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.fitting.manualexplorer import ManualParameterExplorer
from pts.modeling.fitting.geneticexplorer import GeneticParameterExplorer
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration
from pts.core.basics.errors import ConfigurationError

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Required parameters
config.add_required("method", str, "the evolution method (manual or genetic)", False)

# Optional parameters
config.add_optional("remote", str, "the remote host on which to run the parameters exploration", "nancy")
config.add_optional("simulations", int, "the number of simulations to launch in one batch/generation", 100)
config.add_optional("young", "float_tuple", "the range of the FUV luminosity of the young stellar population", (0.0, 4.e16))
config.add_optional("ionizing", "float_tuple", "the range of the FUV luminosity of the ionizing stellar population", (0.0, 5.e10))
config.add_optional("dust", "float_tuple", "the range of the dust mass", (0.5e7, 3.e7))

# Flags
config.add_flag("young_log", "use logarithmic spacing of the young stellar luminosity values")
config.add_flag("ionizing_log", "use logarithmic spacing of the ionizing stellar luminosity values")
config.add_flag("dust_log", "use logarithmic spacing of the dust mass values")
config.add_flag("visualise", "make visualisations")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting explore ...")

# -----------------------------------------------------------------

# Based on the exploration method, create either a ManualParameterExplorer or a GeneticParameterExplorer instance
if config.arguments.method == "manual": explorer = ManualParameterExplorer(config.get_settings())
elif config.arguments.method == "genetic": explorer = GeneticParameterExplorer(config.get_settings())
else: raise ConfigurationError("Wrong option for 'method'")

# Run the parameter exploration
explorer.run()

# -----------------------------------------------------------------
