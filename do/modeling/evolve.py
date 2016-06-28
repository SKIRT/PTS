#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.evolve Evolve the radiative transfer models.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.fitting.basicevolver import BasicModelEvolver
from pts.modeling.fitting.geneticevolver import GeneticModelEvolver
from pts.core.tools import logging, time
from pts.core.basics.errors import ConfigurationError
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Required parameters
config.add_required("method", str, "the evolution method (basic or genetic)", False)

# Optional parameters
config.add_optional("simulations", int, "the number of simulations to launch in one batch/generation", 100)
config.add_optional("remotes", "string_list", "the remote hosts on which to launch the simulations (None means all)", None)
config.add_optional("young", "float_tuple", "the range of the FUV luminosity of the young stellar population", (0.0, 4.e16))
config.add_optional("ionizing", "float_tuple", "the range of the FUV luminosity of the ionizing stellar population", (0.0, 5.e10))
config.add_optional("dust", "float_tuple", "the range of the dust mass", (0.5e7, 3.e7))

# Flags
config.add_flag("relative", "whether the range values are relative to the best (or initial) parameter value")
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
log.start("Starting evolve ...")

# -----------------------------------------------------------------

# Based on the evolution method, create either a BasicModelEvolver or a GeneticModelEvolver
if config.arguments.method == "basic": evolver = BasicModelEvolver(config.get_settings())
elif config.arguments.method == "genetic": evolver = GeneticModelEvolver(config.get_settings())
else: raise ConfigurationError("Wrong option for 'method'")

# Run the evolution
evolver.run()

# -----------------------------------------------------------------
