#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.check_database Check the database.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.evolve.analyse.database import load_database, get_generations, get_runs, get_populations
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools.logging import setup_log

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_required("database", "file_path", "database path")

# Create the configuration
setter = ArgumentConfigurationSetter("check_database")
config = setter.run(definition)

# Set logging
log = setup_log("DEBUG")

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = load_fitting_run(modeling_path, config.fitting_run)

# -----------------------------------------------------------------

# Get generation names
generations = config.generations if config.generations is not None else fitting_run.genetic_generations

# -----------------------------------------------------------------
