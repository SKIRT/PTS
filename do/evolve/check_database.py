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

database = load_database(config.database)

runs = get_runs(database)

print("runs:", runs)

for run in runs:

    generations = get_generations(database, run)

    print(run)
    print("generations:", generations)

    populations = get_populations(database, run)

    print("populations:", populations)

# -----------------------------------------------------------------
