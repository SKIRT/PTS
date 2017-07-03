#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.evolve.solve.extremizer import genetic_definition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Settings
#definition.add_optional("nruns", "positive_integer", "number of runs", 2)
definition.add_optional("ngenerations", "positive_integer", "number of generations", 10)
definition.add_optional("nindividuals", "even_integer", "number of individuals per generation", 100)

# Genetic settings
definition.import_section("genetic", "genetic algorithm settings", genetic_definition)

# Flags
definition.add_flag("plot", "plot", True)

# -----------------------------------------------------------------
