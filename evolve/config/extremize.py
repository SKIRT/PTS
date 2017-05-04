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

min_or_max = ["minimize", "maximize"]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Genetic settings
definition.import_section("genetic", "genetic algorithm settings", genetic_definition)

# Add
definition.sections["genetic"].add_optional("ngenerations", "positive_integer", "number of generations", 10)
definition.sections["genetic"].add_optional("nindividuals", "even_integer", "number of individuals per generation", 100)

# Minimize or maximize
definition.add_required("min_or_max", "string", "minimize or maximize", choices=min_or_max)

# -----------------------------------------------------------------
