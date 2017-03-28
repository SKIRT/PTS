#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.tests.NGC4013.test import fitting_filter_names, free_parameters, free_parameter_labels

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Simulation
definition.add_optional("npackages", "positive_integer", "number of photon packages", int(1e4))

# Free parameters
definition.add_optional("free_parameters", "string_list", "list of free parameters", free_parameter_labels, choices=free_parameters)

# Fitting filters
definition.add_optional("fitting_filters", "filter_list", "filters for the fitting", fitting_filter_names, convert_default=True)

# Genetic algorithm settings
definition.add_optional("ngenerations", "positive_integer", "number of generations", 5)
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 10)
definition.add_optional("mutation_rate", "positive_real", "mutation rate", 0.03)
definition.add_optional("crossover_rate", "positive_real", "crossover rate", 0.65)

# Other
definition.add_optional("kernel_type", "string", "convolution kernel type", "Gaussian")
definition.add_optional("kernel_dimension", "positive_integer", "kernel dimension", 6)

# -----------------------------------------------------------------
