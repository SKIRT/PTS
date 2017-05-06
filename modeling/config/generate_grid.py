#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
#from pts.modeling.fitting.run import get_free_parameter_labels, get_parameter_descriptions

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

default_scale = "logarithmic"
scales = ["linear", "logarithmic"]

# -----------------------------------------------------------------

# Get parameter info: NOT POSSIBLE WITHOUT FITTING RUN INFORMATION
#labels = get_free_parameter_labels(modeling_path)
#descriptions = get_parameter_descriptions(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The number of models
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 80)

# Scales
#definition.add_section("scales", "scales (linear/logarithmic) for the different free parameters")
#for label in get_free_parameter_labels(modeling_path): definition.add_optional(label + "_scale", "string", "scale for " + descriptions[label], default_scale, choices=scales)
definition.add_optional("scales", "string_string_dictionary", "scales (linear/logarithmic) for the different free parameters")

# NOW ASKED IN THE GRIDGENERATOR CLASS ITSELF
# Sampling of grid points
#definition.add_optional("most_sampled_parameters", "string_list", "free parameter(s) which get the most sampling points", choices=labels)
#definition.add_optional("sampling_weights", "weights", "relative sampling for the free parameters: " + ",".join(labels))

# -----------------------------------------------------------------
