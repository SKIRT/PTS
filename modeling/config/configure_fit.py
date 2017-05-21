#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.build.component import get_model_names
from pts.modeling.fitting.component import get_run_names
from pts.modeling.component.component import get_default_fitting_method
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Fitting methods
#default_fitting_method = "genetic"
default_fitting_method = get_default_fitting_method(modeling_path) # as defined in the modeling configuration
fitting_methods = ["genetic", "grid"]

# Default number of digits
default_ndigits = 3

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Name for the fitting run
run_names = get_run_names(modeling_path)
if len(run_names) == 0: definition.add_positional_optional("name", "string", "name for the fitting run", default="run_1")
else: definition.add_required("name", "string", "name for the fitting run", forbidden=run_names)

# Name of the model to use
model_names = get_model_names(modeling_path)
if len(model_names) == 1: definition.add_fixed("model_name", "name of the model to use for the fitting", model_names[0])
else: definition.add_optional("model_name", "string", "name of the model to use for the fitting", choices=model_names)

# NEW: FITTING METHOD
definition.add_optional("fitting_method", "string", "fitting method", default_fitting_method, fitting_methods)

# Add optional
definition.add_optional("parameters", "string_list", "parameters to be used as free parameters during the fitting")
definition.add_optional("descriptions", "string_string_dictionary", "parameter descriptions")
definition.add_optional("types", "string_string_dictionary", "parameter types")
definition.add_optional("units", "string_unit_dictionary", "parameter units")
definition.add_optional("ndigits", "string_integer_dictionary", "number of significant digits of the parameters")
definition.add_optional("filters", "string_list", "fit to the observed data of these filters")
definition.add_optional("genetic", "dictionary", "options for the genetic algorithm")
definition.add_optional("grid", "dictionary", "options for the grid fitting")

# Sections
definition.add_section("ranges", "parameter ranges")

# The default number of significant digits
definition.add_optional("default_ndigits", "positive_integer", "default value for the number of significant digits", default_ndigits)

# -----------------------------------------------------------------
