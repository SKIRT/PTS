#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.component.component import get_default_fitting_method
from pts.modeling.build.suite import ModelSuite
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Default run name
default_run_name = "run_1"

# Fitting methods
default_fitting_method = get_default_fitting_method(modeling_path) # as defined in the modeling configuration
fitting_methods = ["genetic", "grid"]

# Default number of digits
default_ndigits = 3

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# FITTING RUN
if runs.empty: definition.add_positional_optional("name", "string", "name for the fitting run", default=default_run_name)
else: definition.add_required("name", "string", "name for the fitting run", forbidden=runs.names)

# -----------------------------------------------------------------

# MODEL
if suite.no_models: raise ValueError("No models are present (yet)")
elif suite.has_single_model: definition.add_fixed("model_name", "name of the model to use for the fitting", suite.single_model_name)
else: definition.add_optional("model_name", "string", "name of the model to use for the fitting", choices=suite.model_names)

# -----------------------------------------------------------------

# From configuration of other fitting run
if runs.has_any: definition.add_optional("from_run", "string", "adapt configuration of other fitting run", choices=runs.names)
else: definition.add_fixed("from_run", "no previous fitting run", None)

# -----------------------------------------------------------------

# Adapt certain things from the other fitting run
definition.add_flag("adapt_parameters", "adapt the free parameters from the other fitting run")
definition.add_flag("adapt_descriptions", "adapt the descriptions from the other fitting run")
definition.add_flag("adapt_types", "adapt the types from the other fitting run")
definition.add_flag("adapt_units", "adapt the units from the other fitting run")
definition.add_flag("adapt_ranges", "adapt the parameter ranges from the other fitting run")
definition.add_flag("adapt_ndigits", "adapt the ndigits from the other fitting run")
definition.add_flag("adapt_filters", "adapt the filters from the other fitting run")
definition.add_flag("adapt_genetic", "adapt the genetic settings from the other fitting run")
definition.add_flag("adapt_grid", "adapt settings for the grid fitting from the other fitting run")

# Don't use from other fitting run
definition.add_flag("except_parameters", "don't use parameters from other fitting run")
definition.add_flag("except_descriptions", "don't use descriptions from other fitting run")
definition.add_flag("except_types", "don't use types from other fitting run")
definition.add_flag("except_units", "don't use units from other fitting run")
definition.add_flag("except_ranges", "don't use ranges from the other fitting run")
definition.add_flag("except_ndigits", "don't use ndigits from other fitting run")
definition.add_flag("except_filters", "don't use filters from other fitting run")
definition.add_flag("except_genetic", "don't use genetic settings from other fitting run")
definition.add_flag("except_grid", "don't use grid settings from other fitting run")

# -----------------------------------------------------------------

# FITTING METHOD
definition.add_optional("fitting_method", "string", "fitting method", default_fitting_method, fitting_methods)

# -----------------------------------------------------------------

# Add optional
definition.add_optional("parameters", "string_list", "parameters to be used as free parameters during the fitting")
definition.add_optional("descriptions", "string_string_dictionary", "parameter descriptions")
definition.add_optional("types", "string_string_dictionary", "parameter types")
definition.add_optional("units", "string_unit_dictionary", "parameter units")
definition.add_optional("ndigits", "string_integer_dictionary", "number of significant digits of the parameters")
definition.add_optional("filters", "string_list", "fit to the observed data of these filters")
definition.add_optional("genetic", "dictionary", "options for the genetic algorithm")
definition.add_optional("grid", "dictionary", "options for the grid fitting")

# Use free ranges?
definition.add_flag("free_ranges", "use free ranges, instead of using initial values and a certain magnitude range around it", False)

# Sections
definition.add_section("ranges", "parameter ranges")

# The default number of significant digits
definition.add_optional("default_ndigits", "positive_integer", "default value for the number of significant digits", default_ndigits)

# -----------------------------------------------------------------
