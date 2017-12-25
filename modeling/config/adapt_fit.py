#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

properties = ["representation", "filters", "ranges", "genetic", "grid", "units", "types"]

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# The fitting run for which to adapt the configuration
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# -----------------------------------------------------------------

# Dust or stellar
definition.add_positional_optional("properties", "string_list", "properties to adapt", default=properties, choices=properties)

# -----------------------------------------------------------------

# Select certain properties
definition.add_optional("contains", "string", "only adapt properties containing this string in their name")
definition.add_optional("not_contains", "string", "don't adapt properties containing this string in their name")
definition.add_optional("exact_name", "string", "only adapt properties with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't adapt properties with this exact string as their name")
definition.add_optional("startswith", "string", "only adapt properties whose name starts with this string")
definition.add_optional("endswith", "string", "only adapt properties whose name starts with this string")

# -----------------------------------------------------------------

# Save
definition.add_flag("save", "save adapted properties", True)

# -----------------------------------------------------------------
