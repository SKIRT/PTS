#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clone_fitting_run Clone a fitting run (without the generations).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.fitting.refitter import clone_fitting_run

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("fitting_run", "string", "fitting run to clone", choices=runs.names)
definition.add_required("name", "string", "name for the new fitting run", forbidden=runs.names)
config = parse_arguments("clone_fitting_run", definition)

# -----------------------------------------------------------------

# Get original fitting run
fitting_run = runs.load(config.fitting_run)

# -----------------------------------------------------------------

# Clone the fitting run
clone_fitting_run(fitting_run, config.name)

# -----------------------------------------------------------------
