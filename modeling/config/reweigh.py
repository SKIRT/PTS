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

definition = definition.copy()

# -----------------------------------------------------------------

# Name for the reweighing
definition.add_required("name", "string", "name for the reweiging")

# The fitting run for which to adapt the configuration
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# Generations?
definition.add_positional_optional("generations", "string_list", "generation names")

# -----------------------------------------------------------------

wavelength_regimes = ["uv", "optical", "nir", "mir", "fir", "submm"]

# -----------------------------------------------------------------

definition.add_optional("filters", "filter_list", "filters to use for the evaluation (None means default fitting filters)")
definition.add_optional("regimes", "string_list", "wavelength regimes to use", default=wavelength_regimes, choices=wavelength_regimes)

# -----------------------------------------------------------------

# Flags
definition.add_flag("show", "show")
definition.add_flag("plot", "make plots")

# -----------------------------------------------------------------

definition.add_optional("uv", "positive_real", "relative weight for UV bands", 1.)
definition.add_optional("optical", "positive_real", "relative weight for optical bands", 1.)
definition.add_optional("nir", "positive_real", "relative weight for NIR bands", 1.)
definition.add_optional("mir", "positive_real", "relative weight for MIR bands", 1.)
definition.add_optional("fir", "positive_real", "relative weight for FIR bands", 1.)
definition.add_optional("submm", "positive_real", "relative weight for submm bands", 1.)

# -----------------------------------------------------------------

definition.add_flag("only_uv", "only give weight to UV bands")
definition.add_flag("only_optical", "only give weight to optical bands")
definition.add_flag("only_nir", "only give weight to NIR bands")
definition.add_flag("only_mir", "only give weight to MIR bands")
definition.add_flag("only_fir", "only give weight to FIR bands")
definition.add_flag("only_submm", "only give weight to submm bands")

# -----------------------------------------------------------------

definition.add_flag("no_uv", "give no weight to UV bands")
definition.add_flag("no_optical", "give no weight to optical bands")
definition.add_flag("no_nir", "give no weight to NIR bands")
definition.add_flag("no_mir", "give no weight to MIR bands")
definition.add_flag("no_fir", "give no weight to FIR bands")
definition.add_flag("no_submm", "give no weight to submm bands")

# -----------------------------------------------------------------
