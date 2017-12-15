#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

default_generation_method = "genetic"
generation_methods = ["genetic", "grid"]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Positional optional parameter
definition.add_positional_optional("generation_method", "string", "model generation method", default_generation_method, choices=generation_methods)

# -----------------------------------------------------------------

## Remote execution

# Remote hosts
if len(find_host_ids()) > 0: definition.add_optional("remotes", "string_list", "the remote hosts on which to run the parameter exploration", default=find_host_ids(schedulers=False), choices=find_host_ids(schedulers=False))
else: definition.add_fixed("remotes", "remote hosts", [])

# Options
definition.add_flag("attached", "run remote simulations in attached mode")
definition.add_optional("nsimulations", "even_positive_integer", "number of simulations to launch in one batch/generation", 100)
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# Number of points per free parameter
definition.add_optional("npoints", "string_integer_dictionary", "number of grid points for the different free parameters")
definition.add_flag("prompt_npoints", "prompt for the number of grid points for each free parameter")

# -----------------------------------------------------------------

# Flags
definition.add_flag("visualise", "make visualisations")

# Model options
definition.add_optional("npackages_factor", "positive_real", "the factor with which to increase the number of photon packages for the new batch of simulations", 5.)
definition.add_flag("increase_npackages", "increase the number of photon packages with a certain factor", False)

# Use a different wavelength grid or use a different representation
definition.add_flag("refine_spectral", "increase the resolution of the wavelength grid for the new batch of simulations", False)
definition.add_flag("refine_spatial", "increase the spatial resolution of the model for the new batch of simulations", False)
definition.add_flag("highres", "use high-resolution wavelength grids (default is same as previous generation)", None)

# -----------------------------------------------------------------

# Simulation and analysis options
definition.add_flag("selfabsorption", "dust self-absorption", None)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", None)
definition.add_flag("spectral_convolution", "use spectral convolution for evaluating the simulations", None)
definition.add_flag("use_images", "use images for evaluating the simulations", None)

# -----------------------------------------------------------------

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations (for scheduler)", 4)
definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 4)
definition.add_flag("data_parallel", "data parallelization mode", False)

# -----------------------------------------------------------------

# Special options
definition.add_flag("dry", "dry-run (don't actually launch simulations)")

# Advanced
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)

# Launcher settings
definition.add_flag("extract_progress", "extract progress information", False)
definition.add_flag("extract_timeline", "extract timeline information", False)
definition.add_flag("extract_memory", "extract memory information", False)
definition.add_flag("plot_progress", "plot progress information", False)
definition.add_flag("plot_timeline", "plot simulation timeline", False)
definition.add_flag("plot_memory", "plot memory information", False)
definition.add_flag("plot_seds", "plot the SEDs of individual simulations", False)

# Other
definition.add_flag("record_timing", "record timing information in timing table", True)
definition.add_flag("record_memory", "record memory information in memory table", True)

# VERY ADVANCED
definition.add_optional("restart_from_generation", "string", "restart everything from this generation")

# Check recurrence
definition.add_flag("check_recurrence", "check for recurrence of models that have been simulated previously", True)
definition.add_optional("recurrence_rtol", "positive_real", "relative tolerance for recurrence checking", 1e-5)
definition.add_optional("recurrence_atol", "positive_real", "absolute tolerance for recurrence checking", 1e-8)

# -----------------------------------------------------------------

# Deploy on remote hosts
definition.add_flag("deploy", "deploy SKIRT where necessary", False)
definition.add_flag("check_versions", "check versions of SKIRT where necessary", True)
#definition.add_flag("update_dependencies", "update PTS dependencies", False)
definition.add_flag("deploy_clean", "perform clean installs when deploying (use with care!)", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------
