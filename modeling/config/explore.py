#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition
from .explore_analysis import definition as analysis_definition

# -----------------------------------------------------------------

# Set the modeling path
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

modeling_config = environment.modeling_configuration
default_generation_method = modeling_config.fitting_method
generation_methods = ["genetic", "grid"]

# -----------------------------------------------------------------

definition = definition.copy()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Positional optional parameter
definition.add_positional_optional("generation_method", "string", "model generation method", default_generation_method, choices=generation_methods)

# -----------------------------------------------------------------

# Just a test
definition.add_flag("test", "just testing: don't create a generation")

# -----------------------------------------------------------------

## Remote execution

# Remote hosts
all_host_ids = find_host_ids()
has_remotes = len(all_host_ids) > 0
if has_remotes: definition.add_positional_optional("remotes", "string_list", "remote hosts to use", default=modeling_config.fitting_host_ids, choices=all_host_ids)
else: definition.add_fixed("remotes", [])
definition.add_flag("local", "run everything locally")

# Options
definition.add_flag("attached", "run remote simulations in attached mode")
definition.add_optional("nsimulations", "even_positive_integer", "number of simulations to launch in one batch/generation", 100)
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("group_walltime", "real", "the preferred walltime per group job (for schedulers)")

# -----------------------------------------------------------------

# Ranges
definition.add_flag("prompt_ranges", "prompt for the parameter ranges", True)
definition.add_flag("auto_ranges", "determine the parameter ranges automatically based on the previous generation")
definition.add_optional("range_probability", "percentage", "percentage of the total probability to use to select the new ranges automatically", "99.9", convert_default=True)

# -----------------------------------------------------------------

# Number of points per free parameter
definition.add_optional("npoints", "string_integer_dictionary", "number of grid points for the different free parameters")
definition.add_flag("prompt_npoints", "prompt for the number of grid points for each free parameter", False)
#definition.add_flag("auto_npoints", "determine the number of points automatically based on the previous generation")
definition.add_optional("npoints_all", "positive_integer", "number of points for all free parameters")

# -----------------------------------------------------------------

# Flags
definition.add_flag("visualise", "make visualisations")

# Model options
definition.add_optional("npackages_factor", "positive_real", "the factor with which to increase the number of photon packages for the new batch of simulations", 5.)
definition.add_flag("increase_npackages", "increase the number of photon packages with a certain factor", False)
definition.add_flag("adjust_npackages", "adjust the number of packages to the number of dust cells", True)
definition.add_optional("ncells_npackages_factor", "percentage", "relative number of photon packages w.r.t. the number of dust cells (used when adjust_npackages is enabled)", "75", convert_default=True)
definition.add_optional("npackages", "positive_real", "define the number of photon packages explicitly")

# Use a different wavelength grid or use a different representation
definition.add_flag("refine_spectral", "increase the resolution of the wavelength grid for the new batch of simulations", False)
definition.add_flag("refine_spatial", "increase the spatial resolution of the model for the new batch of simulations", False)
definition.add_flag("highres", "use high-resolution wavelength grids (default is same as previous generation)", None)
definition.add_optional("nwavelengths", "positive_integer", "target number of wavelengths of the desired wavelength grid")

# -----------------------------------------------------------------

# Simulation and analysis options
definition.add_flag("selfabsorption", "dust self-absorption", None)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", None)
definition.add_flag("spectral_convolution", "use spectral convolution for evaluating the simulations", None)
definition.add_flag("use_images", "use images for evaluating the simulations", None)
definition.add_flag("fit_not_clipped", "fit to the observed fluxes from the truncated (not clipped) images", False)
definition.add_flag("local_analysis", "perform all analysis locally", False)

# -----------------------------------------------------------------

# The number of parallel processes
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")
definition.add_flag("data_parallel_local", "use data-parallelization", False)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", None)
definition.add_flag("all_sockets", "use all sockets, not just the determined number of 'free' sockets (for remote execution)", False)
definition.add_optional("nsockets", "positive_integer", "use this number of sockets (for remote execution)")
definition.add_flag("allow_multisocket_processes", "allow using multiple sockets per process", False)

# -----------------------------------------------------------------

# Special options
definition.add_flag("keep", "keep remote input and output", True) # for certainty
definition.add_flag("attached", "run SKIRT in attached mode", False)
definition.add_flag("show_progress", "show progress of local simulations", False)
definition.add_flag("dry", "dry-run (don't actually launch simulations)", False)

# Advanced
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)

# -----------------------------------------------------------------

# ANALYSIS
definition.import_settings(analysis_definition)

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
definition.add_flag("deploy_clean", "perform clean installs when deploying (use with care!)", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------

definition.add_flag("plot", "make plots", False)
definition.add_optional("reference_images_colormap", "string", "colormap for plotting the reference images", "magma")
definition.add_flag("write_reference_images", "write the reference images and the plot (only when plotting is enabled)", True)

# -----------------------------------------------------------------
