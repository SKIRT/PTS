#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.component import get_run_names
from pts.core.remote.host import find_host_ids
from pts.modeling.fitting.run import get_fitting_method

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Model generation methods
#generation_methods = ["genetic", "grid", "instinctive"]

# Default model generation method
#default_generation_method = "genetic"

# -----------------------------------------------------------------

# Get the fitting method
fitting_method = get_fitting_method(modeling_path)

# Set generation methods and default generation method, based on the fitting method
if fitting_method == "genetic":

    default_generation_method = "genetic"
    generation_methods = ["genetic", "grid"]

elif fitting_method == "grid":

    default_generation_method = "grid"
    generation_methods = ["grid"]

else: raise ValueError("Fitting method has an invalid value: " + fitting_method + " (must be 'genetic' or 'grid'")

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The fitting run for which to explore the parameter space
run_names = get_run_names(modeling_path)
if len(run_names) == 0: raise RuntimeError("No fitting runs found: first run configure_fit to create a new fitting run")
elif len(run_names) == 1: definition.add_fixed("name", "name of the fitting run", run_names[0])
else: definition.add_required("name", "string", "name of the fitting run", choices=run_names)

# Positional optional parameter
definition.add_positional_optional("generation_method", "string", "model generation method", default_generation_method, choices=generation_methods)

# Optional parameters
if len(find_host_ids()) > 0: definition.add_optional("remotes", "string_list", "the remote hosts on which to run the parameter exploration", default=find_host_ids(schedulers=False), choices=find_host_ids(schedulers=False))
else: definition.add_fixed("remotes", "remote hosts", [])
definition.add_flag("attached", "run remote simulations in attached mode")
definition.add_optional("nsimulations", "even_positive_integer", "the number of simulations to launch in one batch/generation", 100)
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# Flags
definition.add_flag("visualise", "make visualisations")

# Model options
definition.add_optional("npackages_factor", "positive_real", "the factor with which to increase the number of photon packages for the new batch of simulations", 5.)
definition.add_flag("increase_npackages", "increase the number of photon packages with a certain factor", False)
definition.add_flag("refine_spectral", "increase the resolution of the wavelength grid for the new batch of simulations", False)
definition.add_flag("refine_spatial", "increase the spatial resolution of the model for the new batch of simulations", False)
definition.add_flag("selfabsorption", "dust self-absorption", None)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", None)

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations (for scheduler)", 4)
definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 4)
definition.add_flag("data_parallel", "data parallelization mode", False)

# Special options
definition.add_flag("dry", "dry-run (don't actually launch simulations)")

# Advanced
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)

# -----------------------------------------------------------------
