#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
#from pts.core.tools.stringify import stringify_not_list
from pts.core.tools import filesystem as fs
#from pts.modeling.component.component import load_fitting_configuration
from pts.modeling.fitting.component import get_run_names
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Load the fitting configuration
#fitting_configuration = load_fitting_configuration(fs.cwd())

# Free parameter labels
#free_parameter_labels = fitting_configuration.free_parameters

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The fitting run for which to explore the parameter space
run_names = get_run_names(modeling_path)
if len(run_names) == 0: raise RuntimeError("No fitting runs found: first run configure_fit to create a new fitting run")
elif len(run_names) == 1: definition.add_fixed("name", "name of the fitting run", run_names[0])
else: definition.add_required("name", "string", "name of the fitting run", choices=run_names)

# Positional optional parameter
definition.add_positional_optional("generation_method", "string", "the model generation method ('grid', 'instinctive', 'genetic')", "genetic", ["genetic", "grid", "instinctive"])

# Optional parameters
if len(find_host_ids()) > 0: definition.add_optional("remotes", "string_list", "the remote hosts on which to run the parameter exploration", default=find_host_ids(schedulers=False), choices=find_host_ids(schedulers=False))
else: definition.add_fixed("remotes", "remote hosts", [])
definition.add_flag("attached", "run remote simulations in attached mode")
definition.add_optional("nsimulations", "even_positive_integer", "the number of simulations to launch in one batch/generation", 100)
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# Advanced options for the genetic engine
definition.add_optional("crossover_rate", "fraction", "the crossover rate", 0.5)
definition.add_optional("mutation_rate", "fraction", "the mutation rate", 0.5)

# The ranges of the different free parameters (although the absolute ranges are defined in the fitting configuration,
# give the option to refine these ranges for each exploration step (generation))
#for label in free_parameter_labels:

    # Get the default range
    #default_range = fitting_configuration[label + "_range"]
    #ptype, string = stringify_not_list(default_range)

    # Add the optional range setting for this free parameter
    #definition.add_optional(label + "_range", ptype, "the range of " + label, default_range)

# Flags
definition.add_flag("relative", "whether the range values are relative to the best (or initial) parameter value")
definition.add_flag("young_log", "use logarithmic spacing of the young stellar luminosity values")
definition.add_flag("ionizing_log", "use logarithmic spacing of the ionizing stellar luminosity values")
definition.add_flag("dust_log", "use logarithmic spacing of the dust mass values")
definition.add_flag("visualise", "make visualisations")

# Model options
definition.add_optional("npackages_factor", "positive_real", "the factor with which to increase the number of photon packages for the new batch of simulations", 5.)
definition.add_flag("increase_npackages", "increase the number of photon packages with a certain factor", False)
definition.add_flag("refine_wavelengths", "increase the resolution of the wavelength grid for the new batch of simulations", False)
definition.add_flag("refine_dust", "increase the resolution of the dust cell grid for the new batch of simulations", False)
definition.add_flag("selfabsorption", "dust self-absorption", None)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", None)

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations (for scheduler)", 4)
definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 4)
definition.add_flag("data_parallel", "data parallelization mode", False)

# Special options
definition.add_flag("dry", "dry-run (don't actually launch simulations)")

# -----------------------------------------------------------------
