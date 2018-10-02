#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.magic.tools import wavelengths
from pts.modeling.config.component import definition
from pts.core.launch.options import LoggingOptions, SchedulingOptions

# -----------------------------------------------------------------

# Get the analysis runs
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# -----------------------------------------------------------------

# Define the wavelength range
# Only dust absorption and temperature are relevant?
wavelength_range = wavelengths.wavelength_range_before_and_including("MIR")

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# THE ANALYSIS RUN
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_required("run", "string", "name of the analysis run", runs.names)

# -----------------------------------------------------------------

# Optional settings
definition.add_positional_optional("remote", "string", "remote host on which to launch the simulations", choices=find_host_ids())
definition.add_positional_optional("cluster_name", "string", "cluster of the remote host to use for the simulation")
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("group_walltime", "real", "the preferred walltime per group job (for schedulers)")

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("npackages", "real", "number of photon packages per wavelength", 1e7)
definition.add_optional("npackages_seds", "real", "number of photon packages per wavelength for simulations that produce only SEDs", 1e6)
definition.add_flag("selfabsorption", "dust self-absorption", True)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", True)

# -----------------------------------------------------------------

# The number of parallel processes for local execution
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")

# Enable data-parallelization?
definition.add_flag("data_parallel_local", "use data-parallelization", False)
# SET TO FALSE FOR NOW BECAUSE SOMETHING IS BROKEN: SKIRT CRASHES WHEN WRITING ABSORPTION DATA (ONLY IN DATA-PARALLEL MODE)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", False)  # was None

# Specify parallelization
definition.add_optional("parallelization_local", "parallelization", "parallelization scheme for local execution")
definition.add_optional("parallelization_remote", "parallelization", "parallelization scheme for remote execution")

# -----------------------------------------------------------------

definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)

# -----------------------------------------------------------------

# EXTRA SIMULATION OUTPUT
definition.add_flag("temperatures", "get dust temperature data from the simulation", False)
definition.add_flag("emissivities", "get dust emissivity data from the simulation", False)
definition.add_flag("isrf", "get ISRF data from the simulation", False)

# -----------------------------------------------------------------

# If previous launch failed, but input was already uploaded
definition.add_optional("remote_input", "string_string_dictionary", "dictionary of input file paths, the filenames (keys) have to be those defined in the input.dat file of the analysis run")
definition.add_optional("remote_input_path", "string", "remote directory where the uploaded input files can be found")

# -----------------------------------------------------------------

# Keep?
definition.add_flag("keep", "keep remote input and output", True) # for certainty

# -----------------------------------------------------------------
