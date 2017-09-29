#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import load_modeling_environment_cwd

# -----------------------------------------------------------------

# Load the modeling environment and analysis runs
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# Determine remote hosts for simulations
host_ids = find_host_ids()
host_ids_suggestions = find_host_ids(schedulers=False) # suggest to not use schedulers: requires runtime to be estimated

# Determine the remote hosts for heavy computations (making the images)
config = environment.modeling_configuration
other_host_ids = config.host_ids
if other_host_ids is None or len(other_host_ids) == 0: other_host_ids = find_host_ids(schedulers=False)
default_other_host_id = other_host_ids[0]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Remote execution
definition.add_optional("remote", "string", "remote host on which to launch the simulation", choices=host_ids)
definition.add_optional("cluster_name", "string", "cluster of the remote host to use for the simulation")
definition.add_optional("images_remote", "string", "remote host on which to make the observed images", default_other_host_id, choices=other_host_ids)
definition.add_flag("attached", "launch remote executions in attached mode", True)
definition.add_flag("debug_output", "show all simulation output when in debugging mode", False)
definition.add_flag("keep_remote_input", "keep the remote input directory after the simulation is retrieved", True)  # KEEP REMOTE INPUT FOR NEW LAUNCHES AND FOR HEATING LAUNCHES
definition.add_flag("keep_remote_input_and_output", "keep the remote input and output after the simulation is retrieved", False)

# ANALYSIS RUN
if runs.empty: raise RuntimeError("No analysis runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# Parallelization options
#definition.add_optional("nnodes", "integer", "number of nodes to use for the simulations (for scheduler)", 4)
#definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 10)
#definition.add_flag("data_parallel", "data parallelization mode", False)

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e7)
definition.add_flag("selfabsorption", "dust self-absorption", True)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", True)

# -----------------------------------------------------------------

# ANALYSIS OPTIONS

# Remote rebinning/convolution
definition.add_optional("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)
definition.add_optional("convolve_remote_threshold", "data_quantity", "data size threshold for remote convolution", "1. GB", convert_default=True)

# For creating observed images
definition.add_optional("images_unit", "photometric_unit", "unit for the mock observed images", "Jy", convert_default=True)

# Spectral convolution
definition.add_flag("spectral_convolution_fluxes", "enable spectral convolution for creating observed fluxes", True)
definition.add_flag("spectral_convolution_images", "enable spectral convolution for creating observed images", True)

# -----------------------------------------------------------------

# Get temperature data
definition.add_flag("temperatures", "get dust temperature data from the simulation", True)

# Retrieval
definition.add_flag("retrieve_contributions", "retrieve datacubes for the different contributions to the flux", True)

# -----------------------------------------------------------------

# The number of parallel processes for local execution
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")
definition.add_flag("data_parallel_local", "use data-parallelization", False)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", None)

# Specify parallelization
definition.add_optional("parallelization_local", "parallelization", "parallelization scheme for local execution")
definition.add_optional("parallelization_remote", "parallelization", "parallelization scheme for remote execution")

# -----------------------------------------------------------------

# DVANCED: Clear remotes
definition.add_flag("clear_remotes", "clear temporary data and sessions on all remotes (use with care!)", False)

# Deploy SKIRT and PTS
definition.add_flag("deploy", "deploy SKIRT and PTS where necessary", False)
definition.add_flag("check_versions", "check versions of SKIRT and PTS where necessary", True)
definition.add_flag("update_dependencies", "update PTS dependencies", False)
definition.add_flag("deploy_clean", "perform clean installs when deploying (use with care!)", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# If previous launch failed, but input was already uploaded
definition.add_optional("remote_input", "string_string_dictionary", "dictionary of input file paths, the filenames (keys) have to be those defined in the input.dat file of the analysis run")
definition.add_optional("remote_input_path", "string", "remote directory where the uploaded input files can be found")

# When the number of dust cells from the tree is known, so that tree doesn't need to be loaded
definition.add_optional("ncells", "positive_integer", "number of dust cells")

# -----------------------------------------------------------------
