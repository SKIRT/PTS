#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.core.environment import verify_modeling_cwd
from pts.magic.tools import wavelengths

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

# Define the wavelength range
# Only dust absorption and temperature are relevant?
wavelength_range = wavelengths.wavelength_range_before_and_including("MIR")

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# THE ANALYSIS RUN
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run for which to launch the heating simulations", runs.last_name, runs.names)

# -----------------------------------------------------------------

# Optional settings
definition.add_optional("remote", "string", "remote host on which to launch the simulations", choices=find_host_ids())
definition.add_optional("cluster_name", "string", "cluster of the remote host to use for the simulation")
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e7)
definition.add_flag("selfabsorption", "dust self-absorption", True)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", True)

# -----------------------------------------------------------------

# Settings for the wavelength grid
definition.add_section("wg", "options for the wavelength grid")
definition.sections["wg"].add_optional("range", "quantity_range", "the wavelength range", wavelength_range)
definition.sections["wg"].add_optional("npoints", "integer", "the number of wavelength points", 50)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)

# -----------------------------------------------------------------

definition.add_flag("use_grid_resolution", "use resolution of the grid for determining instrument resolution", None) # None means it is asked interactively
definition.add_optional("resolution_reference", "string", "use the map with this name as the reference for the instrument resolution", None) # None means it is asked interactively

# -----------------------------------------------------------------

# The number of parallel processes for local execution
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")
definition.add_flag("data_parallel_local", "use data-parallelization", False)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", None)

# Specify parallelization
#definition.add_optional("parallelization_local", "parallelization", "parallelization scheme for local execution")
#definition.add_optional("parallelization_remote", "parallelization", "parallelization scheme for remote execution")
# SET TO FALSE FOR NOW BECAUSE SOMETHING IS BROKEN: SKIRT CRASHES WHEN WRITING ABSORPTION DATA (ONLY IN DATA-PARALLEL MODE)
definition.add_optional("parallelization_local", "parallelization", "parallelization scheme for local execution", False)
definition.add_optional("parallelization_remote", "parallelization", "parallelization scheme for remote execution", False)

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 15.)

# -----------------------------------------------------------------

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations", 4)
definition.add_flag("data_parallel", "data parallelization mode", False)

# -----------------------------------------------------------------

# EXTRA SIMULATION OUTPUT
definition.add_flag("temperatures", "get dust temperature data from the simulation", False)
definition.add_flag("emissivities", "get dust emissivity data from the simulation", False)
definition.add_flag("isrf", "get ISRF data from the simulation", False)

# -----------------------------------------------------------------
