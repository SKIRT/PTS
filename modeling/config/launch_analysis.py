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
from pts.modeling.analysis.run import AnalysisRuns

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Optional settings
definition.add_optional("remote", "string", "remote host on which to launch the simulation", "nancy", choices=find_host_ids())
definition.add_optional("images_remote", "string", "the remote host on which to make the observed images", "nancy", choices=find_host_ids())

# ANALYSIS RUN
if runs.empty: raise RuntimeError("No analysis runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# Settings for the wavelength grid
definition.add_optional("nwavelengths", "integer", "the number of wavelengths to simulate the best model", 450)

# Settings for the dust grid
definition.add_section("dg", "options for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 1.)
definition.sections["dg"].add_optional("min_level", "integer", "the minimum division level for the tree", 8)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "the maximum mass fraction per cell", 1e-6)

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e7)

# Parallelization options
definition.add_optional("nnodes", "integer", "number of nodes to use for the simulations (for scheduler)", 4)
definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 10)
definition.add_flag("data_parallel", "data parallelization mode", False)

# -----------------------------------------------------------------
