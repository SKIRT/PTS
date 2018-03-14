#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

default_parallelization = "2:1:2" # 2 cores, 1 process, 2 threads per core

# -----------------------------------------------------------------

definition = definition.copy()

# The output directory
definition.add_optional("simulation_path", "directory_path", "simulation directory")

#definition.add_optional("npackages", "positive_integer", "number of photon packages", 1e4)

definition.add_optional("output", "directory_path", "output directory", letter="o")

definition.add_flag("write", "do writing", True)
#definition.add_section("writing", "writing options")
#definition.sections["writing"].add_optional("tree_file_path", "string", "path for the dust grid tree file")

definition.add_flag("plot", "do plotting", True)

# Flags to calculate various aspects of the dust grid quality
definition.add_flag("quality", "get the quality of the dust grid", True)
definition.add_flag("projected_quality", "get the projected quality", True)
definition.add_flag("optical_depth_quality", "get the optical depth quality", True)
definition.add_flag("density_quality", "get the density quality", True)
definition.add_flag("dust_mass_quality", "get the dust mass quality", True)

# -----------------------------------------------------------------

definition.add_optional("tree_path", "string", "filepath for the tree grid file", "tree.dat")

# -----------------------------------------------------------------

# ADVANCED: parallelization
definition.add_optional("parallelization", "parallelization", "parallelization scheme", default_parallelization, convert_default=True)

# -----------------------------------------------------------------

definition.add_flag("plot_dust_cell_distribution", "...", False)

# -----------------------------------------------------------------
