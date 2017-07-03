#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.reporting.reporting import steps

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The output directory
definition.add_optional("simulation_path", "directory_path", "simulation directory")

#definition.add_optional("npackages", "positive_integer", "number of photon packages", 1e4)

definition.add_optional("output", "directory_path", "output directory", letter="o")

definition.add_flag("write", "do writing", True)
#definition.add_section("writing", "writing options")
#definition.sections["writing"].add_optional("tree_file_path", "string", "path for the dust grid tree file")

definition.add_flag("plot", "do plotting", True)

# -----------------------------------------------------------------
