#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.analysis.component import get_analysis_run_names, get_last_run_name
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Positional option
definition.add_positional_optional("run", "string", "name of the analysis run", get_last_run_name(fs.cwd()), get_analysis_run_names(fs.cwd()))

# -----------------------------------------------------------------
