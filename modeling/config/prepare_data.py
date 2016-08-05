#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add required arguments
definition.add_required("image", "string", "the name of the image for which to run the preparation")

# Add optional arguments
definition.add_optional("reference_image", "string", "the name of the reference image")
definition.add_flag("steps", "write the results of intermediate steps")
definition.add_flag("visualise", "make visualisations")

#config.add_section("importation")
#config.add_section("preparation")

# -----------------------------------------------------------------
