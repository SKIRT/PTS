#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration("prepare_data", log_path="log")

# Add required arguments
config.add_required("image", str, "the name of the image for which to run the preparation")

# Add optional arguments
config.add_optional("reference_image", str, "the name of the reference image")
config.add_flag("steps", "write the results of intermediate steps")
config.add_flag("visualise", "make visualisations")

#config.add_section("importation")
#config.add_section("preparation")

# -----------------------------------------------------------------
