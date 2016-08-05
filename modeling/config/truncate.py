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

definition.add_optional("image_names", "string_list", "The image names for which to make the truncated images", ["GALEX FUV", "SDSS g", "IRAC I1", "MIPS 24mu", "Pacs red", "SPIRE PLW"])

definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor to which the disk ellipse should be scaled to obtain the truncation region", "0.5,1.2", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 8)
definition.add_optional("best_factor", "real", "the best estimate for the value of the factor", 0.82)

# -----------------------------------------------------------------
