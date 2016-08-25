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

# The significance level
definition.add_optional("fuv_significance", "real", "significance level of the FUV image below which to cut-off the dust map", 2.5)
definition.add_optional("mips24_significance", "real", "significance level of the MIPS 24 micron image below which to cut-off the dust map", 2.0)
definition.add_optional("pacs70_significance", "real", "significance level of the Pacs 70 micron image below which to cut-off the dust map", 1.0)
definition.add_optional("pacs160_significance", "real", "significance level of the Pacs 160 micron image below which to cut-off the dust map", 2.0)
definition.add_optional("h_significance", "real", "significance level of the 2MASS H image below which to cut-off the dust map", 0.0) # used for SSFR

# Remove holes from the cutoff mask
definition.add_flag("remove_holes", "remove holes from the total cutoff mask", True)

definition.add_flag("make_black_body", "make dust map based on black-body fitting", True)
definition.add_flag("make_emission", "make dust map based on emission", True)
definition.add_flag("make_buat", "make dust map based on Buat", True)
definition.add_flag("make_cortese", "make dust map based on Cortese", True)

definition.add_optional("best_method", "string", "the method of which to use the resul as the final dust map", "cortese")

# -----------------------------------------------------------------
