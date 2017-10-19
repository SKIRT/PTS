#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# The significance level
#definition.add_optional("fuv_significance", "real", "significance level of the FUV image below which to cut-off the dust map", 2.5)
#definition.add_optional("mips24_significance", "real", "significance level of the MIPS 24 micron image below which to cut-off the dust map", 2.0)
#definition.add_optional("pacs70_significance", "real", "significance level of the Pacs 70 micron image below which to cut-off the dust map", 1.0)
#definition.add_optional("pacs160_significance", "real", "significance level of the Pacs 160 micron image below which to cut-off the dust map", 2.0)
#definition.add_optional("h_significance", "real", "significance level of the 2MASS H image below which to cut-off the dust map", 0.0) # used for SSFR

# Remove holes from the cutoff mask
#definition.add_flag("remove_holes", "remove holes from the total cutoff mask", True)

# Flags for enabling/disabling different methods
definition.add_flag("make_black_body", "make dust map based on black-body fitting", True)
definition.add_flag("make_emission", "make dust map based on emission", True)
definition.add_flag("make_attenuation", "make dust map based on attenuation", True)
definition.add_flag("make_hot", "make map of hot dust", True)

definition.add_optional("hot_factor_range", "real_range", "range of factor to create the hot dust maps", "0.2>0.7", convert_default=True)
definition.add_optional("factor_nvalues", "positive_integer", "number of factors", 8)

# Best method
#definition.add_optional("best_method", "string", "the method of which to use the resul as the final dust map", "cortese")

# Sections: different dust map makers
#definition.import_section("cortese", "options for Cortese dust map maker", cortese_definition)
#definition.import_section("buat", "options for Buat dust map maker", buat_definition)
#definition.import_section("black_body", "options for black body dust map maker", bb_definition)
#definition.import_section("emission", "options for emission dust map maker", emission_definition)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------
