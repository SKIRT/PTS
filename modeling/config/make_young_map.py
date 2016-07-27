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
definition = ConfigurationDefinition()

#definition.add_section("cutoff", "options for cutting off the maps at certain noise levels")
#definition.sections["cutoff"].add_optional("reference_path", "string", "...", None)
#definition.sections["cutoff"].add_optional("level", "real", "cutoff when signal < level * uncertainty (ilse: 5)", 3.0)
#definition.sections["cutoff"].add_optional("remove_holes", "boolean", "remove holes from the cutoff mask", True)

#definition.add_section("non_ionizing_stars")
#definition.sections["non_ionizing_stars"].add_optional("fuv_snr_level", float, "cut-off when signal(FUV) < fuv_snr_level * uncertainty(FUV) (Ilse: 10.0)", 0.0)

definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor that denotes the contribution of the old stellar population to the FUV emission", "0.1,0.4", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 4)
definition.add_optional("best_factor", "real", "the best estimate for the value of the factor", 0.2)

# -----------------------------------------------------------------
