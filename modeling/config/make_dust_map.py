#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# Flags for enabling/disabling different methods
definition.add_flag("make_black_body", "make dust map based on black-body fitting", False)
definition.add_flag("make_emission", "make dust map based on emission", False)
definition.add_flag("make_attenuation", "make dust map based on attenuation", True)
definition.add_flag("make_hot", "make map of hot dust", True)

# Old stellar contribution subtraction factor
definition.add_optional("hot_factor_range", "real_range", "range of factor to create the hot dust maps", "0.2>0.7", convert_default=True)
definition.add_optional("factor_nvalues", "positive_integer", "number of factors", 8)

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

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------

# Old stars component
old_components = ["bulge", "disk", "total"]
definition.add_optional("old_component", "string", "old stellar component to use to subtract diffuse emission by evolved stars", "disk")

# -----------------------------------------------------------------
