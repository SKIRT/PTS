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

definition.add_section("cutoff", "options for cutting off the maps at certain noise levels")
definition.sections["cutoff"].add_optional("reference_path", str, "...", None)
definition.sections["cutoff"].add_optional("level", float, "cutoff when signal < level * uncertainty (ilse: 5)", 3.0)
definition.sections["cutoff"].add_optional("remove_holes", bool, "remove holes from the cutoff mask", True)

#definition.add_section("non_ionizing_stars")
#definition.sections["non_ionizing_stars"].add_optional("fuv_snr_level", float, "cut-off when signal(FUV) < fuv_snr_level * uncertainty(FUV) (Ilse: 10.0)", 0.0)

# -----------------------------------------------------------------
