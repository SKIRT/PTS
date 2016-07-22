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
definition.sections["cutoff"].add_optional("reference_path", "string", "...", None)
definition.sections["cutoff"].add_optional("level", "real", "cutoff when signal < level * uncertainty (ilse: 5)", 3.0)
definition.sections["cutoff"].add_optional("remove_holes", "boolean", "remove holes from the cutoff mask", True)

#definition.add_section("old_stars", "old stellar population")
#definition.sections["old_stars"].add_optional("irac_snr_level", float, "cut-off when signal(IRAC) < irac_snr_level * uncertainty(IRAC) (Ilse: 10.0)", 0.0)

# -----------------------------------------------------------------
