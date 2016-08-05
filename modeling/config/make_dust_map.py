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

definition.add_section("cutoff", "options for cutting off the maps at certain noise levels")
definition.sections["cutoff"].add_optional("reference_path", "string", "...", None)
definition.sections["cutoff"].add_optional("level", "real", "cutoff when signal < level * uncertainty (ilse: 5)", 3.0)
definition.sections["cutoff"].add_optional("remove_holes", "boolean", "remove holes from the cutoff mask", True)

#definition.add_section("dust")
#definition.sections["dust"].add_section("ssfr")
#definition.sections["dust"].sections["ssfr"].add_optional("mask_low_fuv_snr", bool, "...", True)
#definition.sections["dust"].sections["ssfr"].add_optional("fuv_snr_level", float, "cut-off when signal(FUV) < fuv_snr_level * uncertainty(FUV)  (Ilse: 10.0)", 0.0)

definition.add_flag("make_black_body", "make dust map based on black-body fitting", True)
definition.add_flag("make_emission", "make dust map based on emission", True)
definition.add_flag("make_buat", "make dust map based on Buat", True)
definition.add_flag("make_cortese", "make dust map based on Cortese", True)

definition.add_optional("best_method", "string", "the method of which to use the resul as the final dust map", "cortese")

# -----------------------------------------------------------------
