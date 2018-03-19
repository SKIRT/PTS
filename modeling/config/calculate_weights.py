#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.fitting.weights import all_regime_names

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The filters
definition.add_positional_optional("filters", "lazy_broad_band_filter_list", "list of filters")

# Regimes
definition.add_optional("regimes", "string_list", "wavelength regimes to use", choices=all_regime_names)

# Use physical regimes -> now the default
definition.add_flag("physical", "use physical regimes", True)

## FOR ORIGINAL REGIMES

# Weights of regimes
definition.add_optional("uv", "positive_real", "default relative weight for UV bands", 1.)
definition.add_optional("optical", "positive_real", "default relative weight for optical bands", 1.)
definition.add_optional("nir", "positive_real", "default relative weight for NIR bands", 1.)
definition.add_optional("mir", "positive_real", "default relative weight for MIR bands", 1.)
definition.add_optional("fir", "positive_real", "defult relative weight for FIR bands", 1.)
definition.add_optional("submm_microwave", "positive_real", "default relative weight for submm/microwave bands", 1.)

# Flags
definition.add_flag("only_uv", "only give weight to UV bands")
definition.add_flag("only_optical", "only give weight to optical bands")
definition.add_flag("only_nir", "only give weight to NIR bands")
definition.add_flag("only_mir", "only give weight to MIR bands")
definition.add_flag("only_fir", "only give weight to FIR bands")
definition.add_flag("only_submm_microwave", "only give weight to submm/microwave bands")

# Flags
definition.add_flag("no_uv", "give no weight to UV bands")
definition.add_flag("no_optical", "give no weight to optical bands")
definition.add_flag("no_nir", "give no weight to NIR bands")
definition.add_flag("no_mir", "give no weight to MIR bands")
definition.add_flag("no_fir", "give no weight to FIR bands")
definition.add_flag("no_submm_microwave", "give no weight to submm/microwave bands")

## FOR PHYSICAL REGIMES

# Weights
definition.add_optional("ionizing", "positive_real", "default relative weight for ionizing stellar regime", 1.)
definition.add_optional("young", "positive_real", "default relative weight for young stellar regime", 1.)
definition.add_optional("evolved", "positive_real", "default relative weight for evolved stellar regime", 1.)
definition.add_optional("mix", "positive_real", "default relative weight for mix regime", 1.)
definition.add_optional("aromatic", "positive_real", "default relative weight for aromatic dust regime", 1.)
definition.add_optional("thermal", "positive_real", "default relative weight for thermal dust regime", 1.)

# Flags
definition.add_flag("only_ionizing", "only give weight to ionizing regime")
definition.add_flag("only_young", "only give weight to young regime")
definition.add_flag("only_evolved", "only give weight to evolved regime")
definition.add_flag("only_mix", "only give weight to mix regime")
definition.add_flag("only_aromatic", "only give weight to aromatic regime")
definition.add_flag("only_thermal", "only give weight to thermal regime")

# Flags
definition.add_flag("no_ionizing", "give no weight to ionizing regime")
definition.add_flag("no_young", "give no weight to young regime")
definition.add_flag("no_evolved", "give no weight to evolved regime")
definition.add_flag("no_mix", "give no weight to mix regime")
definition.add_flag("no_aromatic", "give no weight to aromatic regime")
definition.add_flag("no_thermal", "give no weight to thermal regime")

# Write the weights?
definition.add_flag("write", "write the weights table", True)

# Show?
definition.add_flag("show", "show the weights", True)

# -----------------------------------------------------------------
