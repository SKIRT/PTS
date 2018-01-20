#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.calculate_weights Calculate weights for fitting.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.modeling.fitting.initialization.base import calculate_weights_filters
from pts.modeling.fitting.tables import WeightsTable
from pts.core.tools import sequences
from pts.modeling.fitting.initialization.base import wavelength_regimes

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The filters
definition.add_required("filters", "broad_band_filter_list", "list of filters")

# Regimes
definition.add_optional("regimes", "string_list", "wavelength regimes to use", default=wavelength_regimes, choices=wavelength_regimes)

definition.add_optional("uv", "positive_real", "default relative weight for UV bands", 1.)
definition.add_optional("optical", "positive_real", "default relative weight for optical bands", 1.)
definition.add_optional("nir", "positive_real", "default relative weight for NIR bands", 1.)
definition.add_optional("mir", "positive_real", "default relative weight for MIR bands", 1.)
definition.add_optional("fir", "positive_real", "defult relative weight for FIR bands", 1.)
definition.add_optional("submm_microwave", "positive_real", "default relative weight for submm/microwave bands", 1.)

definition.add_flag("only_uv", "only give weight to UV bands")
definition.add_flag("only_optical", "only give weight to optical bands")
definition.add_flag("only_nir", "only give weight to NIR bands")
definition.add_flag("only_mir", "only give weight to MIR bands")
definition.add_flag("only_fir", "only give weight to FIR bands")
definition.add_flag("only_submm_microwave", "only give weight to submm/microwave bands")

definition.add_flag("no_uv", "give no weight to UV bands")
definition.add_flag("no_optical", "give no weight to optical bands")
definition.add_flag("no_nir", "give no weight to NIR bands")
definition.add_flag("no_mir", "give no weight to MIR bands")
definition.add_flag("no_fir", "give no weight to FIR bands")
definition.add_flag("no_submm_microwave", "give no weight to submm/microwave bands")

# Get configuration
config = parse_arguments("calculate_weights", definition)

# -----------------------------------------------------------------

# Create the table to contain the weights
table = WeightsTable()

# -----------------------------------------------------------------

# Only UV
if config.only_uv:

    if config.no_uv: raise ValueError("Error")
    if config.only_optical: raise ValueError("Error")
    if config.only_nir: raise ValueError("Error")
    if config.only_mir: raise ValueError("Error")
    if config.only_fir: raise ValueError("Error")
    if config.only_submm_microwave: raise ValueError("Error")
    regimes = ["uv"]

# Only optical
elif config.only_optical:

    if config.no_optical: raise ValueError("Error")
    if config.only_uv: raise ValueError("Error")
    if config.only_nir: raise ValueError("Error")
    if config.only_mir: raise ValueError("Error")
    if config.only_fir: raise ValueError("Error")
    if config.only_submm_microwave: raise ValueError("Error")
    regimes = ["optical"]

# Only NIR
elif config.only_nir:

    if config.no_nir: raise ValueError("Error")
    if config.only_uv: raise ValueError("Error")
    if config.only_optical: raise ValueError("Error")
    if config.only_mir: raise ValueError("Error")
    if config.only_fir: raise ValueError("Error")
    if config.only_submm_microwave: raise ValueError("Error")
    regimes = ["nir"]

# Only MIR
elif config.only_mir:

    if config.no_mir: raise ValueError("Error")
    if config.only_uv: raise ValueError("Error")
    if config.only_optical: raise ValueError("Error")
    if config.only_nir: raise ValueError("Error")
    if config.only_fir: raise ValueError("Error")
    if config.only_submm_microwave: raise ValueError("Error")
    regimes = ["mir"]

# Only FIR
elif config.only_fir:

    if config.no_fir: raise ValueError("Error")
    if config.only_uv: raise ValueError("Error")
    if config.only_optical: raise ValueError("Error")
    if config.only_nir: raise ValueError("Error")
    if config.only_mir: raise ValueError("Error")
    if config.only_submm_microwave: raise ValueError("Error")
    regimes = ["fir"]

# Only submm/microwave
elif config.only_submm:

    if config.no_submm: raise ValueError("Error")
    if config.only_uv: raise ValueError("Error")
    if config.only_optical: raise ValueError("Error")
    if config.only_nir: raise ValueError("Error")
    if config.only_mir: raise ValueError("Error")
    if config.only_fir: raise ValueError("Error")
    regimes = ["submm-microwave"]

# Regimes
else: regimes = config.regimes[:]

# Ignore certain regimes?
if config.no_uv: regimes = sequences.removed_item(regimes, "uv")
if config.no_optical: regimes = sequences.removed_item(regimes, "optical")
if config.no_nir: regimes = sequences.removed_item(regimes, "nir")
if config.no_mir: regimes = sequences.removed_item(regimes, "mir")
if config.no_fir: regimes = sequences.removed_item(regimes, "fir")
if config.no_submm_microwave: regimes = sequences.removed_item(regimes, "submm-microwave")

# Check number of regimes
if len(regimes) == 0: raise ValueError("No regimes")

# -----------------------------------------------------------------

# UV weight
if "uv" in regimes: uv_weight = config.uv
else: uv_weight = 0.

# Optical weight
if "optical" in regimes: optical_weight = config.optical
else: optical_weight = 0.

# NIR weight
if "nir" in regimes: nir_weight = config.nir
else: nir_weight = 0.

# MIR weight
if "mir" in regimes: mir_weight = config.mir
else: mir_weight = 0.

# FIR weight
if "fir" in regimes: fir_weight = config.fir
else: fir_weight = 0.

# Submm weight
if "submm-microwave" in regimes: submm_microwave_weight = config.submm_microwave
else: submm_microwave_weight = 0.

# -----------------------------------------------------------------

print("UV", uv_weight)
print("Optical", optical_weight)
print("NIR", nir_weight)
print("MIR", mir_weight)
print("FIR", fir_weight)
print("Submm/Microwave", submm_microwave_weight)

# -----------------------------------------------------------------

# Inform the user
log.info("Calculating the weight to give to each band ...")

# Get the weights
weights = calculate_weights_filters(config.filters, uv=uv_weight, optical=optical_weight, nir=nir_weight, mir=mir_weight, fir=fir_weight, submm_microwave=submm_microwave_weight)

# Add to weights table
for fltr in weights: table.add_point(fltr, weights[fltr])

# -----------------------------------------------------------------

print(table)

table.saveto("weights.dat")

# -----------------------------------------------------------------
