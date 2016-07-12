#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration("explore", log_path="log")

# Positional optional parameter
config.add_positional_optional("generation_method", str, "the model generation method ('grid', 'instinctive', 'genetic')", "genetic", ["genetic", "grid", "instinctive"])

# Optional parameters
config.add_optional("remote", str, "the remote host on which to run the parameters exploration", "nancy")
config.add_optional("simulations", int, "the number of simulations to launch in one batch/generation", 100)
config.add_optional("young", "real_range", "the range of the FUV luminosity of the young stellar population", (0.0, 4.e16))
config.add_optional("ionizing", "real_range", "the range of the FUV luminosity of the ionizing stellar population", (0.0, 5.e10))
config.add_optional("dust", "quantity_range", "the range of the dust mass", (0.5e7, 3.e7))

# Flags
config.add_flag("relative", "whether the range values are relative to the best (or initial) parameter value")
config.add_flag("young_log", "use logarithmic spacing of the young stellar luminosity values")
config.add_flag("ionizing_log", "use logarithmic spacing of the ionizing stellar luminosity values")
config.add_flag("dust_log", "use logarithmic spacing of the dust mass values")
config.add_flag("visualise", "make visualisations")
config.add_flag("dry", "dry-run (don't actually launch simulations)")
config.add_flag("refine_wavelengths", "increase the resolution of the wavelength grid for the new batch of simulations")
config.add_flag("refine_dust", "increase the resolution of the dust cell grid for the new batch of simulations")

# -----------------------------------------------------------------
