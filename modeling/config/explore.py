#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.core.component import get_free_parameter_labels

# -----------------------------------------------------------------

# Determine the path to the free parameters table
free_parameters_path = fs.join(fs.cwd(), "fit", "free_parameters.txt")

# Get the labels of the free parameters
free_parameter_labels = get_free_parameter_labels(free_parameters_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Positional optional parameter
definition.add_positional_optional("generation_method", "string", "the model generation method ('grid', 'instinctive', 'genetic')", "genetic", ["genetic", "grid", "instinctive"])

# Optional parameters
definition.add_optional("remote", "string", "the remote host on which to run the parameters exploration", "nancy")
definition.add_optional("simulations", "integer", "the number of simulations to launch in one batch/generation", 100)
definition.add_optional("young", "real_range", "the range of the FUV luminosity of the young stellar population", "0.0:4.e16", convert_default=True)
definition.add_optional("ionizing", "real_range", "the range of the FUV luminosity of the ionizing stellar population", "0.0:5.e10", convert_default=True)
definition.add_optional("dust", "quantity_range", "the range of the dust mass", "0.5e7:3.e7", convert_default=True)

# Flags
definition.add_flag("relative", "whether the range values are relative to the best (or initial) parameter value")
definition.add_flag("young_log", "use logarithmic spacing of the young stellar luminosity values")
definition.add_flag("ionizing_log", "use logarithmic spacing of the ionizing stellar luminosity values")
definition.add_flag("dust_log", "use logarithmic spacing of the dust mass values")
definition.add_flag("visualise", "make visualisations")
definition.add_flag("dry", "dry-run (don't actually launch simulations)")
definition.add_flag("refine_wavelengths", "increase the resolution of the wavelength grid for the new batch of simulations")
definition.add_flag("refine_dust", "increase the resolution of the dust cell grid for the new batch of simulations")

# -----------------------------------------------------------------
