#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# FITTING RUN
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", runs.names)

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints_range", "integer_range", "the range of the wavelength grid size", "150>500", convert_default=True)
definition.sections["wg"].add_optional("ngrids", "integer", "the number of wavelength grids to generate", 10)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", False)
definition.sections["wg"].add_optional("min_wavelength", "quantity", "minimum wavelength")
definition.sections["wg"].add_optional("max_wavelength", "quantity", "maximum wavelength")

# Add optional arguments
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 2e5)
definition.add_flag("selfabsorption", "enable dust self-absorption")
definition.add_flag("transient_heating", "enable transient heating", True)

# -----------------------------------------------------------------
