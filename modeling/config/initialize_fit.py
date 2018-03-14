#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.modeling.config.component import definition
from pts.modeling.config.calculate_weights import definition as calculate_weights_definition

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

default_npoints_range_basic = "50>250"
default_npoints_range_refined = "50>250"
default_npoints_range_highres = "150>400"

# -----------------------------------------------------------------

default_ngrids_basic = 6
default_ngrids_refined = 5
default_ngrids_highres = 4

# -----------------------------------------------------------------

default_wavelength_range = "0.02 micron > 2000 micron"
default_npackages = 2e5

# -----------------------------------------------------------------

definition = definition.copy()

# Fitting run
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", runs.names)

# -----------------------------------------------------------------

definition.add_flag("regenerate_wavelength_grids", "regenerate the wavelength grids", False)

# -----------------------------------------------------------------

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints_range_basic", "integer_range", "range of the basic wavelength grid size", default_npoints_range_basic, convert_default=True)
definition.sections["wg"].add_optional("npoints_range_refined", "integer_range", "range of the refined wavelength grid size", default_npoints_range_refined, convert_default=True)
definition.sections["wg"].add_optional("npoints_range_highres", "integer_range", "range of the high-resolution wavelength grid size", default_npoints_range_highres, convert_default=True)
definition.sections["wg"].add_optional("ngrids_basic", "integer", "number of basic wavelength grids to generate", default_ngrids_basic)
definition.sections["wg"].add_optional("ngrids_refined", "integer", "number of refined wavelength grids to generate", default_ngrids_refined)
definition.sections["wg"].add_optional("ngrids_highres", "integer", "number of high-resolution wavelength grids to generate", default_ngrids_highres)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)
definition.sections["wg"].add_optional("range", "quantity_range", "range of wavelengths", default_wavelength_range, convert_default=True)

# -----------------------------------------------------------------

# Other simulation settings
definition.add_optional("npackages", "real", "number of photon packages per wavelength", default_npackages)
definition.add_flag("selfabsorption", "enable dust self-absorption", False)
definition.add_flag("transient_heating", "enable transient heating", True)

# -----------------------------------------------------------------

# For calculating weights
definition.import_section("weighing", "wavelength filter weighing", calculate_weights_definition)

# -----------------------------------------------------------------
