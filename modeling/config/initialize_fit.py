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
from pts.modeling.fitting.initialization.base import wavelength_regimes

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

# For fitting weights

# Wavelength regimes
definition.add_optional("regimes", "string_list", "wavelength regimes to use", default=wavelength_regimes, choices=wavelength_regimes)

# Weights for different regimes
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

# -----------------------------------------------------------------