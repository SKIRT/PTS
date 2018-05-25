#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
run_names = environment.analysis_runs.names

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

origins = ["model", "fitting_run"]
dust_grid_types = ["cartesian", "bintree", "octtree"]
default_dust_grid_type = "bintree"

# -----------------------------------------------------------------

definition.add_required("origin", "string", "origin of the analysis model", choices=origins)

# -----------------------------------------------------------------

# Give the analysis run a custom name
definition.add_optional("name", "string", "name for the analysis run") #forbidden=run_names) # don't forbid because --overwrite flag exists

# -----------------------------------------------------------------

# Wavelength grid from fitting run
definition.add_optional("wavelength_grid", "string", "wavelength grid from fitting run (only if origin is fitting run)")

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints", "positive_integer", "range of the wavelength grid size", 450)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)
definition.sections["wg"].add_optional("range", "quantity_range", "range of wavelengths", "0.1 micron > 2000 micron", convert_default=True)

# -----------------------------------------------------------------

# Representation from model building
definition.add_optional("representation", "string", "representation name")

# Dust grid properties
definition.add_section("dg", "settings for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale", "real", "number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 0.5)
definition.sections["dg"].add_optional("bintree_min_level", "integer", "minimum depth level for binary trees", 9)
definition.sections["dg"].add_optional("octtree_min_level", "integer", "minimum depth level for octrees", 3)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "maximum mass fraction in each cell", 0.5e-6)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 10.)

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 25.)
definition.add_optional("radial_factor", "real", "factor with which to multiply the radial extent of the projections", 1.7)

# -----------------------------------------------------------------

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", True)

# -----------------------------------------------------------------

# Whether model has to be adapted
definition.add_flag("adapt", "adapt the parameters of the chosen model (from model suite origin, not from fitting)", True)

# ADVANCED: specify the model name on the command line
definition.add_optional("model_name", "string", "name of the model (from the model suite) to use (only specify when origin is 'model')")

# -----------------------------------------------------------------

# ADVANCED: force overwrite
definition.add_flag("overwrite", "overwrite a possibly existing analysis run with this name (use with care!)", False)

# -----------------------------------------------------------------
