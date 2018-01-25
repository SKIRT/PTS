#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.make_wavelength_grids (Re)make the wavelength grids for a fitting run.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.initialization.galaxy import create_basic_wavelength_grids, create_refined_wavelength_grids, create_highres_wavelength_grids
from pts.modeling.config.initialize_fit import default_npoints_range_basic, default_npoints_range_refined, default_npoints_range_highres
from pts.modeling.config.initialize_fit import default_ngrids_basic, default_ngrids_refined, default_ngrids_highres, default_wavelength_range

# -----------------------------------------------------------------

# Set the modeling path
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run name
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run for which to (re)make the wavelength grids", choices=runs.names)

# Options
definition.add_flag("basic", "make basic wavelength grids", True)
definition.add_flag("refined", "make refined wavelength grids", True)
definition.add_flag("highres", "make high-resolution wavelength grids", True)

# Settings for the wavelength grid generation
definition.add_optional("npoints_range_basic", "integer_range", "range of the basic wavelength grid size", default_npoints_range_basic, convert_default=True)
definition.add_optional("npoints_range_refined", "integer_range", "range of the refined wavelength grid size", default_npoints_range_refined, convert_default=True)
definition.add_optional("npoints_range_highres", "integer_range", "range of the high-resolution wavelength grid size", default_npoints_range_highres, convert_default=True)
definition.add_optional("ngrids_basic", "integer", "number of basic wavelength grids to generate", default_ngrids_basic)
definition.add_optional("ngrids_refined", "integer", "number of refined wavelength grids to generate", default_ngrids_refined)
definition.add_optional("ngrids_highres", "integer", "number of high-resolution wavelength grids to generate", default_ngrids_highres)
definition.add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)
definition.add_optional("range", "quantity_range", "range of wavelengths", default_wavelength_range, convert_default=True)

# Create the configuration
config = parse_arguments("make_wavelength_grids", definition, "(Re)make the wavelength grids for a fitting run")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.run)

# -----------------------------------------------------------------

# Determine the path to the wavelength grids directory
grids_path = fitting_run.wavelength_grids_path

# -----------------------------------------------------------------

# Check whether there are already wavelength grids
if fs.is_directory(grids_path) and not fs.is_empty(grids_path):
    fs.backup_directory(grids_path)
    fs.clear_directory(grids_path)

# -----------------------------------------------------------------

# Basic wavelength grids
if config.basic: create_basic_wavelength_grids()

# Refined wavelength grids
if config.refined: create_refined_wavelength_grids()

# High-resolution wavelength grids
if config.highres: create_highres_wavelength_grids()

# -----------------------------------------------------------------
