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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.initialization.galaxy import create_basic_wavelength_grids, create_refined_wavelength_grids, create_highres_wavelength_grids
from pts.modeling.config.initialize_fit import default_npoints_range_basic, default_npoints_range_refined, default_npoints_range_highres
from pts.modeling.config.initialize_fit import default_ngrids_basic, default_ngrids_refined, default_ngrids_highres, default_wavelength_range
from pts.core.prep.wavelengthgrids import WavelengthGridsTable
from pts.core.tools.stringify import tostr

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

# Backup existing wavelength grids
definition.add_flag("backup", "backup existing wavelength grids", True)

# Create the configuration
config = parse_arguments("make_wavelength_grids", definition, "(Re)make the wavelength grids for a fitting run")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.run)

#print("wavelength range of fitting filters: " + tostr(fitting_run.fitting_wavelength_range))
#print("wavelength range of fitting filters (min and max): " + tostr(fitting_run.absolute_fitting_wavelength_range))

#print("filters: " + tostr(fitting_run.normalization_filters))
#print("wavelengths: " + tostr(fitting_run.normalization_wavelengths))
#print("center wavelengths: " + tostr(fitting_run.normalization_center_wavelengths))
#print("effective wavelengths: " + tostr(fitting_run.normalization_effective_wavelengths))
#print("pivot wavelengths: " + tostr(fitting_run.normalization_pivot_wavelengths))
#print("mean wavelengths: " + tostr(fitting_run.normalization_mean_wavelengths))
#print("peak wavelengths: " + tostr(fitting_run.normalization_peak_wavelengths))

# Show the fitting filters
print("")
print("FITTING FILTERS:")
print("")
for fltr in fitting_run.fitting_filters: print(" - " + tostr(fltr))
print("")

# -----------------------------------------------------------------

# Determine the path to the wavelength grids directory
grids_path = fitting_run.wavelength_grids_path

# -----------------------------------------------------------------

# Check whether there are already wavelength grids
if fs.is_directory(grids_path) and not fs.is_empty(grids_path):
    if config.backup: fs.backup_directory(grids_path)
    fs.clear_directory(grids_path)

# -----------------------------------------------------------------

# Create a wavelength grids table
table = WavelengthGridsTable()

# -----------------------------------------------------------------

# Get fixed wavelengths
fixed_wavelengths = fitting_run.normalization_wavelengths

# -----------------------------------------------------------------

model_definition = fitting_run.model_definition
template_sed_old = model_definition.old_sed
template_sed_young = model_definition.young_sed
template_sed_ionizing = model_definition.ionizing_sed

# -----------------------------------------------------------------

# Make template SEDs (for plotting)
seds = OrderedDict()
seds["old"] = template_sed_old
seds["young"] = template_sed_young
seds["ionizing"] = template_sed_ionizing

# -----------------------------------------------------------------

# Basic wavelength grids
if config.basic:

    # Get the list of the different npoints
    basic_npoints_list = config.npoints_range_basic.linear(config.ngrids_basic)

    # Set paths
    basic_grid_paths = OrderedDict()
    for npoints in basic_npoints_list:

        # Determine path
        dirname = "basic_" + str(npoints)
        path = fs.join(grids_path, dirname)
        if fs.is_directory(path): fs.clear_directory(path)
        else: fs.create_directory(path)

        # Set path
        basic_grid_paths[npoints] = path

    # Generate the grids
    basic_grids = create_basic_wavelength_grids(config.ngrids_basic, config.npoints_range_basic, config.range,
                                               filters=fitting_run.fitting_filters, fixed=fixed_wavelengths,
                                               plot_seds=seds, table=table, out_paths=basic_grid_paths, plot_paths=basic_grid_paths)

# -----------------------------------------------------------------

# Refined wavelength grids
if config.refined:

    # Get the list of the different npoints
    refined_npoints_list = config.npoints_range_refined.linear(config.ngrids_refined)

    # Set paths
    refined_grid_paths = OrderedDict()
    for npoints in refined_npoints_list:

        # Determine path
        dirname = "refined_" + str(npoints)
        path = fs.join(grids_path, dirname)
        if fs.is_directory(path): fs.clear_directory(path)
        else: fs.create_directory(path)

        # Set path
        refined_grid_paths[npoints] = path

    # Generate the grids
    refined_grids = create_refined_wavelength_grids(config.ngrids_refined, config.npoints_range_refined, config.range,
                                                   filters=fitting_run.fitting_filters, fixed=fixed_wavelengths,
                                                   plot_seds=seds, table=table, out_paths=refined_grid_paths, plot_paths=refined_grid_paths)

# -----------------------------------------------------------------

# High-resolution wavelength grids
if config.highres:

    # Get the list of the different npoints
    highres_npoints_list = config.npoints_range_highres.linear(config.ngrids_highres)

    # Set paths
    highres_grid_paths = OrderedDict()
    for npoints in highres_npoints_list:

        # Determine path
        dirname = "highres_" + str(npoints)
        path = fs.join(grids_path, dirname)
        if fs.is_directory(path): fs.clear_directory(path)
        else: fs.create_directory(path)

        # Set path
        highres_grid_paths[npoints] = path

    # Generate the grids
    highres_grids = create_highres_wavelength_grids(config.ngrids_highres, config.npoints_range_highres, config.range,
                                                   filters=fitting_run.fitting_filters, fixed=fixed_wavelengths,
                                                   plot_seds=seds, table=table, out_paths=highres_grid_paths, plot_paths=highres_grid_paths)

# -----------------------------------------------------------------

# Save the table
table_path = fitting_run.wavelength_grids_table_path
table.saveto(table_path)

# -----------------------------------------------------------------
