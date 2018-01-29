#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_wavelength_grid Check a wavelength grid of a fitting run.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.core.tools import sequences
from pts.core.units.parsing import parse_unit as u
from pts.core.basics.log import log
from pts.core.tools import formatting as fmt
from pts.core.simulation.wavelengthgrid import WavelengthGrid

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run name
if runs.empty: raise ValueError("No fitting runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run for which to (re)make the wavelength grids", choices=runs.names)

# The wavelength grid name
definition.add_positional_optional("name", "string", "name of the wavelength grid to check")
definition.add_optional("grid_path", "file_path", "path of the wavelength grid to check")
definition.add_flag("skirt", "input grid is in SKIRT format")

# Flags
definition.add_flag("show", "show the wavelength grid", False)

# Options
definition.add_optional("min_npoints", "positive_integer", "minimum number of points required in filter wavelength range", 8)
definition.add_optional("min_npoints_fwhm", "positive_integer", "minimum number of points required in FWHM filter wavelength range", 5)

# Create the configuration
config = parse_arguments("check_wavelength_grid", definition, "Check a wavelength grid of a fitting run")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.run)

# -----------------------------------------------------------------

# Check whether each grid has a 'SKIRT' version
for name in fitting_run.wavelength_grid_data_names:

    # OK?
    if name in fitting_run.wavelength_grid_names: continue

    # Load the wavelength grid
    grid = fitting_run.get_wavelength_grid_data_grid(name)

    # Determine the new path
    path = fitting_run.get_wavelength_grid_path(name)

    # Save
    grid.to_skirt_input(path)

# -----------------------------------------------------------------

# Grid path is given
if config.grid_path is not None:

    # Set wavelength grid
    grid_name = fs.strip_extension(fs.name(config.grid_path))

    # Load the grid
    if config.skirt: wavelength_grid = WavelengthGrid.from_skirt_input(config.grid_path)
    else: wavelength_grid = WavelengthGrid.from_file(config.grid_path)

# Grid name is given
else:

    # Set wavelength grid name
    grid_name = config.name

    # Load the wavelength grid
    if config.name not in fitting_run.wavelength_grid_names: raise ValueError("Wavelength grid '" + config.name + "' does not exist")
    wavelength_grid = fitting_run.get_wavelength_grid(config.name)

# -----------------------------------------------------------------

# # Set spectral convolution flag
spectral_convolution = "refined" in grid_name or "highres" in grid_name

# -----------------------------------------------------------------

# Get the wavelengths as an array
wavelength_unit = u("micron")
wavelengths = wavelength_grid.wavelengths(wavelength_unit, asarray=True)

# -----------------------------------------------------------------

def check_grid_convolution():

    """
    This function ...
    :return:
    """

    # Wavelengths used for each filter
    wavelengths_for_filters = OrderedDict()

    # Show which wavelengths are used to create filter frames
    if len(wavelengths_for_filters) > 0:
        print("Wavelengths used for filters:")
        print("")
        for fltr in wavelengths_for_filters:
            filter_name = str(fltr)
            wavelength_strings = [str(wavelength) for wavelength in wavelengths_for_filters[fltr]]
            print(" - " + filter_name + ": " + ", ".join(wavelength_strings))
        print("")

# -----------------------------------------------------------------

def check_grid_no_convolution():

    """
    This function ...
    :return:
    """

    # Keep track of the wavelengths that have already been used to
    used_wavelengths = defaultdict(list)

    # Loop over the fitting filters
    for fltr in fitting_run.fitting_filters:

        # Determine the filter wavelength
        filter_wavelength = fltr.wavelength.to(wavelength_unit).value

        # Get the index of the wavelength closest to that of the filter
        index = sequences.find_closest_index(wavelengths, filter_wavelength)

        # Get the actual wavelength
        wavelength = wavelengths[index] * wavelength_unit

        # Get the difference
        difference = abs(filter_wavelength - wavelengths[index])
        reldifference = difference / filter_wavelength

        # Check grid wavelength in FWHM
        in_fwhm = wavelength in fltr.fwhm_range

        print(fltr, filter_wavelength, wavelength, reldifference, in_fwhm)

        wavelength_micron = wavelength.to("micron").value
        if wavelength_micron in used_wavelengths:
            filters = used_wavelengths[wavelength_micron]
            filter_names = [str(f) for f in filters]
            # log.warning("The simulated flux for the wavelength '" + str(wavelength) + "' has already been used to create SED point(s) for the " + ", ".join(filter_names) + " filter(s)")

        # Add the filter for the wavelength
        used_wavelengths[wavelength_micron].append(fltr)

    #print(used_wavelengths)

    # Show which wavelengths are used to create filter frames
    if len(used_wavelengths) > 0:
        print("")
        print("Used wavelengths and corresponding filter(s):")
        print("")
        for wavelength_micron in used_wavelengths:
            filters = used_wavelengths[wavelength_micron]
            filter_names = [str(f) for f in filters]
            nfilters = len(filter_names)
            if nfilters == 1: print(" - " + str(wavelength_micron) + " micron: " + filter_names[0])
            else: print(" - " + str(wavelength_micron) + " micron: " + fmt.bold + ", ".join(filter_names) + fmt.reset)
        print("")

# -----------------------------------------------------------------

# Wavelength grid for use with spectral convolution
if spectral_convolution: check_grid_convolution()

# No spectral convolution
else: check_grid_no_convolution()

# -----------------------------------------------------------------

if config.show and config.grid_path is not None: log.warning("Cannot show")
elif config.show:
    path = fitting_run.get_wavelength_grid_plot_path(config.name)
    fs.open_file(path)

# -----------------------------------------------------------------
