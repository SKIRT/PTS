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
from pts.core.tools.stringify import tostr

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

    print("")
    print(fmt.underlined + fmt.blue + "Filters for convolution:" + fmt.reset)
    print("")

    # Loop over the fitting filters
    for fltr in fitting_run.fitting_filters:

        # Get the wavelength indices in the ranges
        indices_in_minmax = [i for i in range(len(wavelengths)) if wavelengths[i] in fltr.range.to("micron").value]
        indices_in_fwhm = [i for i in range(len(wavelengths)) if wavelengths[i] in fltr.fwhm_range.to("micron").value]

        # Get the number of wavelengths in the ranges
        nwavelengths_in_minmax = len(indices_in_minmax)
        nwavelengths_in_fwhm = len(indices_in_fwhm)

        # SHow checks
        if nwavelengths_in_minmax < config.min_npoints:
            #raise ValueError("Too few wavelengths within the filter wavelength range (" + str(fltr.min.to("micron").value) + " to " + str(fltr.max.to("micron").value) + " micron) for convolution (" + str(nwavelengths_in_minmax) + ")")
            print(fmt.red + " - " + str(fltr) + ": too few wavelengths within the filter wavelength range (" + tostr(fltr.min) + " to " + tostr(fltr.max) + ") for convolution (" + str(nwavelengths_in_minmax) + " instead of " + str(config.min_npoints) + ")" + fmt.reset)
        elif nwavelengths_in_fwhm < config.min_npoints_fwhm:
            #raise ValueError("Too few wavelengths within the filter FWHM wavelength range (" + str(fltr.fwhm_min.to("micron").value) + " to " + str(fltr.fwhm_max.to("micron").value) + " micron) for convolution (" + str(nwavelengths_in_fwhm) + ")")
            print(fmt.red + " - " + str(fltr) + ": too few wavelengths within the filter FWHM wavelength range (" + tostr(fltr.fwhm_min) + " to " + tostr(fltr.fwhm_max) + ") for convolution (" + str(nwavelengths_in_fwhm) + " instead of " + str(config.min_npoints_fwhm) + ")" + fmt.reset)
        else: print(fmt.green + " - " + str(fltr) + ": filter range is sampled well by the wavelengths")

        # Set the wavelengths in the range
        wavelengths_in_minmax = [wavelengths[index] for index in indices_in_minmax]

        # Set the used wavelengths for this filter
        wavelengths_for_filters[fltr] = wavelengths_in_minmax

    # Show which wavelengths are used to create filter frames
    if len(wavelengths_for_filters) > 0:
        print("")
        print(fmt.underlined + fmt.blue + "Wavelengths used for filters:" + fmt.reset)
        print("")
        for fltr in wavelengths_for_filters:
            filter_name = str(fltr)
            wavelength_strings = [str(wavelength) for wavelength in wavelengths_for_filters[fltr]]
            print(" - " + fmt.bold + filter_name + fmt.reset + ": " + ", ".join(wavelength_strings))
        print("")

# -----------------------------------------------------------------

#max_reldifference = 1e-6

# -----------------------------------------------------------------

def check_grid_no_convolution():

    """
    This function ...
    :return:
    """

    # Keep track of the wavelengths that have already been used to
    used_wavelengths = defaultdict(list)

    print("")
    print(fmt.underlined + fmt.blue + "Filters without convolution:" + fmt.reset)
    print("")

    # Loop over the fitting filters
    for fltr in fitting_run.fitting_filters:

        # Determine the filter wavelength
        filter_wavelength = fltr.wavelength.to(wavelength_unit).value

        # Get the index of the wavelength closest to that of the filter
        index = sequences.find_closest_index(wavelengths, filter_wavelength)

        # Get the actual wavelength
        wavelength = wavelengths[index] * wavelength_unit

        # Get the difference
        #difference = abs(filter_wavelength - wavelengths[index])
        #reldifference = difference / filter_wavelength

        # Check grid wavelength in FWHM
        in_fwhm = wavelength in fltr.fwhm_range

        # Check whether the relative difference is smaller than 1e-6
        #close = reldifference < max_reldifference

        # Check grid wavelength in inner range
        in_inner = wavelength in fltr.inner_range

        # Show
        if not in_fwhm: print(fmt.red + " - " + str(fltr) + ": wavelength (" + tostr(wavelength) + ") not in the FWHM range (" + tostr(fltr.fwhm_range) + ") of the filter" + fmt.reset)
        #elif not close: print(fmt.yellow + " - " + str(fltr) + ": wavelength closest to the filter (" + tostr(wavelength) + ") differs more than " + tostr(max_reldifference*100) + "%" + fmt.reset)
        elif not in_inner: print(fmt.yellow + " - " + str(fltr) + ": wavelength (" + tostr(wavelength) + ") not in the inner range (" + tostr(fltr.inner_range) + ") of the filter" + fmt.reset)
        else: print(fmt.green + " - " + str(fltr) + ": wavelength found close to the filter (" + tostr(wavelength) + ")" + fmt.reset)

        wavelength_micron = wavelength.to("micron").value
        if wavelength_micron in used_wavelengths:
            filters = used_wavelengths[wavelength_micron]
            filter_names = [str(f) for f in filters]
            # log.warning("The simulated flux for the wavelength '" + str(wavelength) + "' has already been used to create SED point(s) for the " + ", ".join(filter_names) + " filter(s)")

        # Add the filter for the wavelength
        used_wavelengths[wavelength_micron].append(fltr)

    # Show which wavelengths are used to create filter frames
    if len(used_wavelengths) > 0:
        print("")
        print(fmt.underlined + fmt.blue + "Used wavelengths and corresponding filter(s):" + fmt.reset)
        print("")
        for wavelength_micron in used_wavelengths:
            filters = used_wavelengths[wavelength_micron]
            filter_names = [str(f) for f in filters]
            nfilters = len(filter_names)
            if nfilters == 1: print(fmt.green + " - " + str(wavelength_micron) + " micron: " + filter_names[0] + fmt.reset)
            else: print(fmt.yellow + " - " + str(wavelength_micron) + " micron: " + fmt.bold + ", ".join(filter_names) + fmt.reset)
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
