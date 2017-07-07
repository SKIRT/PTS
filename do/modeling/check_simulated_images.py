#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_simulated_images

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from pts.core.simulation.wavelengthgrid import WavelengthGrid
from pts.magic.core.datacube import DataCube
from pts.magic.core.remote import RemoteDataCube
from pts.core.filter.broad import BroadBandFilter
from pts.core.remote.python import AttachedPythonSession
from pts.core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Do the filter convolution remotely
definition.add_optional("remote", "string", "remote host")

# Specify the number of processes for the remote filter convolution
definition.add_optional("nprocesses", "integer", "number of processes to use for the filter convolution calculation", 8)

# Get configuration
setter = InteractiveConfigurationSetter("check_simulated_images")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting check_simulated_images ...")

# -----------------------------------------------------------------

# Modeling path and galaxy name
modeling_path = config.path
galaxy_name = fs.name(modeling_path)


"""
fit_best_images_path = fs.join(modeling_path, "fit_before_new", "best", "images", "test")

simulated = dict()
for path, name in fs.files_in_path(fit_best_images_path, extension="fits", returns=["path", "name"]):
    if name == galaxy_name + "_earth_total": continue
    frame = Frame.from_file(path)
    if frame.filter is None: raise ValueError("No filter for " + name)
    simulated[str(frame.filter)] = frame
"""


trunc_path = fs.join(modeling_path, "truncated")
prep_path = fs.join(modeling_path, "prep")

observed = dict()

# Loop over the directories in the preparation directory
for path, name in fs.directories_in_path(prep_path, returns=["path", "name"]):

    # Skip Halpha
    if "Halpha" in name: continue

    # Determine the path to the prepared image
    result_path = fs.join(path, "result.fits")

    # Load the frame
    frame = Frame.from_file(result_path)
    if frame.filter is None: raise ValueError("No filter for " + name)

    observed[str(frame.filter)] = frame

# Sort the filter names on wavelength
sorted_filter_names = sorted(observed.keys(), key=lambda key: observed[key].filter.pivotwavelength())

sorted_filter_names = sorted_filter_names[0:6]

##### CAN BE REMOVED ########
"""
plt.figure()

x = []
y = []



for filter_name in sorted_filter_names:

    sum_simulated = simulated[filter_name].sum()
    sum_observed = observed[filter_name].sum()

    ratio = sum_simulated / sum_observed

    wavelength = simulated[filter_name].filter.pivotwavelength()

    x.append(wavelength)
    y.append(ratio)

plt.plot(x, y)

plt.xscale('log')
plt.yscale('log')

plt.show()
"""
##############################

# -----------------------------------------------------------------

# Load wavelength grid
wavelength_grid_path = fs.join(modeling_path, "fit_before_new", "in", "wavelengths_lowres.txt")
wavelength_grid = WavelengthGrid.from_skirt_input(wavelength_grid_path)
wavelengths = wavelength_grid.wavelengths(asarray=True) # list of wavelengths

# Load simulated datacube
datacube_path = fs.join(modeling_path, "fit_before_new", "best", "images", "M81_earth_total.fits")

# Local or remote
session = AttachedPythonSession.from_host_id(config.remote)
if config.remote is not None: datacube = RemoteDataCube.from_file(datacube_path, wavelength_grid, session)
else: datacube = DataCube.from_file(datacube_path, wavelength_grid)

x = []
y = []

# output directory
new_path = fs.join(modeling_path, "fit_before_new", "best", "images", "new")
fs.create_directory(new_path)

### NEW

# Convert datacube to flux (wavelength) density
new_unit = "W / (m2 * arcsec2 * micron)"
wavelength_unit = "micron"
datacube.to_wavelength_density(new_unit, wavelength_unit)

# Created observed images
filters = [BroadBandFilter(filter_name) for filter_name in sorted_filter_names]

# Do the filter convolution
frames = datacube.convolve_with_filters(filters, nprocesses=config.nprocesses)

# Put frames in dictionary
images = dict()
for filter_name, frame in zip(sorted_filter_names, frames): images[filter_name] = frame

# Loop over the frames, convert them to MJy/sr again
for filter_name in images:

    frame = images[filter_name]

    # Get the wavelength
    wavelength = frame.filter.pivot

    # Determine the conversion factor
    conversion_factor = 1.0

    # From W / (m2 * arcsec2 * micron) to W / (m2 * arcsec2 * Hz)
    conversion_factor *= (wavelength ** 2 / speed_of_light).to("micron/Hz").value

    # From W / (m2 * arcsec2 * Hz) to MJy / sr
    # conversion_factor *= (u("W/(m2 * arcsec2 * Hz)") / u("MJy/sr")).to("")
    conversion_factor *= 1e26 * 1e-6 * (u("sr") / u("arcsec2")).to("")

    # Multiply the frame with the conversion factor
    frame *= conversion_factor

    # Set the new unit
    frame.unit = "MJy/sr"

    #print(frame.unit)
    #print(frame.remote.get_python_string("str(" + frame.label + ".filter)"))
    #print(frame.filter)

# Save the frames
for filter_name in images:

    # The frame
    frame = images[filter_name]

    sum_simulated = frame.sum()
    sum_observed = observed[filter_name].sum()

    ratio = sum_simulated / sum_observed

    wavelength = observed[filter_name].filter.pivotwavelength()

    x.append(wavelength)
    y.append(ratio)

    frame_path = fs.join(new_path, filter_name + ".fits")

    # Save the frame
    frame.saveto(frame_path)

# Plotting difference in total flux
plt.plot(x, y)

plt.xscale('log')
plt.yscale('log')

plt.show()

# -----------------------------------------------------------------
