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
from astropy.units import Unit
from astropy import constants

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.simulation.wavelengthgrid import WavelengthGrid
from pts.magic.core.datacube import DataCube
from pts.magic.core.remote import RemoteDataCube
from pts.magic.tools import plotting
from pts.core.basics.filter import Filter

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

definition.add_optional("remote", str, "remote host")

# Get configuration
setter = ArgumentConfigurationSetter("check_simulated_images")
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

modeling_path = config.path

galaxy_name = fs.name(modeling_path)

fit_best_images_path = fs.join(modeling_path, "fit_before_new", "best", "images", "test")

simulated = dict()

for path, name in fs.files_in_path(fit_best_images_path, extension="fits", returns=["path", "name"]):
    
    if name == galaxy_name + "_earth_total": continue
    
    frame = Frame.from_file(path)

    if "2MASS" in str(frame.filter): continue
    if "I4" in str(frame.filter): continue
    if "W3" in str(frame.filter): continue

    if frame.filter is None: raise ValueError("No filter for " + name)

    simulated[str(frame.filter)] = frame

trunc_path = fs.join(modeling_path, "truncated")

observed = dict()

for path, name in fs.files_in_path(trunc_path, extension="fits", returns=["path", "name"]):

    if name == "bulge" or name == "disk" or name == "model": continue

    frame = Frame.from_file(path)

    if "2MASS" in str(frame.filter) or "656_1" in str(frame.filter): continue
    if "I4" in str(frame.filter): continue
    if "W3" in str(frame.filter): continue

    if frame.filter is None: raise ValueError("No filter for " + name)

    observed[str(frame.filter)] = frame

plt.figure()

x = []
y = []

sorted_filter_names = sorted(simulated.keys(), key=lambda key: simulated[key].filter.pivotwavelength())

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

# -----------------------------------------------------------------

# Load wavelength grid
wavelength_grid_path = fs.join(modeling_path, "fit_before_new", "in", "wavelengths_lowres.txt")
wavelength_grid = WavelengthGrid.from_skirt_input(wavelength_grid_path)
wavelengths = wavelength_grid.wavelengths(asarray=True) # list of wavelengths

# Load simulated datacube
datacube_path = fs.join(modeling_path, "fit_before_new", "best", "images", "M81_earth_total.fits")

# Local or remote
if config.remote is not None: datacube = RemoteDataCube.from_file(datacube_path, wavelength_grid, config.remote)
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
filters = [Filter.from_string(filter_name) for filter_name in sorted_filter_names]

# Do the filter convolution
frames = datacube.convolve_with_filters(filters, parallel=True)

# Put frames in dictionary
images = dict()
for filter_name, frame in zip(sorted_filter_names, frames): images[filter_name] = frame

# Loop over the frames, convert them to MJy/sr again
for filter_name in images:

    frame = images[filter_name]

    # Get the wavelength
    wavelength = frame.filter.pivotwavelength() * Unit("micron")

    # Determine the conversion factor
    conversion_factor = 1.0

    # From W / (m2 * arcsec2 * micron) to W / (m2 * arcsec2 * Hz)
    conversion_factor *= (wavelength ** 2 / speed_of_light).to("micron/Hz").value

    # From W / (m2 * arcsec2 * Hz) to MJy / sr
    # conversion_factor *= (Unit("W/(m2 * arcsec2 * Hz)") / Unit("MJy/sr")).to("")
    conversion_factor *= 1e26 * 1e-6 * (Unit("sr") / Unit("arcsec2")).to("")

    # Multiply the frame with the conversion factor
    frame *= conversion_factor

    # Set the new unit
    frame.unit = "MJy/sr"

    print(frame.unit)
    print(frame.remote.get_python_string("str(" + frame.label + ".filter)"))
    print(frame.filter)

# Save the frames
for filter_name in images:

    # The frame
    frame = images[filter_name]

    sum_simulated = frame.sum()
    sum_observed = observed[filter_name].sum()

    ratio = sum_simulated / sum_observed

    wavelength = simulated[filter_name].filter.pivotwavelength()

    x.append(wavelength)
    y.append(ratio)

    frame_path = fs.join(new_path, filter_name + ".fits")

    # Save the frame
    frame.save(frame_path)

# Plotting difference in total flux
plt.plot(x, y)

plt.xscale('log')
plt.yscale('log')

plt.show()

# -----------------------------------------------------------------
