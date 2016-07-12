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
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.units import Unit
from astropy import constants

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.basics.configuration import Configuration
from pts.core.simulation.wavelengthgrid import WavelengthGrid
from pts.magic.core.image import Image
from pts.magic.core.datacube import DataCube
from pts.magic.tools import plotting
from pts.core.basics.filter import Filter

# -----------------------------------------------------------------

# Create the configuration
config = Configuration("check_simulated_images")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting check_simulated_images ...")

# -----------------------------------------------------------------

modeling_path = config.fixed["path"]

galaxy_name = fs.name(modeling_path)

fit_best_images_path = fs.join(modeling_path, "fit", "best", "images", "test")

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

    sum_simulated = np.sum(simulated[filter_name])
    sum_observed = np.sum(observed[filter_name])

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
wavelength_grid_path = fs.join(modeling_path, "fit", "in", "wavelengths_lowres.txt")
wavelength_grid = WavelengthGrid.from_skirt_input(wavelength_grid_path)
wavelengths = wavelength_grid.wavelengths(asarray=True) # list of wavelengths

# Load simulated datacube
datacube_path = fs.join(modeling_path, "fit", "best", "images", "M81_earth_total.fits")
datacube = DataCube.from_file(datacube_path, wavelength_grid)

x = []
y = []

new_path = fs.join(modeling_path, "fit", "best", "images", "new")
fs.create_directory(new_path)

### NEW

# Convert datacube to flux (wavelength) density
datacube.convert_to_fluxdensity("W / (m2 * arcsec2 * micron)")

# Pack datacube into a 3D array
fluxdensities = datacube.asarray()

###

# Loop over filters
for filter_name in sorted_filter_names:

    # Filter wavelength
    fltr = Filter.from_string(filter_name)
    filter_wavelength = fltr.pivotwavelength()

    print(filter_name, filter_wavelength)

    # Index of closest index to wavelength
    index = wavelength_grid.closest_wavelength_index(filter_wavelength)

    # Get frame of closest wavelength to IRAC
    label = "frame"+str(index)
    frame = image.frames[label]

    # Divide by wavelength
    wavelength = wavelength_grid[index]
    #frame /= wavelength.to("micron").value

    print(filter_name, wavelength)

    #frame *= 1./wavelength.to("micron").value
    frame /= wavelength.to("micron").value

    ## -- FILTER CONVOLUTION HERE --

    print(wavelengths)
    data = fltr.convolve(wavelengths, fluxdensities)
    frame = Frame(data)

    # Set the unit of the frame
    frame.unit = "W/(m2 * arcsec2 * micron)"

    ##


    # The speed of light
    speed_of_light = constants.c

    # Determine the conversion factor
    conversion_factor = 1.0
    # From W / (m2 * arcsec2 * micron) to W / (m2 * arcsec2 * Hz)
    conversion_factor *= (wavelength ** 2 / speed_of_light).to("micron/Hz").value
    # From W / (m2 * arcsec2 * Hz) to MJy / sr
    # conversion_factor *= (Unit("W/(m2 * arcsec2 * Hz)") / Unit("MJy/sr")).to("")
    conversion_factor *= 1e26 * 1e-6 * (Unit("sr") / Unit("arcsec2")).to("")

    frame *= conversion_factor

    # Plot the difference
    #plotting.plot_difference(frame, simulated[filter_name])

    sum_simulated = np.sum(frame)
    sum_observed = np.sum(observed[filter_name])

    ratio = sum_simulated / sum_observed

    wavelength = simulated[filter_name].filter.pivotwavelength()

    x.append(wavelength)
    y.append(ratio)

    frame_path = fs.join(new_path, filter_name + ".fits")

    frame.save(frame_path)

    exit()

plt.plot(x, y)

plt.xscale('log')
plt.yscale('log')

plt.show()

# -----------------------------------------------------------------
