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

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

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

galaxy_name = fs.name(config.arguments.path)

fit_best_images_path = fs.join(config.arguments.path, "fit", "best", "images")

simulated = dict()

for path, name in fs.files_in_path(fit_best_images_path, extension="fits", returns=["path", "name"]):
    
    if name == galaxy_name + "_earth_total": continue
    
    frame = Frame.from_file(path)

    if "2MASS" in str(frame.filter): continue
    if "I4" in str(frame.filter): continue
    if "W3" in str(frame.filter): continue

    if frame.filter is None: raise ValueError("No filter for " + name)

    simulated[str(frame.filter)] = frame

trunc_path = fs.join(config.arguments.path, "truncated")

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
