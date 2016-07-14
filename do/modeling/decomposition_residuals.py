#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.decomposition_residuals Calculate the residuals between the 3.6 micron image and the simulated bulge and disk images.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit
from astropy.modeling.models import Gaussian2D

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.modeling.decomposition.decomposition import load_parameters
from pts.magic.core.source import Source
from pts.magic.tools import statistics, plotting, fitting
from pts.magic.basics.geometry import Ellipse
from pts.magic.basics.vector import Extent
from pts.core.basics.configuration import ConfigurationDefinition, ConfigurationReader

# -----------------------------------------------------------------

# Create the configuration
config = Configuration("decomposition_residuals")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting decomposition_residuals ...")

# -----------------------------------------------------------------

components_path = fs.join(config.arguments.path, "components")
truncation_path = fs.join(config.arguments.path, "truncated")
residuals_path = fs.join(config.arguments.path, "residuals")

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the IRAC 3.6 micron image ...")

# Determine the path to the truncated 3.6 micron image
path = fs.join(truncation_path, "IRAC I1.fits")
frame = Frame.from_file(path)

# Convert the frame to Jy/pix
conversion_factor = 1.0
conversion_factor *= 1e6

# Convert the 3.6 micron image from Jy / sr to Jy / pixel
pixelscale = frame.average_pixelscale
pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
conversion_factor /= pixel_factor
frame *= conversion_factor
frame.unit = "Jy"

frame.save(fs.join(residuals_path, "i1_jy.fits"))

# Inform the user
log.info("Loading the bulge image ...")

# Determine the path to the truncated bulge image
bulge_path = fs.join(truncation_path, "bulge.fits")
bulge = Frame.from_file(bulge_path)

# Inform the user
log.info("Loading the disk image ...")

# Determine the path to the truncated disk image
disk_path = fs.join(truncation_path, "disk.fits")
disk = Frame.from_file(disk_path)

# Inform the user
log.info("Loading the model image ...")

# Determine the path to the truncated model image
model_path = fs.join(truncation_path, "model.fits")
model = Frame.from_file(model_path)

# -----------------------------------------------------------------

# Calculate the bulge residual frame
bulge_residual = frame - bulge
bulge_residual_path = fs.join(residuals_path, "bulge_residual.fits")
bulge_residual.save(bulge_residual_path)

#bulge_residual2 = frame - (bulge * 1.3)
#bulge_residual2_path = fs.join(residuals_path, "bulge_residual_1,3.fits")
#bulge_residual2.save(bulge_residual2_path)

# Calculate the disk residual frame
disk_residual = frame - disk
disk_residual_path = fs.join(residuals_path, "disk_residual.fits")
disk_residual.save(disk_residual_path)

# Calculate the model residual frame
model_residual = frame - model
#model_residual = frame - (bulge*1.3)
model_residual = model_residual - disk
model_residual_path = fs.join(residuals_path, "model_residual.fits")
model_residual.save(model_residual_path)

# -----------------------------------------------------------------

exit()

# FWHM of all the images
fwhm = 11.18 * Unit("arcsec")
fwhm_pix = (fwhm / frame.average_pixelscale).to("pix").value
sigma = fwhm_pix * statistics.fwhm_to_sigma

# Get the center pixel of the galaxy
parameters_path = fs.join(components_path, "parameters.dat")
parameters = load_parameters(parameters_path)
center = parameters.center.to_pixel(frame.wcs)

# Create a source around the galaxy center
ellipse = Ellipse(center, 20.0*sigma)
source = Source.from_ellipse(model_residual, ellipse, 1.5)

source.estimate_background("polynomial")

source.plot()

position = source.center
model = source.subtracted.fit_model(position, "Gaussian")

rel_center = center - Extent(source.x_min, source.y_min)
rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
plotting.plot_peak_model(source.cutout, rel_center.x, rel_center.y, rel_model)

model_fwhm_pix = fitting.fwhm(model) * Unit("pix")
model_fwhm = (model_fwhm_pix * frame.average_pixelscale).to("arcsec")

print("Model FWHM: ", model_fwhm)

evaluated_model = source.cutout.evaluate_model(model)

all_residual = Frame(np.copy(model_residual))
all_residual[source.y_slice, source.x_slice] -= evaluated_model
all_residual.save(fs.join(residuals_path, "all_residual.fits"))

model = Gaussian2D(amplitude=0.0087509425805, x_mean=center.x, y_mean=center.y, x_stddev=sigma, y_stddev=sigma)
rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
plotting.plot_peak_model(source.cutout, rel_center.x, rel_center.y, rel_model)

evaluated_model = source.cutout.evaluate_model(model)

all_residual2 = Frame(np.copy(model_residual))
all_residual2[source.y_slice, source.x_slice] -= evaluated_model
all_residual2.save(fs.join(residuals_path, "all_residual2.fits"))

# -----------------------------------------------------------------
