#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.prepare_poisson Prepare the poisson frames (TEMPORARY SCRIPT)

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.filter import Filter
from pts.magic.core.frame import Frame
from pts.magic.core.remote import RemoteFrame
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.misc.extinction import GalacticExtinction
from pts.magic.misc.kernels import AnianoKernels
from pts.magic.core.kernel import ConvolutionKernel

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Get configuration
setter = ArgumentConfigurationSetter("prepare_poisson")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting prepare_poisson ...")

# -----------------------------------------------------------------

modeling_path = fs.cwd()

data_path = fs.join(modeling_path, "data")
prep_path = fs.join(modeling_path, "prep")

galex_images_path = fs.join(data_path, "images", "GALEX")
sdss_images_path = fs.join(data_path, "images", "SDSS")

# -----------------------------------------------------------------

# Path to the reference image
reference_path = fs.join(prep_path, "Pacs red", "initialized.fits")

# Load the reference wcs
reference_wcs = CoordinateSystem.from_file(reference_path)

# Get center coordinate
center_coordinate = reference_wcs.coordinate_range[0]

# -----------------------------------------------------------------

# Create the galactic extinction calculator
extinction = GalacticExtinction(center_coordinate)

# -----------------------------------------------------------------

# The aniano kernels service
aniano = AnianoKernels()

# -----------------------------------------------------------------

remote_host_id = "nancy"

# -----------------------------------------------------------------

# Paths
paths = fs.files_in_path(galex_images_path, extension="fits", contains="poisson") + fs.files_in_path(sdss_images_path, extension="fits", contains="poisson")

# Loop over the GALEX poisson frames
for path in paths:

    # Image name
    name = fs.strip_extension(fs.name(path))

    # Get instrument and band
    galaxy_name, instrument, band, _ = name.split("_")

    # Create the filter
    fltr = Filter.from_instrument_and_band(instrument, band)

    # Load the frame
    poisson = Frame.from_file(path)

    # Get the attenuation
    attenuation = extinction.extinction_for_filter(fltr)

    # CORRECT FOR GALACTIC EXTINCTION
    poisson *= 10**(0.4 * attenuation)

    # CONVERT UNIT to MJy/sr
    poisson *= 1e-6
    poisson /= poisson.pixelarea.to("sr").value
    poisson.unit = "MJy/sr"

    # Make frame remote
    remote_frame = RemoteFrame.from_local(poisson, remote_host_id)

    # Get the kernel path for convolution from this filter to the Pacs red filter
    kernel_file_path = aniano.get_kernel_path(fltr, "Pacs red")
    kernel = ConvolutionKernel.from_file(kernel_file_path)

    # Prepare the kernel
    kernel.prepare_for(poisson)

    # CONVOLVE
    remote_frame.convolve(kernel, allow_huge=True)

    # REBIN
    remote_frame.rebin(reference_wcs)

    # SAVE
    prep_path_band = fs.join(prep_path, instrument + " " + band)
    result_path = fs.join(prep_path_band, "poisson.fits")
    remote_frame.save(result_path)

# -----------------------------------------------------------------
