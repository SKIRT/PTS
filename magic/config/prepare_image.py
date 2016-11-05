#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.config.extract_sources import definition as extraction_definition
from pts.magic.config.subtract_sky import definition as subtraction_definition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add input section
definition.add_section("input", "the names of the input files")
definition.sections["input"].add_required("image", "image_path", "name of the input image")
definition.sections["input"].add_optional("galaxy_region", "file_path", "file with galaxy regions", "galaxies.reg")
definition.sections["input"].add_optional("star_region", "file_path", "file with star regions", "stars.reg")
definition.sections["input"].add_optional("saturation_region", "file_path", "file with regions for saturated stars", "saturation.reg")
definition.sections["input"].add_optional("other_region", "file_path", "file with regions for other contaminating sources", "other_sources.reg")
definition.sections["input"].add_optional("segments", "file_path", "image with segmentation maps (as planes 'galaxies', 'stars' and 'other_sources')", "segments.fits")

definition.add_optional("error_frame_names", "string_list", "the names of error planes to be included in the final error map")

definition.add_flag("write_steps", "write the results of intermediate steps")

# Sky subtraction
definition.add_flag("write_sky_apertures", "write sky apertures and sky annulus")
definition.add_optional("sky_apertures_path", "directory_path", "the path to the directory where the aperture frames and annulus region should be written to")

definition.add_flag("calculate_calibration_uncertainties", "calculate calibration uncertainties", True)

# Source extraction
definition.add_flag("extract_sources", "extract sources", True)
definition.import_section("extraction", "star extraction options", extraction_definition)

# Galactic extinction
definition.add_flag("correct_for_extinction", "correct for galactic extinction", True)
definition.add_optional("attenuation", "real", "the galactic attenuation")

# Unit conversion
definition.add_flag("convert_unit", "convert unit", True)
definition.add_section("unit_conversion", "unit conversion")
definition.sections["unit_conversion"].add_optional("to_unit", "string", "target unit", "MJy/sr")

# Convolution
definition.add_flag("convolve", "convolve")
definition.add_section("convolution", "convolution")
definition.sections["convolution"].add_optional("kernel_path", "file_path", "kernel path")
definition.sections["convolution"].add_optional("kernel_fwhm", "real", "kernel FWHM")
definition.sections["convolution"].add_optional("remote", "string", "remote host")

# Rebinning
definition.add_flag("rebin", "rebin")
definition.add_section("rebinning", "rebinning")
definition.sections["rebinning"].add_optional("reference_path", "file_path", "reference FITS path")

# Sky subtraction
definition.add_flag("subtract_sky", "subtract sky")
definition.add_section("sky_subtraction", "sky subtraction")

# Uncertainties
definition.add_flag("set_uncertainties", "set uncertainties")
definition.add_section("uncertainties", "uncertainties")
definition.sections["uncertainties"].add_optional("calibration_error", "calibration_error", "calibration error")

# Cropping
definition.add_flag("crop", "crop")
definition.add_section("cropping", "cropping")
definition.sections["cropping"].add_optional("limits", "pixel_limits", "pixel limits")

# -----------------------------------------------------------------
