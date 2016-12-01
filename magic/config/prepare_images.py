#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.config.extract import definition as extraction_definition
from pts.magic.config.subtract_sky import definition as subtraction_definition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# The dataset or image
definition.add_positional_optional("dataset", "file_path", "name of the dataset file or image file")

# The maximum FWHM (minimum resolution that has to be retained in the convolved frames)
definition.add_optional("max_fwhm", "quantity", "maximum FWHM for convolution")

# The minimum pixel size (minimum resolution that has to be retained in the rebinned frames)
definition.add_optional("max_pixelscale", "quantity", "maximum pixelscale for rebinning")

# Add input section
#definition.add_section("input", "the names of the input files")
#definition.sections["input"].add_required("image", "image_path", "name of the input image")
#definition.sections["input"].add_optional("galaxy_region", "file_path", "file with galaxy regions", "galaxies.reg")
#definition.sections["input"].add_optional("star_region", "file_path", "file with star regions", "stars.reg")
#definition.sections["input"].add_optional("saturation_region", "file_path", "file with regions for saturated stars", "saturation.reg")
#definition.sections["input"].add_optional("other_region", "file_path", "file with regions for other contaminating sources", "other_sources.reg")
#definition.sections["input"].add_optional("segments", "file_path", "image with segmentation maps (as planes 'galaxies', 'stars' and 'other_sources')", "segments.fits")

# Number of parallel processes
definition.add_optional("nprocesses", "integer", "number of parallel processes for the preparation", 8)

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
definition.sections["unit_conversion"].add_optional("to_unit", "unit", "target unit", "MJy/sr")

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
definition.sections["rebinning"].add_flag("exact", "exact rebinning method", True)

# Sky subtraction
definition.add_flag("subtract_sky", "subtract sky")
#definition.add_section("sky_subtraction", "sky subtraction")
definition.import_section("sky_subtraction", "sky subtraction options", subtraction_definition)

# Error maps
definition.add_flag("create_errormaps", "create error maps")
#definition.add_section("errormaps", "error maps")
#definition.sections["errormaps"].add_optional("calibration_error", "calibration_error", "calibration error")

# Cropping
definition.add_flag("crop", "crop")
definition.add_section("cropping", "cropping")
definition.sections["cropping"].add_optional("limits", "pixel_limits", "pixel limits")

# Writing
definition.add_flag("write", "write")
definition.add_section("writing", "writing options")
definition.add_optional("dataset_path", "string", "path for the output dataset")
definition.add_optional("statistics_path", "string", "path for the statistics file")

# Visualization
#definition.add_flag("visualise", "make visualisations")
#definition.add_section("visualisation")
definition.add_optional("visualisation_path", "directory_path", "directory for saving visualisations")

# -----------------------------------------------------------------
