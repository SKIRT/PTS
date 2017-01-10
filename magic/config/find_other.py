#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Input and output
definition.add_optional("input", "directory_path", "input directory")
definition.add_optional("output", "directory_path", "output directory")

# Flags
definition.add_flag("remove", "remove ...", True)
definition.add_flag("find_apertures", "find apertures", False)
definition.add_flag("remove_apertures", "remove apertures", False)
#definition.add_flag("classify", "classify", True)

definition.add_section("detection", "source detection")

definition.sections["detection"].add_optional("method", "string", "detection method", "segmentation", choices=["segmentation", "peaks", "sextractor"])

definition.sections["detection"].add_section("peaks", "options for peak detection")

definition.sections["detection"].add_section("segmentation", "options for segmentation")

definition.sections["detection"].sections["segmentation"].add_optional("clipping_sigma_level", "real", "sigma level for sigma-clipping", 3.0)
definition.sections["detection"].sections["segmentation"].add_optional("sigma_level", "real", "sigma level of detection", 3.0)

definition.sections["detection"].add_section("sextractor", "options for SExtractor")

definition.sections["detection"].sections["sextractor"].add_optional("zero_point", "real", "zero point in mag arcsec^-2", 20.0)
definition.sections["detection"].sections["sextractor"].add_optional("gain", "real", "gain in [e-/ADU]", 4.0)
definition.sections["detection"].sections["sextractor"].add_optional("pixelscale", "real", "pixelscale in [arcsec/pix]", 1.0)
definition.sections["detection"].sections["sextractor"].add_optional("fwhm", "real", "FWHM of the PSF", 1.0)
definition.sections["detection"].sections["sextractor"].add_optional("input_file", "string", "sextractor input file", "default.sex")
definition.sections["detection"].sections["sextractor"].add_flag("remove_temp", "remove the temporary directory")

definition.sections["detection"].add_section("apertures", "apertures")
definition.sections["detection"].sections["apertures"].add_optional("sigma_level", "real", "sigma level", 4.0)
  

definition.add_flag("dilate", "dilate", False)
definition.add_optional("dilation_factor", "real", "dilation factor", 1.4)

definition.add_section("apertures", "apertures")

definition.add_section("aperture_removal", "aperture removal")

#definition.add_section("classification", "classification")

#definition.sections["classification"].add_section("fitting", "fitting")
#definition.sections["classification"].sections["fitting"].add_optional("use_center_or_peak", "string", "use center or peak", "peak")
#definition.sections["classification"].sections["fitting"].add_optional("model_names", "string_list", "model names", ["Gaussian", "Airy"])
#definition.sections["classification"].sections["fitting"].add_optional("initial_sigma", "real", "initial sigma")
#definition.sections["classification"].sections["fitting"].add_optional("minimum_pixels", "positive_integer", "minimum pixels", 5)
#definition.sections["classification"].sections["fitting"].add_optional("max_model_offset", "real", "maximum model offset", 3.0)
#definition.sections["classification"].sections["fitting"].add_optional("zoom_factor", "real", "zoom factor", 2.0)
#definition.sections["classification"].sections["fitting"].add_optional("background_est_method", "string", "background estimation method", "polynomial")
#definition.sections["classification"].sections["fitting"].add_flag("sigma_clip_background", "sigma clip background", True)
#definition.sections["classification"].sections["fitting"].add_flag("sigma_clip_fwhms", "sigma clip FWHMs", False)
#definition.sections["classification"].sections["fitting"].add_optional("fwhm_sigma_level", "real", "FWHM sigma level for clipping", 3.0)
#definition.sections["classification"].sections["fitting"].add_optional("upsample_factor", "real", "upsample factor", 1.0)
#definition.sections["classification"].sections["fitting"].add_section("debug", "debug")
#definition.sections["classification"].sections["fitting"].sections["debug"].add_flag("model_offset", "model offset")
#definition.sections["classification"].sections["fitting"].sections["debug"].add_flag("success", "success")
#definition.sections["classification"].sections["fitting"].add_flag("fit_if_undetected", "fit if undetected (no source (peak) found)")

# -----------------------------------------------------------------
