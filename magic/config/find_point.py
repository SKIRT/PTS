#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.tools.catalogs import stellar_catalog_descriptions

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

definition.add_flag("use_frame_fwhm", "If possible, avoid the fitting procedure and use the FWHM defined by the frame", True)
definition.add_optional("input_path", "directory_path", "path to the input directory")
definition.add_optional("output_path", "directory_path", "path to the output directory")
definition.add_flag("track_record", "track record", False)
definition.add_flag("plot_track_record_if_exception", True)
definition.add_optional("manual_region", "file_path", "manual star region")
definition.add_flag("remove", "remove stars from the frame", True)
definition.add_flag("find_saturation", "find saturated stars", True)

# -----------------------------------------------------------------

# Settings for fetching the catalogs
definition.add_section("fetching", "fetching")

# THE CATALOGS
default_catalogs = ["2MASS"]
definition.sections["fetching"].add_optional("catalogs", "string_list", "catalogs for point sources", default_catalogs, stellar_catalog_descriptions)

# -----------------------------------------------------------------

definition.sections["fetching"].add_flag("use_statistics_file", "use statistics file")
definition.sections["fetching"].add_optional("statistics_path", "file_path", "statistics file path")

definition.sections["fetching"].add_flag("cross_reference_with_galaxies", "blabla", True)

definition.sections["fetching"].add_section("min_distance_from_galaxy", "minimum distance the star has to be seperated from the galaxy center to be positively identified as a star")
definition.sections["fetching"].sections["min_distance_from_galaxy"].add_optional("principal", "real", "in pixels", 20.0)
definition.sections["fetching"].sections["min_distance_from_galaxy"].add_optional("companion", "real", "in pixels", 15.0)
definition.sections["fetching"].sections["min_distance_from_galaxy"].add_optional("other", "real", "in pixels", 15.0)

definition.add_section("detection", "source detection")

definition.sections["detection"].add_optional("initial_radius", "real", "initial radius (in pixels)", 10.0)
definition.sections["detection"].add_optional("background_est_method", "string", "background estimation method", "polynomial")
definition.sections["detection"].add_flag("sigma_clip_background", "sigma clip background", True)

definition.sections["detection"].add_optional("detection_method", "string", "detection method", "peaks")
definition.sections["detection"].add_optional("minimum_pixels", "integer", "minimum pixels", 5)
definition.sections["detection"].add_optional("sigma_level", "real", "threshold sigmas", 2.0)
definition.sections["detection"].add_optional("peak_offset_tolerance", "real", "peak offset tolerance (in pixels)", 3.0)

definition.sections["detection"].add_optional("min_level", "negative_integer", "minimum level", -2)
definition.sections["detection"].add_optional("max_level", "positive_integer", "maximum level", 2)

definition.sections["detection"].add_optional("scale_factor", "real", "scale factor", 2.0)

definition.sections["detection"].add_optional("background_outer_factor", "real", "background outer factor", 1.5)

definition.sections["detection"].add_flag("always_subtract_background", "always subtract background")

definition.sections["detection"].add_optional("convolution_fwhm", "real", "perform convolution, define the FWHM (in pixels) (for detection_method: 'peaks'", 10.0)

definition.sections["detection"].add_section("debug", "debug")
definition.sections["detection"].sections["debug"].add_flag("zero_peaks_before", "zero peaks before")
definition.sections["detection"].sections["debug"].add_flag("zero_peaks_after", "zero peaks after")
definition.sections["detection"].sections["debug"].add_flag("zero_peaks", "zero peaks")
definition.sections["detection"].sections["debug"].add_flag("one_peak", "one peak")
definition.sections["detection"].sections["debug"].add_flag("more_peaks", "more peaks")
definition.sections["detection"].sections["debug"].add_flag("off_center", "off-center")

definition.add_section("fitting", "fitting")

definition.sections["fitting"].add_optional("use_center_or_peak", "string", "use center of peak for fitting", "peak")
definition.sections["fitting"].add_optional("model_names", "string_list", "model names to use for fitting", ["Gaussian", "Airy"], choices=["Gaussian", "Airy"])
definition.sections["fitting"].add_optional("minimum pixels", "integer", "minimum number of pixels", 5)
definition.sections["fitting"].add_optional("max_model_offset", "real", "maximum model offset", 3.0)
definition.sections["fitting"].add_optional("zoom_factor", "real", "zoom factor", 2.0)
definition.sections["fitting"].add_optional("background_est_method", "string", "background estimation method", "polynomial")
definition.sections["fitting"].add_flag("sigma_clip_background", "sigma clip background", True)

definition.sections["fitting"].add_flag("sigma_clip_fwhms", "sigma clip FWHMs", True)
definition.sections["fitting"].add_optional("fwhm_sigma_level", "real", "FWHM sigma level (clipping)", 3.0)

definition.sections["fitting"].add_optional("upsample_factor", "real", "upsample factor", 1.0)

definition.sections["fitting"].add_section("debug", "debug")
definition.sections["fitting"].sections["debug"].add_flag("model_offset", "model offset")
definition.sections["fitting"].sections["debug"].add_flag("success", "success")

definition.sections["fitting"].add_flag("fit_if_undetected", "fit if undetected (no source (peak) found)")

definition.add_optional("source_psf_sigma_level", "real", "source PSF sigma level", 4.0)
definition.add_optional("source_outer_factor", "real", "source outer factor", 1.6)

definition.add_section("saturation", "saturated stars")

# NEW
definition.sections["saturation"].add_flag("deblend", "apply deblending", False)
definition.sections["saturation"].add_section("deblending", "deblending options")
definition.sections["saturation"].sections["deblending"].add_optional("min_npixels", "positive_integer", "minimum npixels for detection", 1)
definition.sections["saturation"].sections["deblending"].add_optional("contrast", "real", "contrast value", 0.001)
definition.sections["saturation"].sections["deblending"].add_optional("mode", "string", "mode", "exponential", ["exponential", "linear"])
definition.sections["saturation"].sections["deblending"].add_optional("nlevels", "positive_integer", "number of deblending levels", 10)

definition.sections["saturation"].add_flag("only_brightest", "only brightest")
definition.sections["saturation"].add_optional("brightest_method", "string", "brightest method", "percentage", choices=["percentage", "sigma clipping"])
definition.sections["saturation"].add_optional("brightest_level", "real", "for 'percentage': a percentage, for 'sigma clipping': a sigma level", 10.)

definition.sections["saturation"].add_optional("sigmas", "real", "sigmas for segmentation", 15.0)
definition.sections["saturation"].add_optional("background_outer_factor", "real", "background outer factor", 1.2)
definition.sections["saturation"].add_flag("always_subtract_background", "always subtract background", True)
definition.sections["saturation"].add_optional("background_est_method", "string", "background estimation method", "polynomial")
definition.sections["saturation"].add_flag("sigma_clip_background", "sigma clip background", True)
definition.sections["saturation"].add_optional("sigma_level", "real", "sigma level for clipping", 5.0)
definition.sections["saturation"].add_optional("expansion_factor", "real", "expansion factor", 1.5)

definition.sections["saturation"].add_optional("min_pixels", "integer", "minimum connected pixels", 5)
definition.sections["saturation"].add_section("kernel", "kernel")
definition.sections["saturation"].sections["kernel"].add_optional("fwhm", "real", "FWHM", 3.0)
definition.sections["saturation"].sections["kernel"].add_optional("cutoff_level", "real", "cutoff level (in sigmas)", 4.0)

definition.sections["saturation"].add_flag("expand", "expand", True)
definition.sections["saturation"].add_optional("max_expansion_level", "positive_integer", "maximum expansion level", 7)

definition.sections["saturation"].add_flag("allow_overlap", "do not normally allow overlap between the center segment and the background mask of the source")

definition.sections["saturation"].add_optional("interpolation_method", "string", "interpolation method for removing the saturation", "local_mean")
definition.sections["saturation"].add_flag("sigma_cip", "sigma clip", True)
definition.sections["saturation"].add_flag("no_sigma_clip_on_galaxy", "no sigma clipping on galaxy")
definition.sections["saturation"].add_flag("polynomial_on_galaxy", "polynomial on galaxy", True)

definition.sections["saturation"].add_section("debug", "debug")

definition.sections["saturation"].sections["debug"].add_flag("no_segment_before", "no segment before")
definition.sections["saturation"].sections["debug"].add_flag("no_segment_after", "no segment after")
definition.sections["saturation"].sections["debug"].add_flag("no_segment", "no segment")
definition.sections["saturation"].sections["debug"].add_flag("expand", "expand")
definition.sections["saturation"].sections["debug"].add_flag("success", "success")
definition.sections["saturation"].sections["debug"].add_flag("dilated", "dilated")
definition.sections["saturation"].sections["debug"].add_flag("user_expansion", "user expansion")
definition.sections["saturation"].sections["debug"].add_flag("overlap_before", "overlap before")
definition.sections["saturation"].sections["debug"].add_flag("overlap_after", "overlap after")

definition.sections["saturation"].add_flag("dilate", "dilate", True)
definition.sections["saturation"].add_optional("dilation_factor", "real", "dilation factor", 1.4)
definition.sections["saturation"].add_optional("iterations", "positive_integer", "iterations", 5)
definition.sections["saturation"].add_optional("connectiviy", "positive_integer", "connectivity", 2)

definition.sections["saturation"].add_flag("user_expansion", "user expansion")
definition.sections["saturation"].add_optional("user_expansion_factor", "real", "user expansion factor")

definition.sections["saturation"].add_flag("remove_if_not_fitted", "remove if not fitted", True)

definition.sections["saturation"].add_flag("remove_if_undetected", "remove if undetected", False)

definition.sections["saturation"].add_flag("remove_appendages", "remove appendages from overlapping mask", True)

definition.sections["saturation"].add_flag("remove_foreground", "remove foreground stars", True)

definition.sections["saturation"].add_flag("check_centroid", "check centroid")
definition.sections["saturation"].add_optional("max_centroid_offset", "real", "max centroid offset", 10.0)
definition.sections["saturation"].add_optional("max_centroid_ellipticity", "real", "max centroid ellipticity", 0.3)

definition.sections["saturation"].add_section("apertures", "apertures")
definition.sections["saturation"].sections["apertures"].add_optional("sigma_level", "real", "approximate isophotal extent", 4.0)
definition.sections["saturation"].sections["apertures"].add_optional("max_ellipticity", "real", "maximum ellipticity", 0.1)
definition.sections["saturation"].sections["apertures"].add_optional("max_offset", "real", "Maximal offset between the aperture center and star position (in number of pixels) (None=no limit)", 10.0)

definition.sections["saturation"].add_flag("remove_apertures", "remove apertures")
definition.sections["saturation"].add_section("aperture_removal", "aperture_removal")
definition.sections["saturation"].sections["aperture_removal"].add_optional("expansion_factor", "real", "expansion factor", 1.0)
definition.sections["saturation"].sections["aperture_removal"].add_optional("background_outer_factor", "real", "background outer factor", 1.2)
definition.sections["saturation"].sections["aperture_removal"].add_flag("no_sigma_clip_on_galaxy", "no sigma clip on galaxy", True)
definition.sections["saturation"].sections["aperture_removal"].add_flag("sigma_clip", "sigma clip", True)
definition.sections["saturation"].sections["aperture_removal"].add_flag("polynomial_on_galaxy", "polynomial on galaxy", True)
definition.sections["saturation"].sections["aperture_removal"].add_optional("interpolation_method", "string", "interpolation method", "local_mean")

definition.sections["saturation"].add_flag("second_segmentation", "second segmentation", False)
definition.sections["saturation"].add_optional("second_sigma_level", "real", "second sigma level", 1.2)

definition.add_section("region", "region")
definition.sections["region"].add_optional("sigma_level", "real", "sigma level", 5.0)

definition.add_section("fwhm", "calculation of the default FWHM")
definition.sections["fwhm"].add_optional("measure", "string", "measure", "mean", choices=["mean", "max", "median"])
definition.sections["fwhm"].add_optional("scale_factor", "real", "scale factor", 1.0)

# Flags
definition.add_flag("weak", "only do weak search")

# -----------------------------------------------------------------
