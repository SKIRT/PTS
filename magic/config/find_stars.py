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

definition.add_flag("use_frame_fwhm", "If possible, avoid the fitting procedure and use the FWHM defined by the frame", True)

definition.add_optional("input_path", "directory_path", "path to the input directory")

definition.add_optional("output_path", "directory_path", "path to the output directory")

definition.add_flag("track_record", "track record", False)
definition.add_flag("plot_track_record_if_exception", True)

definition.add_optional("manual_region", "file_path", "manual star region")

definition.add_flag("remove", "remove stars from the frame", True)

definition.add_flag("find_saturation", "find saturated stars", True)

definition.add_section("fetching", "fetching")
definition.sections["fetching"].add_flag("use_catalog_file", "use catalog file")
definition.sections["fetching"].add_optional("catalog_path", "file_path", "catalog path")

definition.sections["fetching"].add_flag("use_statistics_file", "use statistics file")
definition.sections["fetching"].add_optional("statistics_path", "file_path", "statistics file path")

definition.sections["fetching"].add_flag("cross_reference_with_galaxies", "blabla", True)

definition.sections["fetching"].add_section("min_distance_from_galaxy", " minimum distance the star has to be seperated from the galaxy center to be positively identified as a star")
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

# Fitting
fitting:
{
  use_center_or_peak: "peak"
  model_names: ["Gaussian", "Airy"]
  initial_sigma: None
  minimum_pixels: 5
  max_model_offset: 3.0
  zoom_factor: 2.0
  background_est_method: "polynomial"
  sigma_clip_background: True

  # sigma-clipping
  sigma_clip_fwhms: True
  fwhm_sigma_level: 3.0

  # Upsample factor
  upsample_factor: 1.0

  # Debug
  debug:
  {
    model_offset: False
    success: False
  }
  
  # Fit if undetected (no source (peak) found)
  fit_if_undetected: False
}

source_psf_sigma_level: 4.0
source_outer_factor: 1.6

# Saturated stars
saturation:
{
  only_brightest: False
  brightest_method: "percentage"  # or: "sigma clipping"
  brightest_level: 10.  # for "percentage": a percentage, for "sigma_clip": a sigma-level
  
  # For segmentation:
  sigmas: 15.0
  background_outer_factor: 1.2
  always_subtract_background: True
  background_est_method: "polynomial"
  sigma_clip_background: True
  sigma_level: 5.0
  expansion_factor: 1.5
  
  # Minimum connected pixels (int)
  min_pixels: 5
  
  # Kernel
  kernel:
  {
    fwhm: 3.0
    cutoff_level: 4.0 # in sigmas
  }
  
  expand: True
  max_expansion_level: 7
  
  # Do not normally allow overlap between the center segment and the background mask of the source
  allow_overlap: False
  
  # For removing the saturation
  interpolation_method: "local_mean"
  sigma_clip: True
  no_sigma_clip_on_galaxy: False
  polynomial_on_galaxy: True
  
  # Debug mode
  debug:
  {
    no_segment_before: False
    no_segment_after: False
    no_segment: False
    expand: False
    success: False
    dilated: False
    
    user_expansion: False
    
    overlap_before: False
    overlap_after: False
  }
  
  dilate: True
  dilation_factor: 1.4
  iterations: 5
  connectivity: 2
  
  # User expansion
  user_expansion: False
  user_expansion_factor: None
  
  # Remove if not fitted
  remove_if_not_fitted: True
  
  # Remove if undetected
  remove_if_undetected: False
  
  # Remove appendages from overlapping mask
  remove_appendages: True
  
  # Remove foreground stars
  remove_foreground: True
  
  check_centroid: True
  max_centroid_offset: 10.0
  max_centroid_ellipticity: 0.3
  
  # Apertures
  apertures:
  {
    sigma_level: 4.0 # approximate isophotal extent
    max_ellipticity: 0.1
  
    # Maximal offset between the aperture center and star position (in number of pixels) (None=no limit)
    max_offset: 10.0
  }
  
  # Aperture removal
  remove_apertures: False
  aperture_removal:
  {
    # Expansion factor
    expansion_factor: 1.0
  
    # Background outer factor
    background_outer_factor: 1.2
  
    # Sigma-clipping
    no_sigma_clip_on_galaxy: True
    sigma_clip: True
  
    # Interpolation
    polynomial_on_galaxy: True
    interpolation_method: "local_mean"
  }

  second_segmentation: False
  second_sigma_level: 1.2
}

# Region
region:
{
  sigma_level: 5.0
}

# Calculation of the default FWHM
fwhm:
{
  measure: "mean"    # other options: "max", "median"
  scale_factor: 1.0
}

# -----------------------------------------------------------------
