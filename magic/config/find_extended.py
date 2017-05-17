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

definition.add_optional("input_path", "directory_path", "path to the input directory")
definition.add_optional("output_path", "directory_path", "path to the output directory")

definition.add_flag("track_record", "track record")
definition.add_flag("plot_track_record_if_exception", "plot track record if exception", True)

definition.add_flag("find_apertures", "find apertures", True)

definition.add_optional("principal_region", "file_path", "path to a region file with a contour for the principal galaxy")

definition.add_flag("remove_apertures", "remove apertures")

definition.add_section("aperture_removal", "aperture removal")

definition.sections["aperture_removal"].add_optional("expansion_factor", "real", 1.0)

definition.add_section("fetching", "fetching")
#definition.sections["fetching"].add_flag("use_catalog_file", "use catalog file")
#definition.sections["fetching"].add_optional("catalog_path", "file_path", "catalog path")

definition.add_section("detection", "detection")

definition.sections["detection"].add_flag("use_d25", "use D25")
definition.sections["detection"].add_optional("d25_expansion_factor", "real", "D25 expansion factor", 1.2)

definition.sections["detection"].add_optional("initial_radius", "real", 20.0)
definition.sections["detection"].add_optional("detection_method", "string", "detection method", "segmentation")

definition.sections["detection"].add_flag("allow_overlap", "Do not normally allow overlap between the center segment and the background mask of the source")

definition.sections["detection"].add_flag("expand", "expand", True)

definition.sections["detection"].add_flag("always_subtract_background", "always subtract background", False)

definition.sections["detection"].add_optional("background_outer_factor", "real", "background outer factor", 1.3)

definition.sections["detection"].add_optional("background_est_method", "string", "polynomial")

definition.sections["detection"].add_flag("sigma_clip_background", "sigma clip background", True)

definition.sections["detection"].add_optional("min_pixels", "integer", "Minimum connected pixels", 5)

definition.sections["detection"].add_section("kernel", "kernel")

definition.sections["detection"].sections["kernel"].add_optional("fwhm", "real", "FWHM", 10.0)
definition.sections["detection"].sections["kernel"].add_optional("cutoff_level", "real", "cutoff_level", 4.0)

definition.sections["detection"].add_optional("sigma_level", "real", "threshold sigmas", 2.0)
definition.sections["detection"].add_optional("expansion_factor", "real", "expansion factor", 1.5)
definition.sections["detection"].add_optional("max_expansion_level", "integer", "maximum expansion level", 4)

definition.sections["detection"].add_section("debug", "debug")
definition.sections["detection"].sections["debug"].add_flag("no_segment_before", "no segment before")
definition.sections["detection"].sections["debug"].add_flag("no_segment_after", "no segment after")
definition.sections["detection"].sections["debug"].add_flag("no_segment", "no segment")
definition.sections["detection"].sections["debug"].add_flag("expand", "expand")
definition.sections["detection"].sections["debug"].add_flag("success", "success")
definition.sections["detection"].sections["debug"].add_flag("dilated", "dilated") 

definition.sections["detection"].sections["debug"].add_flag("user_expansion", "user expansion")
definition.sections["detection"].sections["debug"].add_flag("overlap_before", "overlap before")
definition.sections["detection"].sections["debug"].add_flag("overlap_after", "overlap after")


definition.sections["detection"].add_flag("dilate", "dilate", True)
definition.sections["detection"].add_optional("dilation_factor", "real", "dilation factor", 1.3)
definition.sections["detection"].add_optional("iterations", "integer", "iterations", 5)
definition.sections["detection"].add_optional("connectivity", "integer", "connectivity", 2)

definition.sections["detection"].add_flag("user_expansion", "user expansion")
definition.sections["detection"].add_optional("user_expansion_factor", "real", "user expansion factor")
  
definition.sections["detection"].add_flag("remove_appendages", "remove appendages from overlapping mask")

definition.add_section("region", "region")
definition.sections["region"].add_optional("default_radius", "real", "default radius", 20.0)

definition.add_section("apertures", "apertures")
definition.sections["apertures"].add_optional("sigma_level", "real", "approximate isophotal extent", 4.0)
definition.sections["apertures"].add_optional("max_offset", "real", "maximal offset between the aperture center and galaxy position (in number of pixels) (None=no limit)")

# Flags
definition.add_flag("weak", "weak search", False)

# -----------------------------------------------------------------
