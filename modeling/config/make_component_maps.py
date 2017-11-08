#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from pts.modeling.maps.collection import MapsCollection
from pts.core.remote.host import find_host_ids
from pts.modeling.config.maps import definition
from pts.modeling.maps.components import steps
from pts.magic.core.cutout import interpolation_methods
from pts.modeling.config.build_stars import degeyter_ratio, scalelength_scaleheight_ratios
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.component.galaxy import has_bulge2d_model, has_disk2d_model, get_bulge2d_model, get_disk2d_model
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

default_sigma_level = 3.0

# -----------------------------------------------------------------

# Get the modeling path
modeling_path = verify_modeling_cwd()

# Set paths
maps_path = fs.join(modeling_path, "maps")
maps_components_path = fs.join(maps_path, "components")

# -----------------------------------------------------------------

# Get selection indices
indices = fs.files_in_path(maps_components_path, extension="dat", returns="name", startswith="selection", convert=int, sort=int, convert_split_index=1, convert_split_pattern="_")

# Add the setting
if len(indices) == 0: raise RuntimeError("Could not find any selection file")
elif len(indices) == 1: definition.add_fixed("selection", "selection to use", indices[0])
else: definition.add_optional("selection", "positive_integer", "selection to use", choices=indices)

# -----------------------------------------------------------------

# Steps
definition.add_flag("steps", "save the results of intermediate steps", True)

# Remote
definition.add_optional("remote", "string", "remote host to use for creating the clip masks", choices=find_host_ids(schedulers=False))

# CONVOLUTION
definition.add_flag("convolve", "perform convolution during the creation of the clip masks", True)

# REBINNING
definition.add_optional("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)

# -----------------------------------------------------------------

# LEVELS

# Use previous levels
current_indices = fs.files_in_path(maps_components_path, extension="dat", returns="name", startswith="levels", convert=int, sort=int, convert_split_index=1, convert_split_pattern="_")
if len(current_indices) > 0: definition.add_optional("previous_levels", "positive_integer", "use previous levels", choices=current_indices)
else: definition.add_fixed("previous_levels", "use previous levels", None)

# -----------------------------------------------------------------

# Levels
definition.add_optional("levels", "filter_real_dictionary", "significance levels for the different images")
definition.add_optional("default_level", "real", "default significance level", default_sigma_level)
definition.add_flag("all_levels", "use the default sigma level for all maps")

# CROPPING
definition.add_optional("cropping_factor", "positive_real", "multiply the cropping box with this factor", 1.3)

# Image edge softening
definition.add_optional("softening_start", "real", "relative radius for softening to start (relative to truncation ellipse)", 0.75)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Scale heights
definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=degeyter_ratio, suggestions=scalelength_scaleheight_ratios)
definition.add_optional("young_scaleheight_ratio", "real", "ratio of the young stellar scaleheight to the old stellar scaleheight", 0.5)
definition.add_optional("ionizing_scaleheight_ratio", "real", "ratio of the ionizing scaleheight to the old stellar scaleheight", 0.25)
definition.add_optional("dust_scaleheight_ratio", "real", "ratio of the dust scaleheight to the old stellar scaleheight", 0.5)

# -----------------------------------------------------------------

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# For masking
definition.add_flag("fuzzy_mask", "use fuzzy masks", True)
definition.add_optional("fuzziness", "percentage", "relative fuzziness edge width", "50", convert_default=True)
definition.add_optional("fuzzy_min_significance_offset", "positive_real", "minimum significance offset from start of fuzzy edge to maximum (peak) significance (in sigma levels)", 1.)

# -----------------------------------------------------------------

# ADVANCED
definition.add_optional("rerun", "string", "rerun the map processing (for all maps) from this step", choices=steps)
definition.add_optional("rerun_old", "string", "rerun the map processing (for all old stellar maps) from this step", choices=steps)
definition.add_optional("rerun_young", "string", "rerun the map processing (for all young stellar maps) from this step", choices=steps)
definition.add_optional("rerun_ionizing", "string", "rerun the map processing (for all ionizing stellar maps) from this step", choices=steps)
definition.add_optional("rerun_dust", "string", "rerun the map processing (for all dust maps) from this step", choices=steps)

# ADVANCED
# to save space
definition.add_flag("remove_other", "remove maps, masks and intermediate results for maps other than those that are selected", False)
definition.add_flag("remove_other_old", "remove other old stellar maps", False)
definition.add_flag("remove_other_young", "remove other young stellar maps", False)
definition.add_flag("remove_other_ionizing", "remove other ionizing stellar maps", False)
definition.add_flag("remove_other_dust", "remove other dust maps", False)

# -----------------------------------------------------------------

# REDEPROJECT

definition.add_flag("redeproject", "redeproject all maps")
definition.add_flag("redeproject_old", "redeproject the old stellar maps")
definition.add_flag("redeproject_young", "redeproject the young stellar maps")
definition.add_flag("redeproject_ionizing", "redeproject the ionizing stellar maps")
definition.add_flag("redeproject_dust", "redeproject the dust maps")

# REDEPROJECT WITH SKIRT

definition.add_flag("redeproject_skirt", "redeproject all maps with SKIRT")
definition.add_flag("redeproject_skirt_old", "redeproject the old stellar maps with SKIRT")
definition.add_flag("redeproject_skirt_young", "redeproject the young stellar maps with SKIRT")
definition.add_flag("redeproject_skirt_ionizing", "redeproject the ionizing stellar maps with SKIRT")
definition.add_flag("redeproject_skirt_dust", "redeproject the dust maps with SKIRT")

# REPROJECT

definition.add_flag("reproject", "reproject all maps")
definition.add_flag("reproject_old", "reproject the old stellar maps")
definition.add_flag("reproject_young", "reproject the young stellar maps")
definition.add_flag("reproject_ionizing", "reproject the ionizing stellar maps")
definition.add_flag("reproject_dust", "reproject the dust maps")

# -----------------------------------------------------------------

#default_interpolation_method = "pts"
default_interpolation_method = "kernel"

# -----------------------------------------------------------------

default_core_region_factor = 0.06

# -----------------------------------------------------------------

# INTERPOLATION OF CORE OF THE MAPS
definition.add_flag("interpolate_old", "interpolate core region of old stellar maps", True)
definition.add_flag("interpolate_young", "interpolate core region of young stellar maps", False)
definition.add_flag("interpolate_ionizing", "interpolate core region of ionizing stellar maps", False)
definition.add_flag("interpolate_dust", "interpolate core region of dust maps", False)

# Interpolate negatives WITH DILATION
definition.add_flag("interpolate_old_negatives", "interpolate negatives in old stellar maps", True)
definition.add_flag("interpolate_young_negatives", "interpolate negatives in young stellar maps", True)
definition.add_flag("interpolate_ionizing_negatives", "interpolate negatives in ionizing stellar maps", True)
definition.add_flag("interpolate_dust_negatives", "interpolate negatives in dust maps", True)

# Central ellipse factor
definition.add_optional("negatives_central_ellipse_factor", "real", "factor for the central ellipse for considering negatives", 0.4)

# Dilation radius
default_negatives_dilation_radius = 10
definition.add_optional("old_negatives_dilation_radius", "real", "old negatives dilation radius", default_negatives_dilation_radius)
definition.add_optional("young_negatives_dilation_radius", "real", "young negatives dilation radius", default_negatives_dilation_radius)
definition.add_optional("ionizing_negatives_dilation_radius", "real", "ionizing negatives dilation radius", default_negatives_dilation_radius)
definition.add_optional("dust_negatives_dilation_radius", "real", "dust negatives dilation radius", default_negatives_dilation_radius)

# Interpolation core
definition.add_optional("old_core_region_factor", "real", "interpolation core boundary for the old stellar maps, relative to the truncation ellipse", default=default_core_region_factor)
definition.add_optional("young_core_region_factor", "real", "interpolation core boundary for the young stellar maps, relative to the truncation ellipse", default=default_core_region_factor)
definition.add_optional("ionizing_core_region_factor", "real", "interpolation core boundary for the ionizing stellar maps, relative to the truncation ellipse", default=default_core_region_factor)
definition.add_optional("dust_core_region_factor", "real", "interpolation core boundary for the dust maps, relative to the truncation ellipse", default=default_core_region_factor)

# Interpolation settings
definition.add_optional("source_outer_factor", "real", "outer factor", 1.4)
definition.add_optional("interpolation_method", "string", "interpolation method", default_interpolation_method, choices=interpolation_methods)
definition.add_flag("sigma_clip", "apply sigma clipping before interpolation", True)

# -18 deg # suggestion for M81 for offset old
if has_disk2d_model(modeling_path) and has_bulge2d_model(modeling_path):

    # Get the difference in position angle
    bulge = get_bulge2d_model(modeling_path)
    disk = get_disk2d_model(modeling_path)
    default_old_offset = bulge.position_angle - disk.position_angle

# No models: no offset
else: default_old_offset = Angle(0.0, "deg")

# ALSO FOR INTERPOLATION
definition.add_optional("interpolation_angle_offset_old", "angle", "offset of angle of ellipse for interpolation w.r.t. angle of truncation ellipse", default_old_offset)
definition.add_optional("interpolation_angle_offset_young", "angle", "offset of angle of ellipse for interpolation w.r.t. angle of truncation ellipse", "0 deg", convert_default=True)
definition.add_optional("interpolation_angle_offset_ionizing", "angle", "offset of angle of ellipse for interpolation w.r.t. angle of truncation ellipse", "0 deg", convert_default=True)
definition.add_optional("interpolation_angle_offset_dust", "angle", "offset of angle of ellipse for interpolation w.r.t. angle of truncation ellipse", "0 deg", convert_default=True)

# MORE FOR INTERPOLATION
definition.add_optional("interpolation_softening_start", "real", "relative radius for softening to start (relative to interpolation ellipse)", 0.65)
definition.add_optional("interpolation_softening_end", "real", "relative radius for softening to end (relative to interpolation ellipse", 1.2)

# INTERPOLATION SMOOTHING
default_smoothing_factor = 5.
definition.add_optional("old_interpolation_smoothing_factor", "real", "smoothing factor for interpolation of old stellar maps", default_smoothing_factor)
definition.add_optional("young_interpolation_smoothing_factor", "real", "smoothing factor for interpolation of young stellar maps", default_smoothing_factor)
definition.add_optional("ionizing_interpolation_smoothing_factor", "real", "smoothing factor for interpolation of ionizing stellar maps", default_smoothing_factor)
definition.add_optional("dust_interpolation_smoothing_factor", "real", "smoothing factor for interpolation of dust maps", default_smoothing_factor)

# INTERPOLATE IN CUTOUT (FOR SPEED?)
definition.add_flag("interpolate_in_cutout", "interpolate in cutouts", False)

# -----------------------------------------------------------------

# Clear results
definition.add_flag("clear_results", "clear previous results")
definition.add_flag("clear_results_old", "clear previous results of old stellar maps")
definition.add_flag("clear_results_young", "clear previous results of young stellar maps")
definition.add_flag("clear_results_ionizing", "clear previous results of ionizing stellar maps")
definition.add_flag("clear_results_dust", "clear previous results of dust maps")

# CLEAR
definition.add_flag("clear_all", "clear all previous results and steps")
definition.add_flag("clear_old", "clear all previous results and steps of old stellar maps")
definition.add_flag("clear_young", "clear all previous results and steps of young stellar maps")
definition.add_flag("clear_ionizing", "clear all previous results and steps of ionizing maps")
definition.add_flag("clear_dust", "clear all previous results and steps of dust maps")

# -----------------------------------------------------------------

# FOR DEPROJECTION
definition.add_optional("downsample_factor", "positive_real", "downsample factor for rendering the deprojected maps", 2.)

# FOR PROJECTION
definition.add_optional("scale_heights", "positive_real", "scale heights", 10.)

# -----------------------------------------------------------------

# Plot
formats = ["pdf", "png"]
default_format = "pdf"
definition.add_flag("plot", "make plots", True)
definition.add_optional("plotting_format", "string", "plotting format", default=default_format, choices=formats)

# -----------------------------------------------------------------
