#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids
from pts.modeling.config.maps import definition
from pts.modeling.core.environment import verify_modeling_cwd
from pts.magic.core.rgba import alpha_methods

# -----------------------------------------------------------------

#relative_sigma_levels = [0.2, 0.3, 0.4, 0.5, 0.75, 0.85, 1.]
relative_sigma_levels = [0.5, 0.75, 1.]
default_relative_sigma_level = 1.0

# -----------------------------------------------------------------

scales = ["log", "sqrt"]
default_colour = "jet"
default_interval = "pts"

# -----------------------------------------------------------------

default_mask_color = "black"

# -----------------------------------------------------------------

# Set the modeling path
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

# Show?
definition.add_flag("show", "show the page", False)

# Remote
definition.add_optional("remote", "string", "remote host to use for creating the clip masks", choices=find_host_ids(schedulers=False))

# Flags
definition.add_flag("convolve", "perform convolution during the creation of the clip masks", False)

# CROPPING
definition.add_optional("cropping_factor", "positive_real", "multiply the cropping box with this factor", 1.3)

# REBINNING
definition.add_optional("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)

# Flags
definition.add_flag("add_old", "add old stellar maps", True)
definition.add_flag("add_young", "add young stellar maps", True)
definition.add_flag("add_ionizing", "add ionizing stellar maps", True)
definition.add_flag("add_dust", "add dust maps", True)

# ONLY
map_types = ["old", "young", "ionizing", "dust"]
definition.add_optional("only", "string_list", "only add these types of maps", choices=map_types)

# Sigma levels
definition.add_positional_optional("sigma_levels", "ascending_real_list", "different sigma levels for which to generate significance masks", relative_sigma_levels)
definition.add_optional("default_sigma_level", "real", "default sigma level", default_relative_sigma_level)

# Additional levels for some images
definition.add_optional("additional_levels", "filter_real_list_dictionary", "additional relative sigma levels for some images")

# DEFAULT SIGMA LEVELS FROM LEVELS FILE

# Get indices of level files
current_indices = fs.files_in_path(maps_components_path, extension="dat", returns="name", startswith="levels", convert=int, sort=int, convert_split_index=1, convert_split_pattern="_")

# There are current level files
if len(current_indices) == 0: definition.add_fixed("default_levels_from", "use default levels from levels file", None)

# Just one level file: use it
elif len(current_indices) == 1: definition.add_fixed("default_levels_from", "use default levels from levels file", current_indices[0])

# More than one
else: definition.add_optional("default_levels_from", "positive_integer", "use default levels from levels file", choices=current_indices)

# Flags
definition.add_flag("replot", "replot already existing figures", False)
definition.add_flag("replot_old", "replot already exising old stellar map plots", False)
definition.add_flag("replot_young", "replot already existing young stellar map plots", False)
definition.add_flag("replot_ionizing", "replot already existing ionizing stellar map plots", False)
definition.add_flag("replot_dust", "replot already existing dust map plots")
definition.add_flag("replot_image_masks", "replot the image masks")

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Image
definition.add_optional("image_width", "positive_integer", "width of the image")
definition.add_optional("image_height", "positive_integer", "height of the image", 300)

# -----------------------------------------------------------------

# Write data
definition.add_flag("write_data", "write the data in the form of FITS files", False)

# Clear data?
definition.add_flag("clear_data", "clear all data", False)
definition.add_flag("clear_old_data", "clear all data from old stellar maps", False)
definition.add_flag("clear_young_data", "clear all data from young stellar maps", False)
definition.add_flag("clear_ionizing_data", "clear all data from ionizing stellar maps", False)
definition.add_flag("clear_dust_data", "clear all data from dust maps", False)

# -----------------------------------------------------------------

# For masks
definition.add_optional("mask_colour", "string", "colour for the mask", default=default_mask_color)

# For masking
fuzzy_mask = False
definition.add_flag("fuzzy_mask", "use fuzzy masks", fuzzy_mask)
definition.add_optional("fuzziness", "percentage", "relative fuzziness edge width", "50", convert_default=True)
definition.add_optional("fuzzy_min_significance_offset", "positive_real", "minimum significance offset from start of fuzzy edge to maximum (peak) significance (in sigma levels)", 1.)

# Dilate masks?
definition.add_flag("dilate_masks", "dilate regular masks", True)
definition.add_flag("dilate_fuzzy_masks", "dilate alpha masks", True)

# Dilation relative radii
definition.add_optional("relative_dilation_radius_old", "positive_real", "dilation radius relative to old stellar map xsize", 1./60.)
definition.add_optional("relative_dilation_radius_young", "positive_real", "dilation radius relative to young stellar map xsize", 1./25.)
definition.add_optional("relative_dilation_radius_ionizing", "positive_real", "dilation radius relative to ionizing stellar map xsize", 1./25.)
definition.add_optional("relative_dilation_radius_dust", "positive_real", "dilation radius relative to dust map xsize", 1./50.)

# Soften regular masks
definition.add_flag("soften_masks", "soften regular masks", True)
definition.add_optional("relative_softening_radius_old", "positive_real", "softening radius relative to old stellar map xsize", 1./25.)
definition.add_optional("relative_softening_radius_young", "positive_real", "softening radius relative to young stellar map xsize", 1./15)
definition.add_optional("relative_softening_radius_ionizing", "positive_real", "softening radius relative to ionizing stellar map xsize", 1./20)
definition.add_optional("relative_softening_radius_dust", "positive_real", "softening radius relative to dust map xsize", 1./15)

# -----------------------------------------------------------------

# For PNG
default_alpha_method = "combined"
default_peak_alpha = 1.5
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, choices=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", default_peak_alpha)

# -----------------------------------------------------------------

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# -----------------------------------------------------------------

# ADVANCED
definition.add_flag("reclip_from_masks", "reclip from saved masks", False)
definition.add_optional("data_from", "string", "get the data from this directory (relative to maps html directory of absolute)")
definition.add_flag("remove_other_data", "remove data that is not needed for the selected levels", False)
definition.add_flag("resoften_masks", "resoften masks", False)

# -----------------------------------------------------------------

# Compactness
definition.add_optional("old_compactness_factor", "positive_real", "relative fraction of the truncation ellipse to take as the boundary of the old stellar maps", 1.)
definition.add_optional("young_compactness_factor", "positive_real", "relative fraction of the truncation ellipse to take as the boundary of the young stellar maps", 1.)
definition.add_optional("ionizing_compactness_factor", "positive_real", "relative fraction of the truncation ellipse to take as the boundary of the ionizing stellar maps", 1.)
definition.add_optional("dust_compactness_factor", "positive_real", "relative fraction of the truncation ellipse to take as the boundary of the dust maps", 0.75)

# -----------------------------------------------------------------

# Plotting for debugging
definition.add_flag("plot_clipping_old", "plot clipping steps for old stellar maps")
definition.add_flag("plot_clipping_young", "plot clipping steps for young stellar maps")
definition.add_flag("plot_clipping_ionizing", "plot clipping steps for ionizing stellar maps")
definition.add_flag("plot_clipping_dust", "plot clipping steps for dust maps")

# -----------------------------------------------------------------

# Ignore filters
definition.add_optional("ignore_filters", "filter_list", "ignore these filters for making the clip maps")

# -----------------------------------------------------------------
