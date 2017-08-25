#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.modeling.maps.collection import MapsCollection
from pts.core.remote.host import find_host_ids
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

default_sigma_level = 3.0

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# Create the maps collection
collection = MapsCollection.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

# Get maps
old_map_paths = collection.get_old_stellar_disk_map_paths()
young_map_paths = collection.get_young_map_paths(flatten=True)
ionizing_map_paths = collection.get_ionizing_map_paths(flatten=True)
dust_map_paths = collection.get_not_hot_dust_map_paths(flatten=True)

# Get map names
old_map_names = old_map_paths.keys()
young_map_names = young_map_paths.keys()
ionizing_map_names = ionizing_map_paths.keys()
dust_map_names = dust_map_paths.keys()

# -----------------------------------------------------------------

# Steps
definition.add_flag("steps", "save the results of intermediate steps", True)

# Remote
definition.add_optional("remote", "string", "remote host to use for creating the clip masks", choices=find_host_ids(schedulers=False))

# CONVOLUTION
definition.add_flag("convolve", "perform convolution during the creation of the clip masks", True)

# REBINNING
definition.add_optional("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)

# Selections
definition.add_optional("old", "string_list", "selected old stellar maps", choices=old_map_names)
definition.add_optional("young", "string_list", "selected young stellar maps", choices=young_map_names)
definition.add_optional("ionizing", "string_list", "selected ionizing stellar maps", choices=ionizing_map_names)
definition.add_optional("dust", "string_list", "selected dust maps", choices=dust_map_names)

# Anti-selections
definition.add_optional("not_old", "string_list", "ignore old stellar maps", choices=old_map_names)
definition.add_optional("not_young", "string_list", "ignore young stellar maps", choices=young_map_names)
definition.add_optional("not_ionizing", "string_list", "ignore ionizing stellar maps", choices=ionizing_map_names)
definition.add_optional("not_dust", "string_list", "ignore dust maps", choices=dust_map_names)

# Random selections
definition.add_optional("random_old", "positive_integer", "select x random old stellar maps")
definition.add_optional("random_young", "positive_integer", "select x random young stellar maps")
definition.add_optional("random_ionizing", "positive_integer", "select x random ionizing stellar maps")
definition.add_optional("random_dust", "positive_integer", "select x random dust maps")
definition.add_optional("random", "positive_integer", "select x maps for old stars, young stars, ionizing stars and dust")

# Flags
definition.add_flag("all_old", "select all old stellar maps")
definition.add_flag("all_young", "select all young stellar maps")
definition.add_flag("all_ionizing", "select all ionizing stellar maps")
definition.add_flag("all_dust", "select all dust maps")
definition.add_flag("all", "select all maps")

# Levels
definition.add_optional("levels", "string_real_dictionary", "significance levels for the different images")
definition.add_optional("default_level", "real", "default significance level", default_sigma_level)
definition.add_flag("all_levels", "use the default sigma level for all maps")

# CROPPING
definition.add_optional("cropping_factor", "positive_real", "multiply the cropping box with this factor", 1.3)

# Image edge softening
definition.add_optional("softening_start", "real", "relative radius for softening to start (relative to truncation ellipse)", 0.75)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Scale heights
degeyter_ratio = 8.26
#mosenkov_ratio = None
#definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=mosenkov_ratio, suggestions=[mosenkov_ratio, degeyter_ratio])
definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=degeyter_ratio, suggestions=[degeyter_ratio])
definition.add_optional("young_scaleheight_ratio", "real", "ratio of the young stellar scaleheight to the old stellar scaleheight", 0.5)
definition.add_optional("ionizing_scaleheight_ratio", "real", "ratio of the ionizing scaleheight to the old stellar scaleheight", 0.25)
definition.add_optional("dust_scaleheight_ratio", "real", "ratio of the dust scaleheight to the old stellar scaleheight", 0.5)

# -----------------------------------------------------------------

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# -----------------------------------------------------------------
