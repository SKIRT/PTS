#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.maps.collection import MapsCollection
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

sigma_levels = [0.5, 0.75, 1., 1.5, 2., 2.5, 3., 3.5, 4.]
default_sigma_level = 1.0

# -----------------------------------------------------------------

scales = ["log", "sqrt"]
default_colour = "jet"
default_interval = "pts"

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

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")
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

# Sigma levels
definition.add_positional_optional("sigma_levels", "real_list", "different sigma levels for which to generate significance masks", sigma_levels)
definition.add_optional("default_sigma_level", "real", "default sigma level", default_sigma_level)

# Flags
definition.add_flag("replot", "replot already existing figures", True)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Image
definition.add_optional("image_width", "positive_integer", "width of the image")
definition.add_optional("image_height", "positive_integer", "height of the image", 300)

# For PNG
definition.add_optional("scale", "string", "scaling", "log", scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("colours", "string", "colour or colour scale", "red")
definition.add_flag("alpha", "use alpha", True)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.)

# -----------------------------------------------------------------
