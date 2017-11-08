#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.maps.collection import MapsCollection
from pts.modeling.config.maps import definition
from pts.modeling.core.environment import verify_modeling_cwd
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Get the modeling path
modeling_path = verify_modeling_cwd()

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

# Get number of maps
nold_maps = len(old_map_names)
nyoung_maps = len(young_map_names)
nionizing_maps = len(ionizing_map_names)
ndust_maps = len(dust_map_names)

# -----------------------------------------------------------------

maps_path = fs.join(modeling_path, "maps")
maps_components_path = fs.join(maps_path, "components")

# -----------------------------------------------------------------

# SELECTION

# AUTO-SELECT??
definition.add_flag("auto", "make selections automatically based on the preferred modeling guidelines", False)

# Auto selection options
default_max_nnegatives = 0.1
definition.add_optional("hot_dust_max_nnegatives", "positive_real", "maximum relative number of negatives in a central ellipse for selection of hot dust map", default_max_nnegatives)
definition.add_optional("young_max_nnegatives", "positive_real", "maximum relative number of negatives in a central ellipse for selection of young stellar map", default_max_nnegatives)

# Selections
definition.add_optional("old", "string_list", "selected old stellar maps", choices=old_map_names)
definition.add_optional("young", "string_list", "selected young stellar maps", choices=young_map_names)
definition.add_optional("ionizing", "string_list", "selected ionizing stellar maps", choices=ionizing_map_names)
definition.add_optional("dust", "string_list", "selected dust maps", choices=dust_map_names)

# Selections with indices
definition.add_optional("old_indices", "integer_list", "selected old stellar maps", choices=range(nold_maps))
definition.add_optional("young_indices", "integer_list", "selected young stellar maps", choices=range(nyoung_maps))
definition.add_optional("ionizing_indices", "integer_list", "selected ionizing stellar maps", choices=range(nionizing_maps))
definition.add_optional("dust_indices", "integer_list", "selected dust maps", choices=range(ndust_maps))

# Anti-selections
definition.add_optional("not_old", "string_list", "ignore old stellar maps", choices=old_map_names)
definition.add_optional("not_young", "string_list", "ignore young stellar maps", choices=young_map_names)
definition.add_optional("not_ionizing", "string_list", "ignore ionizing stellar maps", choices=ionizing_map_names)
definition.add_optional("not_dust", "string_list", "ignore dust maps", choices=dust_map_names)

# Anti-selections with indices
definition.add_optional("not_old_indices", "integer_list", "ignore old stellar maps", choices=range(nold_maps))
definition.add_optional("not_young_indices", "integer_list", "ignore young stellar maps", choices=range(nyoung_maps))
definition.add_optional("not_ionizing_indices", "integer_list", "ignore ionizing stellar maps", choices=range(nionizing_maps))
definition.add_optional("not_dust_indices", "integer_list", "ignore dust maps", choices=range(ndust_maps))

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

# -----------------------------------------------------------------
