#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.maps_similarity

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
#from skimage.measure import structural_similarity as ssim
from skimage.measure import compare_ssim

# Import the relevant PTS classes and modules
from pts.core.tools.logging import setup_log
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.magic.core.frame import Frame
from pts.modeling.maps.collection import MapsCollection
from pts.core.tools import formatting as fmt
from pts.magic.core.list import NamedFrameList

# -----------------------------------------------------------------

maps = ["old", "young", "ionizing", "dust"]

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_optional("maps", "string_list", "maps for which to show the similarity", choices=maps, default=maps)

# Create the configuration
setter = ArgumentConfigurationSetter("maps_similarity")
config = setter.run(definition)

# -----------------------------------------------------------------

log = setup_log(level="DEBUG")

# -----------------------------------------------------------------

# Create the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# Load the maps collection
collection = MapsCollection(environment.maps_path)

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M81
m81_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "M81")

# Determine path to maps directory
maps_path = fs.join(m81_data_path, "maps")

# -----------------------------------------------------------------

old_filename = "old_stars.fits"
young_filename = "young_stars.fits"
ionizing_filename = "ionizing_stars.fits"
dust_filename = "dust.fits"

# -----------------------------------------------------------------

# The paths
map_paths = dict()
map_paths["old"] = fs.join(maps_path, old_filename)
map_paths["young"] = fs.join(maps_path, young_filename)
map_paths["ionizing"] = fs.join(maps_path, ionizing_filename)
map_paths["dust"] = fs.join(maps_path, dust_filename)

# -----------------------------------------------------------------

def mse(image_a, image_b):

    """
    This function ...
    :param image_a:
    :param image_b:
    :return:
    """

    # the 'Mean Squared Error' between the two images is the
    # sum of the squared difference between the two images;
    # NOTE: the two images must have the same dimension
    err = np.sum((image_a - image_b) ** 2)
    err /= float(image_a.shape[0] * image_a.shape[1])

    # return the MSE, the lower the error, the more "similar"
    # the two images are
    return err

# -----------------------------------------------------------------

# Loop over the different kinds of maps
for which_map in config.maps:

    # Get the map path
    path = map_paths[which_map]

    # Load the map
    the_map = Frame.from_file(path)

    # Set psf filter
    the_map.psf_filter = "Pacs red"

    #print("REFERENCE WCS:", the_map.wcs)

    # Load the paths to the created maps
    #if which_map == "dust": paths = get_dust_map_paths(modeling_path)
    #elif which_map == "old": paths = get_old_stellar_map_paths(modeling_path)
    #elif which_map == "young": paths = get_young_stellar_map_paths(modeling_path)
    #elif which_map == "ionizing": paths = get_ionizing_stellar_map_paths(modeling_path)
    #else: raise ValueError("Invalid map type: '" + which_map + "'")

    # Loop over the maps
    if which_map == "dust": maps = collection.get_dust_maps(flatten=True)
    elif which_map == "old": maps = collection.get_old_maps(flatten=True)
    elif which_map == "young": maps = collection.get_young_maps(flatten=True)
    elif which_map == "ionizing": maps = collection.get_ionizing_maps(flatten=True)
    else: raise ValueError("Invalid map type: '" + which_map + "'")

    similarities = dict()

    # No maps
    if len(maps) == 0:
        log.error("No " + which_map + " maps (yet)")
        continue

    # Loop over the maps
    for name in maps:

        # Get the map
        comparison_map = maps[name]

        #print(comparison_map.wcs.is_celestial)
        #print(comparison_map.wcs.has_celestial)

        if not comparison_map.wcs.is_celestial:
            log.warning("The " + name + " " + which_map + " dust map doesn't have a celestial WCS: skipping ...")
            continue

        # Debugging
        log.debug("Bringing the reference and comparison image to the same resolution ...")

        # Rebin to same pixel grid
        frames = NamedFrameList(reference=the_map, comparison=comparison_map)
        #frames.set_uniform_properties()

        #frames.show_coordinate_systems()

        frames.convolve_and_rebin()

        # Normalize both frames
        frames.normalize()

        # REplace nans
        frames.replace_nans(0.0)

        # Get the data
        # SAME DTYPE
        reference = frames["reference"].data.astype("float")
        comparison = frames["comparison"].data.astype("float")

        #print(reference.dtype, comparison.dtype)

        # Calculate Chi squared ?

        # Calculate Mean Squared error
        mean_sq = mse(reference, comparison)

        # Calcualte Structural Similarity Index
        struct_sim = compare_ssim(reference, comparison)

        # Add to dict
        similarities[name] = (mean_sq, struct_sim)

    # Ordered names, sorted on SSIM
    ordered_names = sorted(similarities.keys(), key=lambda name: similarities[name][1], reverse=True)

    print("")
    print(fmt.green + fmt.underlined + which_map + " maps:" + fmt.reset)
    print("")

    # Show
    with fmt.print_in_columns(4, indent="  ") as print_row:

        # Print
        for name in ordered_names: print_row(name, ":", str(similarities[name][0]), str(similarities[name][1]))

    print("")

# -----------------------------------------------------------------
