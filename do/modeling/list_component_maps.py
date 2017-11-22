#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.list_component_maps List the component maps.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.magic.tools.info import get_image_info_from_header_file
from pts.core.basics.composite import SimplePropertyComposite
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
collection = environment.static_maps_collection
selection = environment.static_maps_selection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Get configuration
config = parse_arguments("list_component_maps", definition)

# -----------------------------------------------------------------

info_kwargs = dict()
info_kwargs["path"] = False
info_kwargs["name"] = False
info_kwargs["xsize"] = True
info_kwargs["ysize"] = True
info_kwargs["psf_filter"] = False
info_kwargs["filesize"] = True
info_kwargs["filter"] = False
info_kwargs["wavelength"] = False
info_kwargs["unit"] = False

# -----------------------------------------------------------------

def get_highest_pixelscale_filter(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    highest_pixelscale = None
    highest_pixelscale_filter = None

    for fltr in filters:

        image_name = str(fltr)

        # Show origin with highest pixelscale
        #path = environment.get_frame_path_for_filter(fltr)
        path = environment.get_frame_path(image_name) # FASTER

        info = get_image_info_from_header_file(image_name, path, name=False, filter=False, wavelength=False, unit=False, fwhm=False, psf_filter=False, xsize=True, ysize=True, filesize=False, pixelscale=True)
        pixelscale = info["Pixelscale"]

        if highest_pixelscale is None or pixelscale > highest_pixelscale:
            highest_pixelscale = pixelscale
            highest_pixelscale_filter = fltr

    return highest_pixelscale_filter

# -----------------------------------------------------------------

def show_info(filepath, origins=None):

    """
    This function ...
    :param filepath:
    :param origins:
    :return:
    """

    # Get map info
    info = get_image_info_from_header_file(name, filepath, **info_kwargs)

    # Make property composite
    info = SimplePropertyComposite.from_dict(info)

    # Show the properties
    print(info.to_string(line_prefix="   ", bullet="*"))

    # Show origins
    print("    * " + fmt.bold + "origins" + fmt.reset + ": " + tostr(origins, delimiter=", "))

    # Show origin with highest pixelscale
    highest_pixelscale_filter = get_highest_pixelscale_filter(origins)

    print("    * " + fmt.bold + "pixelgrid reference" + fmt.reset + ": " + str(highest_pixelscale_filter))

    # Show origin with poorest resolution

# -----------------------------------------------------------------

print("")

# Old
print(fmt.blue + "OLD STELLAR DISK" + fmt.reset)
print("")
for name in selection.old_map_names:

    print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
    print("")

    # Get path
    path = selection.old_map_paths[name]

    # Get origins
    origins = collection.get_old_stellar_disk_origins()[name]

    # Show the info
    show_info(path, origins)
    print("")

# -----------------------------------------------------------------

# Young
print(fmt.blue + "YOUNG STELLAR DISK" + fmt.reset)
print("")
for name in selection.young_map_names:

    print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
    print("")

    # Get path
    path = selection.young_map_paths[name]

    # Get origins
    origins = collection.young_origins_flat[name]

    # Show the info
    show_info(path, origins)
    print("")

# -----------------------------------------------------------------

# Ionizing
print(fmt.blue + "IONIZING STELLAR DISK" + fmt.reset)
print("")
for name in selection.ionizing_map_names:

    print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
    print("")

    # Get path
    path = selection.ionizing_map_paths[name]

    # Get origins
    origins = collection.ionizing_origins_flat[name]

    # Show the info
    show_info(path, origins)
    print("")

# -----------------------------------------------------------------

# Dust
print(fmt.blue + "DUST DISK" + fmt.reset)
print("")
for name in selection.dust_map_names:

    print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
    print("")

    # Get path
    path = selection.dust_map_paths[name]

    # Get origins
    origins = collection.dust_origins_flat[name]

    # Show the info
    show_info(path, origins)
    print("")

# -----------------------------------------------------------------
