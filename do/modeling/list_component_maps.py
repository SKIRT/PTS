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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.magic.tools.info import get_image_info_from_header_file
from pts.core.basics.composite import SimplePropertyComposite
from pts.core.tools.stringify import tostr
from pts.core.tools.serialization import load_dict
from pts.core.tools import filesystem as fs

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

def get_highest_fwhm_filter(filters):

    """
    Thisf unction ...
    :param filters:
    :return:
    """

    highest_fwhm = None
    highest_fwhm_filter = None

    for fltr in filters:

        image_name = str(fltr)

        path = environment.get_frame_path(image_name)

        info = get_image_info_from_header_file(image_name, path, name=False, filter=False, wavelength=False, unit=False, fwhm=True, psf_filter=True, xsize=False, ysize=False, filesize=False, pixelscale=False)
        fwhm = info["FWHM"]

        if highest_fwhm is None or fwhm > highest_fwhm:
            highest_fwhm = fwhm
            highest_fwhm_filter = fltr

    return highest_fwhm_filter

# -----------------------------------------------------------------

def show_info(filepath, origins=None, steps_path=None):

    """
    This function ...
    :param filepath:
    :param origins:
    :param steps_path:
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
    highest_fwhm_filter = get_highest_fwhm_filter(origins)
    print("    * " + fmt.bold + "PSF reference" + fmt.reset + ": " + str(highest_fwhm_filter))

    # Load clip info
    clip_info_path = fs.join(steps_path, "clipped.dat")
    clip_info = load_dict(clip_info_path)

    # Show clip image references
    clip_origins = clip_info["origins"]
    print("    * " + fmt.bold + "clipping image references" + fmt.reset + ": " + tostr(clip_origins, delimiter=", "))

    # Show levels
    clipping_path = fs.join(steps_path, "clipping")
    levels = OrderedDict()
    for origin in clip_origins:
        mask_filename = "mask__" + tostr(origin).replace(" ", "_")
        mask_filename = fs.find_file_in_path(clipping_path, extension="fits", startswith=mask_filename, returns="name")
        level = float(mask_filename.split("__")[2])
        levels[origin] = level
    print("    * " + fmt.bold + "signal-to-noise levels" + fmt.reset + ":")
    for origin in levels: print("       " + str(origin) + ": " + str(levels[origin]))

    # Show
    softened = clip_info["soften"]
    print("    * " + fmt.bold + "softened clip mask" + fmt.reset + ": " + str(softened))

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

    # Get steps path
    steps_path = selection.get_old_steps_path_for_map(name)

    # Show the info
    show_info(path, origins, steps_path)
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

    # Get steps path
    steps_path = selection.get_young_steps_path_for_map(name)

    # Show the info
    show_info(path, origins, steps_path)
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

    # Get steps path
    steps_path = selection.get_ionizing_steps_path_for_map(name)

    # Show the info
    show_info(path, origins, steps_path)
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

    # Get steps path
    steps_path = selection.get_dust_steps_path_for_map(name)

    # Show the info
    show_info(path, origins, steps_path)
    print("")

# -----------------------------------------------------------------
