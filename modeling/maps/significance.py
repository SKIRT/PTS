#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.significance Contains the SignificanceMaskCreator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...core.tools import sequences, types
from ...core.tools import filesystem as fs
from ...magic.core.mask import intersection
from ...magic.core.list import NamedFrameList
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class SignificanceMaskCreator(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SignificanceMaskCreator, self).__init__(*args, **kwargs)

        # The origins
        self.origins_dict = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the origin dicts
        self.load_origins()

        # 3. Load the significance maps
        self.load_significance()

        # 4. Make the significance maps
        self.make_masks()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This functino ...
        :param kwargs: 
        :return: 
        """

        # Call the setup of the base class
        super(SignificanceMaskCreator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_origins(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the origins of each map ...")

        # Get the origins
        self.origins_dict = self.get_origins_sub_names()

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        filters = []
        for sub_name in self.origins_dict:

            # Get
            origins = self.origins_dict[sub_name]

            # Check if dicts of dict
            if types.is_dictionary_of_dictionaries(origins):

                # Loop over the methods
                for method in origins:

                    # Get origins
                    origins_method = origins[method]

                    # Loop over the map names
                    for map_name in origins_method:

                        # Get the origins of the map
                        origins_map = origins_method[map_name]

                        # Extend the filters
                        sequences.extend_unique(filters, origins_map)

            # Not a dict of dicts, but a dict of strings
            else:

                # Extend the filters
                for map_name in origins:

                    # Get the origins of the map
                    origins_map = origins[map_name]

                    # Extend the filters
                    sequences.extend_unique(filters, origins_map)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    def load_significance(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the significance maps ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def map_paths_and_origins(self):

        """
        This function ...
        :return: 
        """

        # Loop over the maps subdirectories
        for sub_name in self.origins_dict:

            # Set
            maps_sub_path = fs.join(self.maps_path, sub_name)

            # Get
            origins = self.origins_dict[sub_name]

            # Check if dicts of dict
            if types.is_dictionary_of_dictionaries(origins):

                # Loop over the methods
                for method in origins:

                    # Get origins
                    origins_method = origins[method]

                    # Loop over the map names
                    for map_name in origins_method:

                        # Set the map path
                        map_path = fs.join(maps_sub_path, method, map_name + ".fits")

                        # Get the origins of the map
                        origins_map = origins_method[map_name]

                        # Yield
                        yield map_path, origins_map

            # Not a dict of dicts, but a dict of strings
            else:

                # Loop over the map names
                for map_name in origins:

                    # Set the map path
                    map_path = fs.join(maps_sub_path, map_name + ".fits")

                    # Get the origins of the map
                    origins_map = origins[map_name]

                    # Yield
                    yield map_path, origins_map

    # -----------------------------------------------------------------

    def make_masks(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the masks ...")

        # Loop over the unique
        for path, origins in self.map_paths_and_origins:

            # Get frame list
            frames = self.dataset.get_framelist_for_filters(origins)

            # Get error map list
            errors = self.dataset.get_errormaplist_for_filters(origins)

            # Convolve to same resolution
            frames.convolve_to_highest_resolution()

            # Convolve
            frames.convolve_to_highest_fwhm()
            errors.convolve_to_highest_fwhm()

            # Rebin the frames to the same pixelgrid
            frames.rebin_to_highest_pixelscale()
            errors.rebin_to_highest_pixelscale()

            # Convert the frames to the same unit
            frames.convert_to_same_unit(unit="Jy")
            errors.convert_to_same_unit(unit="Jy")

            # Get the significance maps
            significances = NamedFrameList()
            for name in frames.names:
                frame = frames[name]
                errormap = errors[name]
                significances.append(frame / errormap, name=name)

            # The masks for this map
            cutoff_masks = dict()

            # Loop over the significance level combinations for the different origins
            for sigma_levels in sequences.combinations(self.config.sigma_levels, len(frames), repeat=True):

                # Create the masks
                masks = []
                combination_names = []
                for index in range(len(significances)):
                    sigma_level = sigma_levels[index]
                    significance = significances[index]
                    string = significances.names[index] + str(sigma_level)
                    combination_names.append(string)
                    mask = significance > sigma_level
                    masks.append(mask)

                # Determine name for the mask
                mask_name = "_".join(combination_names)

                # Create intersection mask
                mask = intersection(*masks)

                # Fill holes
                mask.fill_holes()

                # Invert
                mask.invert()

                # Add the mask
                cutoff_masks[mask_name] = mask

            # Write the cutoff masks
            directory_path = fs.directory_of(path)
            masks_directory_name = fs.strip_extension(fs.name(path)) + "_masks"
            masks_directory_path = fs.create_directory_in(directory_path, masks_directory_name)

            # Loop over the masks
            for mask_name in cutoff_masks:

                # Determine the path
                mask_path = fs.join(masks_directory_path, mask_name + ".fits")

                # Write the mask
                cutoff_masks[mask_name].saveto(mask_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This funtion ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
