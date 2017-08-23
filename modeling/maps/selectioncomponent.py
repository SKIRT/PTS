#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.selectioncomponent Contains the MapsSelectionComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from .component import MapsComponent
from ...core.basics.log import log
from ...magic.core.mask import intersection
from ...core.tools import sequences
from ...core.tools.stringify import tostr
from ...core.tools.utils import lazyproperty
from ...core.basics.containers import hashdict
from ...core.tools import types
from ...core.filter.filter import parse_filter
from ...core.basics.configuration import prompt_string_list

# -----------------------------------------------------------------

class MapsSelectionComponent(MapsComponent):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MapsSelectionComponent, self).__init__(*args, **kwargs)

        # The selections
        self.old_selection = None
        self.young_selection = None
        self.ionizing_selection = None
        self.dust_selection = None

        # The significance levels
        self.levels = None

        # The maps
        self.old_maps = dict()
        self.young_maps = dict()
        self.ionizing_maps = dict()
        self.dust_maps = dict()

        # Cached clip masks
        self.clip_masks = dict()

    # -----------------------------------------------------------------

    @property
    def has_old_selection(self):

        """
        This function ...
        :return:
        """

        return self.old_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_young_selection(self):

        """
        This function ...
        :return:
        """

        return self.young_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_ionizing_selection(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_dust_selection(self):

        """
        Thisf unctino ...
        :return:
        """

        return self.dust_selection is not None

    # -----------------------------------------------------------------

    def prompt_old(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the old stellar maps ...")

        # Ask for the old stellar maps
        self.old_selection = prompt_string_list("old_maps", "selection of old stellar disk maps", choices=self.old_map_names, all_default=True)

    # -----------------------------------------------------------------

    def prompt_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the young stellar maps ...")

        # Ask for the young stellar maps
        self.young_selection = prompt_string_list("young_maps", "selection of young stellar disk maps", choices=self.young_map_names, all_default=True)

    # -----------------------------------------------------------------

    def prompt_ionizing(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the ionizing stellar maps ...")

        # Ask for the ionizing stellar map
        self.ionizing_selection = prompt_string_list("ionizing_maps", "selection of ionizing stellar disk maps", choices=self.ionizing_map_names, all_default=True)

    # -----------------------------------------------------------------

    def prompt_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the dust maps ...")

        # Ask for the dust map to use
        self.dust_selection = prompt_string_list("dust_maps", "selection of dust disk maps", choices=self.dust_map_names, all_default=True)

    # -----------------------------------------------------------------

    @property
    def has_levels(self):

        """
        This function ...
        :return:
        """

        return self.levels is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_old_stellar_disk_map_paths()

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_origins(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_old_stellar_disk_origins()

    # -----------------------------------------------------------------

    @property
    def old_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_young_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_origins(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_young_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def young_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_paths(self):

        """
        Thisn function ...
        :return:
        """

        # Get map paths
        return self.maps_collection.get_ionizing_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_origins(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_ionizing_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_paths(self):

        """
        Thisn function ...
        :return:
        """

        # Get paths
        # map_paths = self.maps_collection.get_dust_map_paths(flatten=True)
        return self.maps_collection.get_not_hot_dust_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_origins(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_not_hot_dust_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def dust_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_paths.keys()

    # -----------------------------------------------------------------

    def correct_map(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        # Replace infinities
        the_map.replace_infs(0.0)

        # Replace NaNs
        the_map.replace_nans(0.0)

        # Replace negative values
        the_map.replace_negatives(0.0)

    # -----------------------------------------------------------------

    def get_truncation_mask_for_map(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        # Create the truncation ellipse in pixel coordinates
        ellipse = self.truncation_ellipse.to_pixel(the_map.wcs)

        # Create mask
        mask = ellipse.to_mask(the_map.xsize, the_map.ysize).inverse()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def truncate_map(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        # Get the truncation mask
        mask = self.get_truncation_mask_for_map(the_map)

        # Mask
        the_map[mask] = 0.0

    # -----------------------------------------------------------------

    def crop_map(self, the_map, factor):

        """
        This function ...
        :param the_map:
        :param factor:
        :return:
        """

        the_map.crop_to(self.truncation_box, factor=factor)

    # -----------------------------------------------------------------

    def clip_map(self, the_map, origins, convolve=True, remote=None, npixels=1, connectivity=8, rebin_remote_threshold=None):

        """
        This function ...
        :param the_map:
        :param origins:
        :param convolve:
        :param remote:
        :param npixels:
        :param connectivity:
        :param rebin_remote_threshold:
        :return:
        """

        # Create the clip mask
        mask = self.get_clip_mask(origins, wcs=the_map.wcs, convolve=convolve, remote=remote, npixels=npixels, connectivity=connectivity, rebin_remote_threshold=rebin_remote_threshold)

        # Clip
        the_map[mask] = 0.0

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def make_clipped_maps(self, the_map, origins, levels_dict, convolve=True, remote=None, rebin_remote_threshold=None,
                          npixels=1, connectivity=8, present=None):

        """
        This function ...
        :param the_map:
        :param origins:
        :param levels_dict:
        :param convolve:
        :param remote:
        :param rebin_remote_threshold:
        :param npixels:
        :param connectivity:
        :return:
        """

        # Make the masks
        masks = self.make_clip_masks(origins, levels_dict, wcs=the_map.wcs, convolve=convolve, remote=remote,
                                     rebin_remote_threshold=rebin_remote_threshold, npixels=npixels, connectivity=connectivity, present=present)

        # The maps
        maps = dict()

        # Make the maps
        for levels in masks:

            # Debugging
            log.debug("Clipping the map for sigma levels [" + tostr(levels) + "] ...")

            # Get the mask
            mask = masks[levels]

            # Make clipped copy
            clipped_map = the_map.copy()
            clipped_map[mask] = 0.0

            # Add to the dictionary
            maps[levels] = clipped_map

        # Return the dictionary
        return maps

    # -----------------------------------------------------------------

    def get_clip_mask(self, origins, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8, rebin_remote_threshold=None):

        """
        This function ...
        :param origins:
        :param wcs:
        :param convolve:
        :param remote:
        :param npixels:
        :param connectivity:
        :param rebin_remote_threshold:
        :return:
        """

        # Check if cached
        if wcs is None and tuple(origins) in self.clip_masks: return self.clip_masks[tuple(origins)]
        else:

            # Make the mask
            mask = self.make_clip_mask(origins, self.levels, wcs=wcs, convolve=convolve, remote=remote, npixels=npixels, connectivity=connectivity, rebin_remote_threshold=rebin_remote_threshold)

            # Cache the mask
            if wcs is None: self.clip_masks[tuple(origins)] = mask

            # Return the mask
            return mask

    # -----------------------------------------------------------------

    def make_clip_mask(self, origins, levels, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8, rebin_remote_threshold=None):

        """
        This function ...
        :param origins:
        :param levels:
        :param wcs:
        :param convolve:
        :param remote:
        :param npixels:
        :param connectivity:
        :param rebin_remote_threshold:
        :return:
        """

        # Debugging
        log.debug("Making a clip mask for the filters: " + tostr(origins) + " ...")

        # Get frame list
        frames = self.dataset.get_framelist_for_filters(origins)

        # Get error map list
        errors = self.dataset.get_errormaplist_for_filters(origins)

        # Convolve
        if convolve:

            frames.convolve_to_highest_fwhm(remote=remote)
            errors.convolve_to_highest_fwhm(remote=remote)

        # WCS is specified: rebin to this WCS
        if wcs is not None:

            frames.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold)
            errors.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold)

        # Otherwise, rebin to the highest pixelscale WCS
        else:

            # Rebin the frames to the same pixelgrid
            frames.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold)
            errors.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold)

        # NOT REALLY NECESSARY!!
        # Convert the frames to the same unit
        # frames.convert_to_same_unit(unit="Jy")
        # errors.convert_to_same_unit(unit="Jy")
        # frames.convert_to_same_unit()
        # errors.convert_to_same_unit()

        # # Get the significance maps
        # significances = NamedFrameList()
        # for name in frames.names:
        #     frame = frames[name]
        #     errormap = errors[name]
        #     significances.append(frame / errormap, name=name)
        #
        # # Create the masks
        # masks = []
        # combination_names = []
        #
        # for index in range(len(significances)):
        #     sigma_level = sigma_levels[index]
        #     significance = significances[index]
        #     string = significances.names[index] + str(sigma_level)
        #     combination_names.append(string)
        #     mask = significance > sigma_level
        #     masks.append(mask)

        masks = []
        for name in frames.names:
            frame = frames[name]
            errormap = errors[name]
            level = levels[frame.filter]
            mask = frame > level * errormap
            masks.append(mask)

        # Determine name for the mask
        # mask_name = "_".join(combination_names)

        # Create intersection mask
        mask = intersection(*masks)

        # Only keep central
        #mask = mask.central()
        mask = mask.largest(npixels=npixels, connectivity=connectivity)

        # Fill holes
        mask.fill_holes()

        # Invert
        mask.invert()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def make_clip_masks(self, origins, levels_dict, wcs=None, convolve=True, remote=None, rebin_remote_threshold=None,
                        npixels=1, connectivity=8, present=None):

        """
        Thisn function ...
        :param origins:
        :param levels_dict:
        :param wcs:
        :param convolve:
        :param remote:
        :parma rebin_remote_threshold:
        :param npixels:
        :param connectivity:
        :param present:
        :return:
        """

        # Debugging
        log.debug("Making clip masks for the filters: " + tostr(origins) + " ...")

        # Get frame list, IN ORDER OF ORIGINS
        frames = self.dataset.get_framelist_for_filters(origins)

        # Get error map list, IN ORDER OR ORIGINS
        errors = self.dataset.get_errormaplist_for_filters(origins)

        # Convolve
        if convolve:

            frames.convolve_to_highest_fwhm(remote=remote)
            errors.convolve_to_highest_fwhm(remote=remote)

        # WCS is specified: rebin to this WCS
        if wcs is not None:

            frames.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold)
            errors.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold)

        # Otherwise, rebin to the highest pixelscale WCS
        else:

            # Rebin the frames to the same pixelgrid
            frames.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold)
            errors.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold)

        # Get the number of frames
        names = frames.names
        nframes = len(names)

        # Get list of level lists, sorted on the order of the frames

        #frame_levels = [levels[name] for name in names]

        frame_levels = []

        # LEVELS ARE GIVEN AS A SEQUENCE: SAME FOR EACH FRAME
        if types.is_sequence(levels_dict):
            for _ in range(nframes): frame_levels.append(levels_dict)

        # LEVELS ARE GIVEN AS A DICTIONARY: DEFINED FOR EACH IMAGE NAME
        elif types.is_dictionary(levels_dict):
            #levels_dict = {parse_filter(fltr): levels for fltr, levels in levels_dict.items()}
            #for fltr in frames.filters: frame_levels.append(levels_dict[fltr])
            for name in names: frame_levels.append(levels_dict[name]) # name = origin

        # Invalid
        else: raise ValueError("Levels must be specified as list (same for each image) or dictionary keyed on filter")

        # Container to keep the set of masks
        #masks_levels = []
        masks_levels = dict()

        # Loop over each level combination
        # Loop over the significance level combinations for the different origins
        for sigma_levels in sequences.lists_combinations(*frame_levels):

            # Create dictionary that says which sigma level was used for which frame
            levels_dict = hashdict({name: level for name, level in zip(names, sigma_levels)})

            # Check
            if present is not None and levels_dict in present: continue

            # Debugging
            log.debug("Making clip mask for sigma levels [" + tostr(levels_dict) + "] ...")

            masks = []

            #for name in frames.names:
            for index in range(nframes):

                name = names[index]
                frame = frames[name]
                errormap = errors[name]
                level = sigma_levels[index]
                mask = frame > level * errormap
                masks.append(mask)

            # Create intersection mask
            mask = intersection(*masks)

            # Only keep central
            #mask = mask.central()
            mask = mask.largest(npixels=npixels, connectivity=connectivity)

            # Fill holes
            mask.fill_holes()

            # Invert
            mask.invert()

            # Add item to the list
            #item = (levels_dict, mask)
            #masks_levels.append(item)
            masks_levels[levels_dict] = mask

        # Return the masks
        return masks_levels

    # -----------------------------------------------------------------

    def soften_map(self, the_map, ellipse, softening_range):

        """
        This function ...
        :param the_map:
        :param ellipse:
        :param softening_range:
        :return:
        """

        # Get ellipse
        ellipse = ellipse.to_pixel(the_map.wcs)

        # Soften edges
        the_map.soften_edges(ellipse, softening_range)

# -----------------------------------------------------------------