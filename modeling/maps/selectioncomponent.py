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
import numpy as np
from abc import ABCMeta

# Import the relevant PTS classes and modules
from .component import MapsComponent
from ...core.basics.log import log
from ...magic.core.mask import intersection
from ...magic.core.alpha import product
from ...core.tools import sequences
from ...core.tools.stringify import tostr
from ...core.tools.utils import lazyproperty
from ...core.basics.containers import hashdict
from ...core.tools import types
from ...core.tools import filesystem as fs
from ...core.basics.configuration import prompt_string_list
from ...magic.core.mask import Mask
from ...magic.core.alpha import AlphaMask
from ...core.basics.range import RealRange

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

    @lazyproperty
    def old_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_old_stellar_disk_methods()

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

    @lazyproperty
    def young_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_young_methods(flatten=True)

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

    @lazyproperty
    def ionizing_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_ionizing_methods(flatten=True)

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

    @lazyproperty
    def dust_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_not_hot_dust_methods(flatten=True)

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

        # Debugging
        log.debug("Correcting map ...")

        # Infinities
        self.replace_infinities(the_map)

        # NaNs
        self.replace_nans(the_map)

        # Negatives
        self.replace_negatives(the_map)

    # -----------------------------------------------------------------

    def replace_infinities(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        # Debugging
        log.debug("Replacing infinities in the map ...")

        # Replace infinities
        infs = the_map.replace_infs(0.0)

        # Show number of infinities
        ninfs = np.sum(infs)
        log.debug("Map contained " + str(ninfs) + " inifinities")

    # -----------------------------------------------------------------

    def replace_nans(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        # Debugging
        log.debug("Replacing NaNs in the map ...")

        # Replace NaNs
        nans = the_map.replace_nans(0.0)

        # Show number of NaNs
        nnans = np.sum(nans)
        log.debug("Map contained " + str(nnans) + " NaN values")

    # -----------------------------------------------------------------

    def replace_negatives(self, the_map):

        """
        Thisj function ...
        :param the_map:
        :return:
        """

        # Debugging
        log.debug("Replacing negatives in the map ...")

        # Replace negative values
        negatives = the_map.replace_negatives(0.0)

        # Show number of negatives
        nnegatives = np.sum(negatives)
        log.debug("Map contained " + str(nnegatives) + " negative values")

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

    def clip_map(self, the_map, origins, convolve=True, remote=None, npixels=1, connectivity=8,
                 rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None):

        """
        This function ...
        :param the_map:
        :param origins:
        :param convolve:
        :param remote:
        :param npixels:
        :param connectivity:
        :param rebin_remote_threshold:
        :param fuzzy
        :param fuzziness:
        :param fuzziness_offset:
        :param output_path:
        :return:
        """

        # Create the clip mask
        mask = self.get_clip_mask(origins, wcs=the_map.wcs, convolve=convolve, remote=remote, npixels=npixels,
                                  connectivity=connectivity, rebin_remote_threshold=rebin_remote_threshold,
                                  fuzzy=fuzzy, fuzziness=fuzziness, fuzziness_offset=fuzziness_offset, output_path=output_path)

        # Clip
        #the_map[mask] = 0.0

        # Apply the mask
        if fuzzy: the_map.apply_alpha_mask(mask)
        else: the_map.apply_mask(mask)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def make_clipped_maps(self, name, the_map, origins, levels_dict, convolve=True, remote=None, rebin_remote_threshold=None,
                          npixels=1, connectivity=8, present=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1.,
                          return_masks=False):

        """
        This function ...
        :param name:
        :param the_map:
        :param origins:
        :param levels_dict:
        :param convolve:
        :param remote:
        :param rebin_remote_threshold:
        :param npixels:
        :param connectivity:
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
        :param return_masks:
        :return:
        """

        # Make the masks
        masks = self.make_clip_masks(name, origins, levels_dict, wcs=the_map.wcs, convolve=convolve, remote=remote,
                                     rebin_remote_threshold=rebin_remote_threshold, npixels=npixels,
                                     connectivity=connectivity, present=present, fuzzy=fuzzy, fuzziness=fuzziness,
                                     fuzziness_offset=fuzziness_offset)

        # The maps
        maps = dict()

        # Make the maps
        for levels in masks:

            # Debugging
            log.debug("Clipping the map for sigma levels [" + tostr(levels) + "] ...")

            # Get the mask
            mask = masks[levels]

            # Make clipped copy of the map
            clipped_map = the_map.copy()
            if fuzzy: clipped_map.apply_alpha_mask(mask)
            else: clipped_map.apply_mask(mask)

            # Add to the dictionary
            maps[levels] = clipped_map

        # Return the dictionary
        if return_masks: return maps, masks
        else: return maps

    # -----------------------------------------------------------------

    def get_clip_mask(self, origins, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8,
                      rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None):

        """
        This function ...
        :param origins:
        :param wcs:
        :param convolve:
        :param remote:
        :param npixels:
        :param connectivity:
        :param rebin_remote_threshold:
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
        :param output_path:
        :return:
        """

        # Check if cached
        if wcs is None and tuple(origins) in self.clip_masks: return self.clip_masks[tuple(origins)]
        else:

            # Make the mask
            mask = self.make_clip_mask(origins, self.levels, wcs=wcs, convolve=convolve, remote=remote,
                                       npixels=npixels, connectivity=connectivity,
                                       rebin_remote_threshold=rebin_remote_threshold, fuzzy=fuzzy, fuzziness=fuzziness,
                                       fuzziness_offset=fuzziness_offset, output_path=output_path)

            # Cache the mask
            if wcs is None: self.clip_masks[tuple(origins)] = mask

            # Return the mask
            return mask

    # -----------------------------------------------------------------

    def make_clip_mask(self, origins, levels, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8,
                       rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None):

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
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
        :param output_path:
        :return:
        """

        # Debugging
        log.debug("Making a clip mask for the filters: " + tostr(origins) + " ...")

        # Get frame list
        frames = self.dataset.get_framelist_for_filters(origins)

        # Get error map list
        errors = self.dataset.get_errormaplist_for_filters(origins)

        # Construct levels list, in the order of the frames
        levels_list = [levels[frames[name].filter] for name in frames.names]

        # Check whether the resulting mask is already saved
        if output_path is not None and has_resulting_mask(output_path, levels_list, fuzzy=fuzzy):

            # Return the mask
            return load_resulting_mask(output_path, levels_list, fuzzy=fuzzy)

        # Check whether the raw combination mask is already saved
        if output_path is not None and has_combination_mask(output_path, levels_list, fuzzy=fuzzy):

            # Load
            mask = load_combination_mask(output_path, levels_list, fuzzy=fuzzy)

        # Combination mask still has to be created
        else:

            # Set names to ignore for convolution and rebinning because we already have the mask
            ignore = []
            if output_path:
                for name in frames.names:
                    frame = frames[name]
                    level = levels[frame.filter]
                    if has_single_mask(output_path, name, level, fuzzy=fuzzy):
                        log.success("The clip mask for the '" + name + "' image with a sigma level of " + str(level) + " has already been created")
                        ignore.append(name)

            # Convolve
            if convolve:

                frames.convolve_to_highest_fwhm(remote=remote, ignore=ignore)
                errors.convolve_to_highest_fwhm(remote=remote, ignore=ignore)

            # WCS is specified: rebin to this WCS
            if wcs is not None:

                frames.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)
                errors.rebin_to_wcs(wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)

            # Otherwise, rebin to the highest pixelscale WCS
            else:

                # Rebin the frames to the same pixelgrid
                frames.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)
                errors.rebin_to_highest_pixelscale(remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)

            masks = []
            for name in frames.names:

                # Get frame, errormap and level
                frame = frames[name]
                errormap = errors[name]
                level = levels[frame.filter]

                # Already calculated mask
                if has_single_mask(output_path, name, level, fuzzy=fuzzy):
                    # Load the mask and continue
                    mask = load_single_mask(output_path, name, level, fuzzy=fuzzy)
                    masks.append(mask)
                    continue

                # Create significance map
                significance = frame / errormap

                # Create the mask
                if fuzzy: mask = self.create_fuzzy_mask_for_level(significance, level, fuzziness=fuzziness, offset=fuzziness_offset)
                else: mask = self.create_mask_for_level(significance, level)

                # Save the mask
                if output_path is not None: mask.saveto(get_single_mask_path(output_path, name, level, fuzzy=fuzzy))

                # Add the mask
                masks.append(mask)

            # Create intersection mask
            if fuzzy: mask = product(*masks)
            else: mask = intersection(*masks)

            # Save the mask
            if output_path is not None:
                mask_base_path = fs.join(output_path, "mask__combination__" + "_".join(repr(level) for level in levels_list))
                if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
                else: mask_path = mask_base_path + ".fits"
                mask.saveto(mask_path)

        # Only keep largest patch
        mask = mask.largest(npixels=npixels, connectivity=connectivity)

        # Fill holes
        mask.fill_holes()

        # Set the WCS
        if wcs is not None: mask.wcs = wcs

        # Invert FOR NORMAL MASKS: WE HAVE TO SET PIXELS TO ZERO THAT ARE NOT ON THE MASK
        if not fuzzy: mask.invert()

        # Save the mask
        if output_path is not None:
            mask_base_path = fs.join(output_path, "mask__result__" + "_".join(repr(level) for level in levels_list))
            if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
            else: mask_path = mask_base_path + ".fits"
            mask.saveto(mask_path)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def make_clip_masks(self, name, origins, levels_dict, wcs=None, convolve=True, remote=None, rebin_remote_threshold=None,
                        npixels=1, connectivity=8, present=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1.):

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
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
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
            if present is not None and levels_dict in present:
                log.success("The clipped '" + name + "' map for sigma levels [" + tostr(levels_dict) + "] is already present")
                continue

            # Debugging
            log.debug("Making clip mask for sigma levels [" + tostr(levels_dict) + "] ...")

            masks = []

            #for name in frames.names:
            for index in range(nframes):

                name = names[index]
                frame = frames[name]
                errormap = errors[name]
                level = sigma_levels[index]

                # Create significance map
                significance = frame / errormap

                # Create the mask
                if fuzzy: mask = self.create_fuzzy_mask_for_level(significance, level, fuzziness=fuzziness, offset=fuzziness_offset)
                else: mask = self.create_mask_for_level(significance, level)

                # Add the mask
                masks.append(mask)

            # Create intersection mask
            if fuzzy:  mask = product(*masks)
            else: mask = intersection(*masks)

            # Only keep largest patch
            mask = mask.largest(npixels=npixels, connectivity=connectivity)

            # Fill holes
            mask.fill_holes()

            # Invert FOR NORMAL MASKS: WE HAVE TO SET PIXELS TO ZERO THAT ARE NOT ON THE MASK
            if not fuzzy: mask.invert()

            # Add item to the list
            #item = (levels_dict, mask)
            #masks_levels.append(item)
            masks_levels[levels_dict] = mask

        # Return the masks
        return masks_levels

    # -----------------------------------------------------------------

    def create_mask_for_level(self, significance, level):

        """
        This function ...
        :param significance:
        :param level:
        :return:
        """

        # Create the mask
        mask = Mask(significance > level)

        # Only keep largest patch
        #mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

        # Fill holes
        #mask.fill_holes()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def create_fuzzy_mask_for_level(self, significance, level, fuzziness, offset=1.0):

        """
        This function ...
        :param significance:
        :param level:
        :param fuzziness:
        :param offset:
        :return:
        """

        # Determine the maximum significance
        max_significance = np.nanmax(significance)

        # Debugging
        log.debug("Maximal significance: " + str(max_significance))

        # Construct value range
        lower_relative = 1. - fuzziness  # example: 1. - 0.1
        upper_relative = 1. + fuzziness

        # ADAPT IF THE RANGE IS TOO BROAD
        if level * upper_relative > max_significance - offset:
            log.warning("Changing the upper relative sigma level for the fuzzy edge from " + str(upper_relative) + " to 1")
            upper_relative = 1.

        value_range = RealRange.around(level, lower_relative, upper_relative)

        # Debugging
        log.debug("Sigma level range for fuzzy edge: " + str(value_range))

        # Check maximum of the range
        if value_range.max + offset > max_significance: raise ValueError("This should not happen")

        # Create the mask
        mask = AlphaMask.between(significance, value_range)

        # Only keep largest patch
        #mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

        # Fill holes
        #mask.fill_holes()

        # Return the mask
        return mask

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
        return the_map.soften_edges(ellipse, softening_range)

    # -----------------------------------------------------------------

    @property
    def has_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_colour_maps

    # -----------------------------------------------------------------

    @property
    def colour_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.colour_has_methods

    # -----------------------------------------------------------------

    @property
    def colour_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.colour_map_methods

    # -----------------------------------------------------------------

    @property
    def colour_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.colour_map_names

    # -----------------------------------------------------------------

    def colour_map_names_for_method(self, method):

        """
        This function...
        :param method:
        :return:
        """

        return self.static_collection.colour_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def colour_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.colour_origins

    # -----------------------------------------------------------------

    @property
    def has_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_ssfr_maps

    # -----------------------------------------------------------------

    @property
    def ssfr_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ssfr_has_methods

    # -----------------------------------------------------------------

    @property
    def ssfr_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ssfr_map_methods

    # -----------------------------------------------------------------

    @property
    def ssfr_map_names_no_methods(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_collection.ssfr_map_names

    # -----------------------------------------------------------------

    def ssfr_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.static_collection.ssfr_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def ssfr_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ssfr_origins

    # -----------------------------------------------------------------

    @property
    def has_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_tir_maps

    # -----------------------------------------------------------------

    @property
    def tir_has_methods(self):

        """
        Thisf ucntion ...
        :return:
        """

        return self.static_collection.tir_has_methods

    # -----------------------------------------------------------------

    @property
    def tir_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.tir_map_methods

    # -----------------------------------------------------------------

    @property
    def tir_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.tir_map_names

    # -----------------------------------------------------------------

    def tir_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.static_collection.tir_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def tir_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.tir_origins

    # -----------------------------------------------------------------

    @property
    def has_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_attenuation_maps

    # -----------------------------------------------------------------

    @property
    def attenuation_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.attenuation_has_methods

    # -----------------------------------------------------------------

    @property
    def attenuation_map_methods(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_collection.attenuation_map_methods

    # -----------------------------------------------------------------

    @property
    def attenuation_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.attenuation_map_names

    # -----------------------------------------------------------------

    def attenuation_map_names_for_method(self, method):

        """
        Thisfunction ...
        :param method:
        :return:
        """

        return self.static_collection.attenuation_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def attenuation_origins(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_collection.attenuation_origins

    # -----------------------------------------------------------------

    @property
    def has_old_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_old_maps

    # -----------------------------------------------------------------

    @property
    def old_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.old_has_methods

    # -----------------------------------------------------------------

    @property
    def old_map_methods(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_collection.old_map_methods

    # -----------------------------------------------------------------

    @property
    def old_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.old_map_names

    # -----------------------------------------------------------------

    def old_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.static_collection.old_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def old_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.old_origins

    # -----------------------------------------------------------------

    @property
    def has_young_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_young_maps

    # -----------------------------------------------------------------

    @property
    def young_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.young_has_methods

    # -----------------------------------------------------------------

    @property
    def young_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.young_map_methods

    # -----------------------------------------------------------------

    @property
    def young_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.young_map_names

    # -----------------------------------------------------------------

    def young_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.static_collection.young_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def young_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.young_origins

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Static because maps collection will not change when in a selection component
        return self.static_collection.has_ionizing_maps

    # -----------------------------------------------------------------

    @property
    def ionizing_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ionizing_has_methods

    # -----------------------------------------------------------------

    @property
    def ionizing_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ionizing_map_methods

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ionizing_map_names

    # -----------------------------------------------------------------

    def ionizing_map_names_for_method(self, method):

        """
        This function ....
        :param method:
        :return:
        """

        return self.static_collection.ionizing_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def ionizing_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ionizing_origins

    # -----------------------------------------------------------------

    @property
    def has_dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_dust_maps

    # -----------------------------------------------------------------

    @property
    def dust_has_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.dust_has_methods

    # -----------------------------------------------------------------

    @property
    def dust_map_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.dust_map_methods

    # -----------------------------------------------------------------

    @property
    def dust_map_names_no_methods(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.dust_map_names

    # -----------------------------------------------------------------

    def dust_map_names_for_method(self, method):

        """
        This unction ...
        :param method:
        :return:
        """

        return self.static_collection.dust_map_names_for_method(method)

    # -----------------------------------------------------------------

    @property
    def dust_origins(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.dust_origins

# -----------------------------------------------------------------

def get_single_mask_path(output_path, name, level, fuzzy):

    """
    This function ...
    :param output_path:
    :param name:
    :param level:
    :param fuzzy:
    :return:
    """

    # Determine the path
    mask_base_path = fs.join(output_path, "mask__" + name.replace(" ", "_") + "__" + repr(level))
    if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
    else: mask_path = mask_base_path + ".fits"

    # Return the path
    return mask_path

# -----------------------------------------------------------------

def has_single_mask(output_path, name, level, fuzzy):

    """
    This function ...
    :param output_path:
    :param name:
    :param level:
    :param fuzzy:
    :return:
    """

    # Get the path
    mask_path = get_single_mask_path(output_path, name, level, fuzzy=fuzzy)

    # Check
    return fs.is_file(mask_path)

# -----------------------------------------------------------------

def load_single_mask(output_path, name, level, fuzzy):

    """
    This function ...
    :param output_path:
    :param name:
    :param level:
    :param fuzzy:
    :return:
    """

    # Get the path
    mask_path = get_single_mask_path(output_path, name, level, fuzzy=fuzzy)

    # Check
    if fuzzy: return AlphaMask.from_file(mask_path)
    else: return Mask.from_file(mask_path)

# -----------------------------------------------------------------

def get_combination_mask_path(output_path, levels, fuzzy):

    """
    This function ....
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine path
    mask_base_path = fs.join(output_path, "mask__combination__" + "_".join(repr(level) for level in levels))
    if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
    else: mask_path = mask_base_path + ".fits"

    # Return
    return mask_path

# -----------------------------------------------------------------

def has_combination_mask(output_path, levels, fuzzy):

    """
    This function ...
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine path
    mask_path = get_combination_mask_path(output_path, levels, fuzzy)

    # Check
    return fs.is_file(mask_path)

# -----------------------------------------------------------------

def load_combination_mask(output_path, levels, fuzzy):

    """
    This function ...
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine path
    mask_path = get_combination_mask_path(output_path, levels, fuzzy)

    # Load and return
    if fuzzy: return AlphaMask.from_file(mask_path)
    else: return Mask.from_file(mask_path)

# -----------------------------------------------------------------

def get_resulting_mask_path(output_path, levels, fuzzy):

    """
    This function ...
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine path
    mask_base_path = fs.join(output_path, "mask__result__" + "_".join(repr(level) for level in levels))
    if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
    else: mask_path = mask_base_path + ".fits"

    # Return
    return mask_path

# -----------------------------------------------------------------

def has_resulting_mask(output_path, levels, fuzzy):

    """
    This function ...
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine the path
    mask_path = get_resulting_mask_path(output_path, levels, fuzzy)

    # Check
    return fs.is_file(mask_path)

# -----------------------------------------------------------------

def load_resulting_mask(output_path, levels, fuzzy):

    """
    This function ...
    :param output_path:
    :param levels:
    :param fuzzy:
    :return:
    """

    # Determine the path
    mask_path = get_resulting_mask_path(output_path, levels, fuzzy)

    # Load and return
    if fuzzy: return AlphaMask.from_file(mask_path)
    else: return Mask.from_file(mask_path)

# -----------------------------------------------------------------
