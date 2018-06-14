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
import copy
from abc import ABCMeta

# Import the relevant PTS classes and modules
from .component import MapMakingComponent
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
from ...magic.basics.mask import MaskBase
from ...magic.core.alpha import AlphaMask, load_mask_or_alpha_mask
from ...core.basics.range import RealRange
from ...magic.tools import plotting
from ...core.basics import containers
from ...magic.region.region import PixelRegion, SkyRegion

# -----------------------------------------------------------------

plots_name = "plots"

# -----------------------------------------------------------------

map_plot_filename = "map"
edgeon_plot_filename = "edgeon"
deprojected_plot_filename = "deprojected"
deprojected_skirt_plot_filename = "deprojected_skirt"

# -----------------------------------------------------------------

class MapsSelectionComponent(MapMakingComponent):

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

        # Plots directories
        self.old_plots_path = None
        self.young_plots_path = None
        self.ionizing_plots_path = None
        self.dust_plots_path = None

        # The deprojected maps
        self.old_deprojected = dict()
        self.young_deprojected = dict()
        self.ionizing_deprojected = dict()
        self.dust_deprojected = dict()

        # The deprojected maps with SKIRT
        self.old_deprojected_skirt = dict()
        self.young_deprojected_skirt = dict()
        self.ionizing_deprojected_skirt = dict()
        self.dust_deprojected_skirt = dict()

        # Edgeon maps with SKIRT
        self.old_edgeon_skirt = dict()
        self.young_edgeon_skirt = dict()
        self.ionizing_edgeon_skirt = dict()
        self.dust_edgeon_skirt = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup of the base class
        super(MapsSelectionComponent, self).setup(**kwargs)

        # Plots directories
        self.old_plots_path = fs.create_directory_in(self.old_component_maps_path, plots_name)
        self.young_plots_path = fs.create_directory_in(self.young_component_maps_path, plots_name)
        self.ionizing_plots_path = fs.create_directory_in(self.ionizing_component_maps_path, plots_name)
        self.dust_plots_path = fs.create_directory_in(self.dust_component_maps_path, plots_name)

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
        ninfs = self.replace_infinities(the_map)

        # NaNs
        nnans = self.replace_nans(the_map)

        # Negatives
        nnegatives = self.replace_negatives(the_map)

        # Return
        return ninfs, nnans, nnegatives

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

        # Return the number of infinities
        return ninfs

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

        # Return the number of nans
        return nnans

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

        # Return the number of negatives
        return nnegatives

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

        # Return the number of leftover pixels in the map
        return np.sum(np.logical_not(mask))

    # -----------------------------------------------------------------

    def crop_map(self, the_map, factor, mask=None):

        """
        This function ...
        :param the_map:
        :param factor:
        :param mask:
        :return:
        """

        x_min, x_max, y_min, y_max = the_map.crop_to(self.truncation_box, factor=factor)
        if mask is not None: mask.crop(x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    def clip_map(self, the_map, origins, convolve=True, remote=None, npixels=1, connectivity=8,
                 rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None,
                 dilate=True, dilate_fuzzy=True, dilation_radius=20, dilation_nbins=20, soften=True,
                 softening_radius=15, softening_nbins=10, relative_softening_radius=None, boundary=None, plot=False,
                 relative_dilation_radius=None):

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
        :param dilate:
        :param dilate_fuzzy:
        :param dilation_radius:
        :param dilation_nbins:
        :param soften:
        :param softening_radius:
        :param softening_nbins:
        :param relative_softening_radius:
        :param boundary:
        :param plot:
        :param relative_dilation_radius:
        :return:
        """

        # Create the clip mask
        mask = self.get_clip_mask(origins, wcs=the_map.wcs, convolve=convolve, remote=remote, npixels=npixels,
                                  connectivity=connectivity, rebin_remote_threshold=rebin_remote_threshold,
                                  fuzzy=fuzzy, fuzziness=fuzziness, fuzziness_offset=fuzziness_offset,
                                  output_path=output_path, dilate=dilate, dilate_fuzzy=dilate_fuzzy,
                                  dilation_radius=dilation_radius, dilation_nbins=dilation_nbins, soften=soften,
                                  softening_radius=softening_radius, softening_nbins=softening_nbins,
                                  relative_softening_radius=relative_softening_radius, boundary=boundary, plot=plot,
                                  relative_dilation_radius=relative_dilation_radius)

        # Plot the mask
        if plot: plotting.plot_mask(mask, title="clipping mask")

        # Apply the mask
        if isinstance(mask, AlphaMask): the_map.apply_alpha_mask(mask)
        elif isinstance(mask, MaskBase): the_map.apply_mask(mask)
        else: raise ValueError("Invalid mask: must be regular mask or alpha mask")

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def make_clipped_maps(self, name, the_map, origins, levels_dict, convolve=True, remote=None, rebin_remote_threshold=None,
                          npixels=1, connectivity=8, present=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1.,
                          return_masks=False, current=None, current_masks=None, dilate=True, dilate_fuzzy=True,
                          dilation_radius=20, dilation_nbins=20, soften=True, softening_radius=15, softening_nbins=10,
                          resoften_current_masks=False, relative_softening_radius=None, boundary=None, plot=False,
                          relative_dilation_radius=None):

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
        :param present:
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
        :param return_masks:
        :param current:
        :param current_masks:
        :param dilate:
        :param dilate_fuzzy:
        :param dilation_radius:
        :param dilation_nbins:
        :param soften:
        :param softening_radius:
        :param softening_nbins:
        :param resoften_current_masks:
        :param relative_softening_radius:
        :param boundary:
        :param plot:
        :param relative_dilation_radius:
        :return:
        """

        # Check
        if current is not None and current_masks is not None:
            if not containers.equal_keys(current, current_masks): raise ValueError("levels dictionaries should be the same")

        # Make the masks
        masks = self.make_clip_masks(name, origins, levels_dict, wcs=the_map.wcs, convolve=convolve, remote=remote,
                                     rebin_remote_threshold=rebin_remote_threshold, npixels=npixels,
                                     connectivity=connectivity, present=present, fuzzy=fuzzy, fuzziness=fuzziness,
                                     fuzziness_offset=fuzziness_offset, current_masks=current_masks, dilate=dilate,
                                     dilate_fuzzy=dilate_fuzzy, dilation_radius=dilation_radius,
                                     dilation_nbins=dilation_nbins, soften=soften, softening_radius=softening_radius,
                                     softening_nbins=softening_nbins, resoften_current_masks=resoften_current_masks,
                                     relative_softening_radius=relative_softening_radius, boundary=boundary, plot=plot,
                                     relative_dilation_radius=relative_dilation_radius)

        # The maps
        if current is not None: maps = copy.copy(current) # shallow copy
        else: maps = dict()

        # Make the maps
        for levels in masks:

            # Already in the maps dictionary
            if levels in maps: continue

            # Debugging
            log.debug("Clipping the map for sigma levels [" + tostr(levels) + "] ...")

            # Get the mask
            mask = masks[levels]

            # Make clipped copy of the map
            clipped_map = the_map.copy()
            #if fuzzy: clipped_map.apply_alpha_mask(mask)
            #else: clipped_map.apply_mask(mask)
            # CAN ALSO BE 'FUZZY' when regular mask is softened
            if isinstance(mask, AlphaMask): clipped_map.apply_alpha_mask(mask)
            elif isinstance(mask, MaskBase): clipped_map.apply_mask(mask)
            else: raise ValueError("Invalid mask: must be regular mask or alpha mask")

            # Plot
            if plot: plotting.plot_frame(clipped_map, title="clipped '" + name + "' map")

            # Add to the dictionary
            maps[levels] = clipped_map

        # Return the dictionary
        if return_masks: return maps, masks
        else: return maps

    # -----------------------------------------------------------------

    def get_clip_mask(self, origins, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8,
                      rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None,
                      dilate=True, dilate_fuzzy=True, dilation_radius=20, dilation_nbins=20, soften=True,
                      softening_radius=15, softening_nbins=10, relative_softening_radius=None, boundary=None,
                      plot=False, relative_dilation_radius=None):

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
        :param dilate:
        :param dilate_fuzzy:
        :param dilation_radius:
        :param dilation_nbins:
        :param soften:
        :param softening_radius:
        :param softening_nbins:
        :param relative_softening_radius:
        :param boundary:
        :param plot:
        :param relative_dilation_radius:
        :return:
        """

        # Check if cached
        if wcs is None and tuple(origins) in self.clip_masks: return self.clip_masks[tuple(origins)]
        else:

            # Make the mask
            mask = self.make_clip_mask(origins, self.levels, wcs=wcs, convolve=convolve, remote=remote,
                                       npixels=npixels, connectivity=connectivity,
                                       rebin_remote_threshold=rebin_remote_threshold, fuzzy=fuzzy, fuzziness=fuzziness,
                                       fuzziness_offset=fuzziness_offset, output_path=output_path,
                                       dilate=dilate, dilate_fuzzy=dilate_fuzzy,
                                       dilation_radius=dilation_radius, dilation_nbins=dilation_nbins, soften=soften,
                                       softening_radius=softening_radius, softening_nbins=softening_nbins,
                                       relative_softening_radius=relative_softening_radius, boundary=boundary,
                                       plot=plot, relative_dilation_radius=relative_dilation_radius)

            # Cache the mask
            if wcs is None: self.clip_masks[tuple(origins)] = mask

            # Return the mask
            return mask

    # -----------------------------------------------------------------

    def make_clip_mask(self, origins, levels, wcs=None, convolve=True, remote=None, npixels=1, connectivity=8,
                       rebin_remote_threshold=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1., output_path=None,
                       dilate=True, dilate_fuzzy=True, dilation_radius=5, dilation_max_radius=20,
                       dilation_nbins=20, soften=True, softening_radius=15, softening_nbins=10,
                       relative_softening_radius=None, boundary=None, plot=False, relative_dilation_radius=None):

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
        :param dilate:
        :param dilate_fuzzy:
        :param dilation_radius:
        :param dilation_max_radius:
        :param dilation_nbins:
        :param soften:
        :param softening_radius:
        :param softening_nbins:
        :param relative_softening_radius:
        :param boundary:
        :param plot:
        :param relative_dilation_radius:
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
            if mask.wcs is None and wcs is None: raise ValueError("Coordinate system of mask not defined")
            elif wcs is None: wcs = mask.wcs

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

            # All can be ignored?
            ignore_all_frames = sequences.same_contents(ignore, frames.names)
            if ignore_all_frames: log.debug("All frames can be ignored")

            # Some frames should be convolved and rebinned
            if not ignore_all_frames:

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

                # Get the common WCS of all the frames
                first_not_ignored_name = sequences.find_first_not_in(frames.names, ignore)
                if wcs is None: wcs = frames[first_not_ignored_name].wcs #wcs = frames.wcs #wcs = frames[0].wcs

                # Debugging
                log.debug("Coordinate system of the frames:")
                if log.is_debug: print(wcs)

            # Set softening radius from relative softening radius
            if relative_softening_radius is not None:
                softening_radius = relative_softening_radius * float(wcs.xsize)
                # Debuggig
                log.debug("The softening radius is " + str(softening_radius) + " pixels")

            # Set dilation radius from relative dilation radius
            if relative_dilation_radius is not None:
                dilation_radius = relative_dilation_radius * float(wcs.xsize)
                # Debugging
                log.debug("The dilation radius is " + str(dilation_radius) + " pixels")

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

            # Plot
            if plot: plotting.plot_mask(mask, title="raw combination mask")

            # Save the mask
            if output_path is not None:
                mask_base_path = fs.join(output_path, "mask__combination__" + "_".join(repr(level) for level in levels_list))
                if fuzzy: mask_path = mask_base_path + "__fuzzy.fits"
                else: mask_path = mask_base_path + ".fits"
                mask.saveto(mask_path)

        # Create the boundary mask
        if boundary is not None:
            if isinstance(boundary, PixelRegion): boundary_mask = boundary.to_mask(frames[0].xsize, frames[0].ysize, invert=True)
            elif isinstance(boundary, SkyRegion): boundary_mask = boundary.to_pixel(wcs).to_mask(frames[0].xsize, frames[0].ysize, invert=True)
            elif isinstance(boundary, MaskBase): boundary_mask = boundary.inverse()
            else: raise ValueError("Invalid value for 'boundary': must be region or mask")
        else: boundary_mask = None

        # Only keep largest patch
        mask = mask.largest(npixels=npixels, connectivity=connectivity)

        # Plot
        if plot: plotting.plot_mask(mask, title="largest")

        # Alpha mask?
        if fuzzy:

            # First dilate if requested
            if dilate_fuzzy: mask.disk_dilate(radius=dilation_radius, nbins=dilation_nbins, max_radius=dilation_max_radius)

            # Remove parts outside boundary
            if boundary_mask is not None: mask.unmask(boundary_mask)

            # Now fill holes AFTER dilation, some holes may already be filled
            mask.fill_holes()

        # Regular mask?
        else:

            # First dilate if requested
            if dilate: mask.disk_dilate(radius=dilation_radius)

            # Plot
            if plot: plotting.plot_mask(mask, title="dilated")

            # Remove parts outside boundary
            if boundary_mask is not None: mask.unmask(boundary_mask)

            # Plot
            if plot: plotting.plot_mask(mask, title="after boundary cut-off")

            # First fill holes
            mask.fill_holes(connectivity=4)

            # Now create softened (fuzzy) mask based on dilation if requested
            if soften:
                mask = mask.softened(softening_radius, nbins=softening_nbins)
                if plot: plotting.plot_mask(mask, title="after softening")

            # Invert FOR NORMAL MASKS: WE HAVE TO SET PIXELS TO ZERO THAT ARE NOT ON THE MASK
            else: mask.invert()

        # Set the mask's coordinate system
        mask.wcs = wcs

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
                        npixels=1, connectivity=8, present=None, fuzzy=False, fuzziness=0.5, fuzziness_offset=1.,
                        current_masks=None, dilate=True, dilate_fuzzy=True, dilation_radius=5, dilation_max_radius=20,
                        dilation_nbins=20, soften=True, softening_radius=20, softening_nbins=10, resoften_current_masks=False,
                        relative_softening_radius=None, boundary=None, plot=False, relative_dilation_radius=None):

        """
        Thisn function ...
        :param name:
        :param origins:
        :param levels_dict:
        :param wcs:
        :param convolve:
        :param remote:
        :param rebin_remote_threshold:
        :param npixels:
        :param connectivity:
        :param present:
        :param fuzzy:
        :param fuzziness:
        :param fuzziness_offset:
        :param current_masks:
        :param dilate:
        :param dilate_fuzzy:
        :param dilation_radius:
        :param dilation_max_radius:
        :param dilation_nbins:
        :param soften:
        :param softening_radius:
        :param softening_nbins:
        :param resoften_current_masks:
        :param relative_softening_radius:
        :param boundary:
        :param plot:
        :param relative_dilation_radius:
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

        # Get the common WCS of all the frames
        if wcs is None: wcs = frames[0].wcs

        # Create the boundary mask
        if boundary is not None:
            if isinstance(boundary, PixelRegion): boundary_mask = boundary.to_mask(frames[0].xsize, frames[0].ysize, invert=True)
            elif isinstance(boundary, SkyRegion): boundary_mask = boundary.to_pixel(wcs).to_mask(frames[0].xsize, frames[0].ysize, invert=True)
            elif isinstance(boundary, MaskBase): boundary_mask = boundary.inverse()
            else: raise ValueError("Invalid value for 'boundary': must be region or mask")
        else: boundary_mask = None

        # Set softening radius from relative softening radius
        if relative_softening_radius is not None:
            softening_radius = relative_softening_radius * float(wcs.xsize)
            # Debuggig
            log.debug("The softening radius is " + str(softening_radius) + " pixels")

        # Set dilation radius from relative dilation radius
        if relative_dilation_radius is not None:
            dilation_radius = relative_dilation_radius * float(wcs.xsize)
            # Debugging
            log.debug("The dilation radius is " + str(dilation_radius) + " pixels")

        # Get the number of frames
        names = frames.names
        nframes = len(names)

        # Get list of level lists, sorted on the order of the frames
        frame_levels = []

        # LEVELS ARE GIVEN AS A SEQUENCE: SAME FOR EACH FRAME
        if types.is_sequence(levels_dict):
            for _ in range(nframes): frame_levels.append(levels_dict)

        # LEVELS ARE GIVEN AS A DICTIONARY: DEFINED FOR EACH IMAGE NAME
        elif types.is_dictionary(levels_dict):
            for name in names: frame_levels.append(levels_dict[name]) # name = origin

        # Invalid
        else: raise ValueError("Levels must be specified as list (same for each image) or dictionary keyed on filter")

        # Container to keep the set of masks
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

            # Check if already created
            if current_masks is not None and levels_dict in current_masks:

                # Get the mask
                mask = current_masks[levels_dict]
                #plotting.plot_mask(mask, title="current mask")

                # Re-soften?
                if resoften_current_masks:

                    # Get actual mask
                    if isinstance(mask, AlphaMask): mask = mask.opaque
                    elif isinstance(mask, MaskBase): mask = mask.inverse()
                    else: raise ValueError("Invalid mask object")
                    #plotting.plot_mask(mask, title="opaque or inverse")

                    # First dilate if requested
                    if dilate: mask.disk_dilate(radius=dilation_radius)

                    # Remove parts outside boundary
                    if boundary_mask is not None: mask.unmask(boundary_mask)

                    # Fill holes first
                    mask.fill_holes(connectivity=4)
                    #plotting.plot_mask(mask, title="filled holes")

                    # Soften the mask
                    mask = mask.softened(softening_radius, nbins=softening_nbins)
                    #plotting.plot_mask(mask, title="softened")

                    # Set the mask's coordinate system
                    mask.wcs = wcs

                # Set the mask and continue
                masks_levels[levels_dict] = mask
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

            # Plot
            if plot: plotting.plot_mask(mask, title="largest")

            # Alpha mask?
            if fuzzy:

                # First dilate if requested
                if dilate_fuzzy: mask.disk_dilate(radius=dilation_radius, nbins=dilation_nbins, max_radius=dilation_max_radius)

                # Remove parts outside boundary
                if boundary_mask is not None: mask.unmask(boundary_mask)

                # Now fill holes AFTER dilation, some holes may already be filled
                mask.fill_holes()

            # Regular mask?
            else:

                # First dilate if requested
                if dilate: mask.disk_dilate(radius=dilation_radius)

                # Plot
                if plot: plotting.plot_mask(mask, title="dilated")

                # Remove parts outside boundary
                if boundary_mask is not None: mask.unmask(boundary_mask)

                # Plot
                if plot: plotting.plot_mask(mask, title="after boundary cut-off")

                # First fill holes
                mask.fill_holes(connectivity=4)

                # Now create softened (fuzzy) mask based on dilation if requested
                if soften:
                    mask = mask.softened(softening_radius, nbins=softening_nbins)
                    if plot: plotting.plot_mask(mask, title="after softening")

                # Invert FOR NORMAL MASKS: WE HAVE TO SET PIXELS TO ZERO THAT ARE NOT ON THE MASK
                else: mask.invert()

            # Set the mask's coordinate system
            mask.wcs = wcs

            # Add the mask to the dictionary
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

    def old_plotting_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.old_plots_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def young_plotting_path_for_map(self, name, create=True):

        """
        Thisf unction ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.young_plots_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def ionizing_plotting_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.ionizing_plots_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def dust_plotting_path_for_map(self, name, create=True):

        """
        Thisf unction ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.dust_plots_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def plot_components_maps(self, format="pdf"):

        """
        Thisf unction ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # Old
        self.plot_old_component_maps(format=format)

        # Young
        self.plot_young_component_maps(format=format)

        # Ionizing
        self.plot_ionizing_component_maps(format=format)

        # Dust
        self.plot_dust_component_maps(format=format)

    # -----------------------------------------------------------------

    @property
    def old_scale(self):

        """
        Thisf unction ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    @property
    def old_cmap(self):

        """
        This function ...
        :return:
        """

        return "afmhot"

    # -----------------------------------------------------------------

    @property
    def old_color(self):

        """
        This function ...
        :return:
        """

        return "orange"

    # -----------------------------------------------------------------

    def old_map_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.old_plotting_path_for_map(name), map_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_old_component_maps(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting old stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Plotting the '" + name + "' old stellar component map ...")

            # Determine path
            path = self.old_map_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.old_maps[name], path=path, format=format, interval=minmax,
                                             scale=self.old_scale, cmap=self.old_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent,
                                             truncate_outside=self.truncation_ellipse)

    # -----------------------------------------------------------------

    @property
    def young_scale(self):

        """
        This function ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    @property
    def young_cmap(self):

        """
        This function ...
        :return:
        """

        return "plasma"

    # -----------------------------------------------------------------

    @property
    def young_color(self):

        """
        This function ...
        :return:
        """

        return "lime_green"

    # -----------------------------------------------------------------

    def young_map_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.young_plotting_path_for_map(name), map_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_young_component_maps(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting young stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Plotting the '" + name + "' young stellar component map ...")

            # Determine path
            path = self.young_map_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.young_maps[name], path=path, format=format, interval=minmax,
                                             scale=self.young_scale, cmap=self.young_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent,
                                             truncate_outside=self.truncation_ellipse)

    # -----------------------------------------------------------------

    @property
    def ionizing_scale(self):

        """
        Thisfunction ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    @property
    def ionizing_cmap(self):

        """
        This function ...
        :return:
        """

        return "summer"

    # -----------------------------------------------------------------

    @property
    def ionizing_color(self):

        """
        Thisnfunction ...
        :return:
        """

        #return "sky_blue"
        return "turquoise"

    # -----------------------------------------------------------------

    def ionizing_map_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.ionizing_plotting_path_for_map(name), map_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_ionizing_component_maps(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting ionizing stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Plotting the '" + name + "' ionizing stellar component map ...")

            # Determine path
            path = self.ionizing_map_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.ionizing_maps[name], path=path, format=format, interval=minmax,
                                             scale=self.ionizing_scale, cmap=self.ionizing_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent,
                                             truncate_outside=self.truncation_ellipse)

    # -----------------------------------------------------------------

    @property
    def dust_scale(self):

        """
        Thisf unction ...
        :return:
        """

        return "linear"

    # -----------------------------------------------------------------

    @property
    def dust_cmap(self):

        """
        Thisf unction ...
        :return:
        """

        return "inferno"

    # -----------------------------------------------------------------

    @property
    def dust_color(self):

        """
        Thisnf unction ...
        :return:
        """

        return "red"

    # -----------------------------------------------------------------

    def dust_map_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.dust_plotting_path_for_map(name), map_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_dust_component_maps(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        Thisf unction ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting dust component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Plotting the '" + name + "' dust component map ...")

            # Determine path
            path = self.dust_map_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.dust_maps[name], path=path, format=format, interval=minmax,
                                             scale=self.dust_scale, cmap=self.dust_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent,
                                             truncate_outside=self.truncation_ellipse)

    # -----------------------------------------------------------------

    def plot_components_radial_profiles(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the radial profiles of the maps ...")

        # Old
        self.plot_old_component_radial_profiles(format=format)

        # Young
        self.plot_young_component_radial_profiles(format=format)

        # Ionizing
        self.plot_ionizing_component_radial_profiles(format=format)

        # Dust
        self.plot_dust_component_radial_profiles(format=format)

    # -----------------------------------------------------------------

    def plot_old_component_radial_profiles(self, format="pdf"):

        """
        Thisf unction ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting radial profiles of old stellar component maps ...")

    # -----------------------------------------------------------------

    def plot_young_component_radial_profiles(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting radial profiles of young stellar component maps ...")

    # -----------------------------------------------------------------

    def plot_ionizing_component_radial_profiles(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting radial profiles of ionizing stellar component maps ...")

    # -----------------------------------------------------------------

    def plot_dust_component_radial_profiles(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting radial profiles of dust component maps ...")

    # -----------------------------------------------------------------

    def plot_components_masks(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the component map masks ...")

        # Old
        self.plot_old_component_masks(format=format)

        # Young
        self.plot_young_component_masks(format=format)

        # Ionizing
        self.plot_ionizing_component_masks(format=format)

        # Dust
        self.plot_dust_component_masks(format=format)

    # -----------------------------------------------------------------

    def plot_old_component_masks(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting old stellar component map masks ...")

    # -----------------------------------------------------------------

    def plot_young_component_masks(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting young stellar component map masks ...")

    # -----------------------------------------------------------------

    def plot_ionizing_component_masks(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting ionizing stellar component map masks ...")

    # -----------------------------------------------------------------

    def plot_dust_component_masks(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting dust component map masks ...")

    # -----------------------------------------------------------------

    def plot_components_deprojected(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected component maps ...")

        # Old
        self.plot_old_deprojected(format=format)

        # Young
        self.plot_young_deprojected(format=format)

        # Ionizing
        self.plot_ionizing_deprojected(format=format)

        # Dust
        self.plot_dust_deprojected(format=format)

    # -----------------------------------------------------------------

    def old_deprojected_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.old_plotting_path_for_map(name), deprojected_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_old_deprojected(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting deprojected old stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.old_deprojected:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected old stellar component map ...")

            # Determine path
            path = self.old_deprojected_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.old_deprojected[name], path=path, format=format, interval=minmax,
                                             scale=self.old_scale, cmap=self.old_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def young_deprojected_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.young_plotting_path_for_map(name), deprojected_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_young_deprojected(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting deprojected young stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.young_deprojected:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected young stellar component map ...")

            # Determine path
            path = self.young_deprojected_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.young_deprojected[name], path=path, format=format, interval=minmax,
                                             scale=self.young_scale, cmap=self.young_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def ionizing_deprojected_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.ionizing_plotting_path_for_map(name), deprojected_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_ionizing_deprojected(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting deprojected ionizing stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.ionizing_deprojected:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected ionizing stellar component map ...")

            # Determine path
            path = self.ionizing_deprojected_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.ionizing_deprojected[name], path=path, format=format, interval=minmax,
                                             scale=self.ionizing_scale, cmap=self.ionizing_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def dust_deprojected_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.dust_plotting_path_for_map(name), deprojected_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_dust_deprojected(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting deprojected dust component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.dust_deprojected:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected dust component map ...")

            # Determine path
            path = self.dust_deprojected_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.dust_deprojected[name], path=path, format=format, interval=minmax,
                                             scale=self.dust_scale, cmap=self.dust_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def plot_components_deprojected_skirt(self, format="pdf"):

        """
        This function ...
        :parma format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected with SKIRT component maps ...")

        # Old
        self.plot_old_deprojected_skirt(format=format)

        # Young
        self.plot_young_deprojected_skirt(format=format)

        # Ionizing
        self.plot_ionizing_deprojected_skirt(format=format)

        # Dust
        self.plot_dust_deprojected_skirt(format=format)

    # -----------------------------------------------------------------

    def old_deprojected_skirt_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.old_plotting_path_for_map(name), deprojected_skirt_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_old_deprojected_skirt(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        Thisf unction ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected with SKIRT old stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.old_deprojected_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected with SKIRT old stellar component map ...")

            # Determine path
            path = self.old_deprojected_skirt_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.old_deprojected_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.old_scale, cmap=self.old_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def young_deprojected_skirt_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.young_plotting_path_for_map(name), deprojected_skirt_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_young_deprojected_skirt(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        Thisf unction ...
        :param show_axes:
        :param transparent
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected with SKIRT young stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.young_deprojected_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected with SKIRT young stellar component map ...")

            # Determine path
            path = self.young_deprojected_skirt_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.young_deprojected_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.young_scale, cmap=self.young_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def ionizing_deprojected_skirt_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.ionizing_plotting_path_for_map(name), deprojected_skirt_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_ionizing_deprojected_skirt(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected with SKIRT ionizing stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.ionizing_deprojected_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected with SKIRT ionizing stellar component map ...")

            # Determine path
            path = self.ionizing_deprojected_skirt_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.ionizing_deprojected_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.ionizing_scale, cmap=self.ionizing_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def dust_deprojected_skirt_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.dust_plotting_path_for_map(name), deprojected_skirt_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_dust_deprojected_skirt(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        Thisf ucntion ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the deprojected with SKIRT dust component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.dust_deprojected_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' deprojected with SKIRT dust component map ...")

            # Determine path
            path = self.dust_deprojected_skirt_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.dust_deprojected_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.dust_scale, cmap=self.dust_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def plot_components_edgeon(self, format="pdf"):

        """
        This function ...
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the edgeon component maps ...")

        # Old
        self.plot_old_edgeon(format=format)

        # Young
        self.plot_young_edgeon(format=format)

        # Ionizing
        self.plot_ionizing_edgeon(format=format)

        # Dust
        self.plot_dust_edgeon(format=format)

    # -----------------------------------------------------------------

    def old_edgeon_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.old_plotting_path_for_map(name), edgeon_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_old_edgeon(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting edge-on old stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.old_edgeon_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' edge-on old stellar component map ...")

            # Determine path
            path = self.old_edgeon_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.old_edgeon_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.old_scale, cmap=self.old_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def young_edgeon_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.young_plotting_path_for_map(name), edgeon_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_young_edgeon(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting edge-on young stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.young_edgeon_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' edge-on young stellar component map ...")

            # Determine path
            path = self.young_edgeon_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.young_edgeon_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.young_scale, cmap=self.young_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def ionizing_edgeon_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.ionizing_plotting_path_for_map(name), edgeon_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_ionizing_edgeon(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        This function ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting edge-on ionizing stellar component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.ionizing_edgeon_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' edge-on ionizing stellar component map ...")

            # Determine path
            path = self.ionizing_edgeon_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.ionizing_edgeon_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.ionizing_scale, cmap=self.ionizing_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def dust_edgeon_plot_for_map(self, name, format="pdf"):

        """
        This function ...
        :param name:
        :param format:
        :return:
        """

        return fs.join(self.dust_plotting_path_for_map(name), edgeon_plot_filename + "." + format)

    # -----------------------------------------------------------------

    def plot_dust_edgeon(self, show_axes=False, transparent=True, colorbar=True, format="pdf"):

        """
        Thisf unction ...
        :param show_axes:
        :param transparent:
        :param colorbar:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting edge-on dust component maps ...")

        # Interval
        minmax = "pts"

        # Loop over the maps
        for name in self.dust_edgeon_skirt:

            # Debugging
            log.debug("Plotting the '" + name + "' edge-on dust component map ...")

            # Determine path
            path = self.dust_edgeon_plot_for_map(name, format=format)

            # Plot the map
            vmin, vmax = plotting.plot_frame(self.dust_edgeon_skirt[name], path=path, format=format, interval=minmax,
                                             scale=self.dust_scale, cmap=self.dust_cmap, colorbar=colorbar,
                                             show_axes=show_axes, transparent=transparent)

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
    #if fuzzy: return AlphaMask.from_file(mask_path, no_wcs=True)
    #else: return Mask.from_file(mask_path, no_wcs=True)
    return load_mask_or_alpha_mask(mask_path, no_wcs=True)

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
    #if fuzzy: return AlphaMask.from_file(mask_path)
    #else: return Mask.from_file(mask_path)
    return load_mask_or_alpha_mask(mask_path)

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
    #if fuzzy: return AlphaMask.from_file(mask_path)
    #else: return Mask.from_file(mask_path)
    return load_mask_or_alpha_mask(mask_path)

# -----------------------------------------------------------------
