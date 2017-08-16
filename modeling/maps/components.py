#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.components Contains the ComponentMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...core.basics.configuration import prompt_string_list, prompt_real
from ...core.tools.utils import lazyproperty
from ...magic.core.frame import Frame
from ...core.tools import filesystem as fs
from ...core.tools import numbers
from ...core.basics.range import RealRange
from ...core.tools.stringify import tostr
from ...magic.core.mask import intersection

# -----------------------------------------------------------------

masks_name = "masks"
deprojected_name = "deprojected"

# -----------------------------------------------------------------

class ComponentMapsMaker(MapsComponent):

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
        super(ComponentMapsMaker, self).__init__(*args, **kwargs)

        # Paths
        self.old_masks_path = None
        self.young_masks_path = None
        self.ionizing_masks_path = None
        self.dust_masks_path = None
        self.old_deprojection_path = None
        self.young_deprojection_path = None
        self.ionizing_deprojection_path = None
        self.dust_deprojection_path = None

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

        # The clip masks
        self.old_masks = dict()
        self.young_masks = dict()
        self.ionizing_masks = dict()
        self.dust_masks = dict()

        # The deprojected maps
        self.old_deprojected = dict()
        self.young_deprojected = dict()
        self.ionizing_deprojected = dict()
        self.dust_deprojected = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Prompt
        self.prompt()

        # 3. Load the maps
        self.load_maps()

        # 4. Process the maps
        self.process_maps()

        # 5. Deproject the maps
        self.deproject()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsMaker, self).setup(**kwargs)

        # Masks directories
        self.old_masks_path = fs.create_directory_in(self.old_component_maps_path, masks_name)
        self.young_masks_path = fs.create_directory_in(self.young_component_maps_path, masks_name)
        self.ionizing_masks_path = fs.create_directory_in(self.ionizing_component_maps_path, masks_name)
        self.dust_masks_path = fs.create_directory_in(self.dust_component_maps_path, masks_name)

        # Deprojected directories
        self.old_deprojection_path = fs.create_directory_in(self.old_component_maps_path, deprojected_name)
        self.young_deprojection_path = fs.create_directory_in(self.young_component_maps_path, deprojected_name)
        self.ionizing_deprojection_path = fs.create_directory_in(self.ionizing_component_maps_path, deprojected_name)
        self.dust_deprojection_path = fs.create_directory_in(self.dust_component_maps_path, deprojected_name)

        # Check
        if self.config.old is not None and self.config.not_old is not None: raise ValueError("Cannot specify both 'old' and 'not_old'")
        if self.config.young is not None and self.config.not_young is not None: raise ValueError("Cannot specify both 'young' and 'not_young'")
        if self.config.ionizing is not None and self.config.not_ionizing is not None: raise ValueError("Cannot specify both 'ionizing' and 'not_ionizing'")
        if self.config.dust is not None and self.config.not_dust is not None: raise ValueError("Cannot specifiy both 'dust' and 'not_dust'")

        # Set selections
        if self.config.old is not None: self.old_selection = self.config.old
        if self.config.young is not None: self.young_selection = self.config.young
        if self.config.ionizing is not None: self.ionizing_selection = self.config.ionizing
        if self.config.dust is not None: self.dust_selection = self.config.dust

        # All maps?
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # Check
        if self.config.all_old and self.config.not_old is not None: raise ValueError("Cannot specify 'not_old' with 'all_old' enabled")
        if self.config.all_young and self.config.not_young is not None: raise ValueError("Cannot specify 'not_young' with 'all_young' enabled")
        if self.config.all_ionizing and self.config.not_ionizing is not None: raise ValueError("Cannot specify 'not_ionizing' with 'all_ionizing' enabled")
        if self.config.all_dust and self.config.not_dust is not None: raise ValueError("Cannot specify 'not_dust' with 'all_dust' enabled")

        # All
        if self.config.all_old: self.old_selection = self.old_map_names
        if self.config.all_young: self.young_selection = self.young_map_names
        if self.config.all_ionizing: self.ionizing_selection = self.ionizing_map_names
        if self.config.all_dust: self.dust_selection = self.dust_map_names

        # Levels
        if self.config.levels is not None: self.levels = self.config.levels

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

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for user input ...")

        # Old
        if not self.has_old_selection: self.prompt_old()

        # Young
        if not self.has_young_selection: self.prompt_young()

        # Ionizing
        if not self.has_ionizing_selection: self.prompt_ionizing()

        # Dust
        if not self.has_dust_selection: self.prompt_dust()

        # Levels
        if not self.has_levels: self.prompt_levels()

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

    @lazyproperty
    def old_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.old_selection: origins.update(self.old_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.young_selection: origins.update(self.young_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.ionizing_selection: origins.update(self.ionizing_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.dust_selection: origins.update(self.dust_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        origins.update(self.old_selection_origins)
        origins.update(self.young_selection_origins)
        origins.update(self.ionizing_selection_origins)
        origins.update(self.dust_selection_origins)
        return list(origins)

    # -----------------------------------------------------------------

    def prompt_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the significance levels ...")

        # Initialize the levels dictionary
        self.levels = dict()

        # Loop over the images (origins)
        for fltr in self.selection_origins:

            # Set filter name
            filter_name = str(fltr)
            id = filter_name.lower().replace(" ", "_")

            # Prompt for the level
            level = prompt_real(id + "_level", "sigma level for the " + filter_name + " image", default=self.config.default_level)

            # Set the level
            self.levels[fltr] = level

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected maps ...")

        # Load the old stellar maps
        self.load_old()

        # Load the young stellar maps
        self.load_young()

        # Load the ionizing stellar maps
        self.load_ionizing()

        # Load the dust maps
        self.load_dust()

    # -----------------------------------------------------------------

    def load_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected old stellar maps ...")

        # Load
        for name in self.old_selection:

            # Debugging
            log.debug("Loading the '" + name + "' old stellar map ...")

            # Load
            self.old_maps[name] = Frame.from_file(self.old_map_paths[name])

    # -----------------------------------------------------------------

    def load_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected young stellar maps ...")

        # Load
        for name in self.young_selection:

            # Debugging
            log.debug("Loading the '" + name + "' young stellar map ...")

            # Load
            self.young_maps[name] = Frame.from_file(self.young_map_paths[name])

    # -----------------------------------------------------------------

    def load_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected ionizing stellar maps ...")

        # Load
        for name in self.ionizing_selection:

            # Debugging
            log.debug("Loading the '" + name + "' ionizing stellar map ...")

            # Load
            self.ionizing_maps[name] = Frame.from_file(self.ionizing_map_paths[name])

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected dust maps ...")

        # Load
        for name in self.dust_selection:

            # Debugging
            log.debug("Loading the '" + name + "' dust map ...")

            # Load
            self.dust_maps[name] = Frame.from_file(self.dust_map_paths[name])

    # -----------------------------------------------------------------

    def process_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the maps ...")

        # 1. Truncate
        self.truncate_maps()

        # 2. Clip
        self.clip_maps()

        # 3. Soften the edges
        self.soften_edges()

    # -----------------------------------------------------------------

    def truncate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the maps ...")

        # Old
        self.truncate_old_maps()

        # Young
        self.truncate_young_maps()

        # Ionizing
        self.truncate_ionizing_maps()

        # Dust
        self.truncate_dust_maps()

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

    def truncate_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Truncating the '" + name + "' old stellar map ...")

            # Get the map
            old_map = self.old_maps[name]

            # Get the truncation mask
            mask = self.get_truncation_mask_for_map(old_map)

            # Mask
            old_map[mask] = 0.0

    # -----------------------------------------------------------------

    def truncate_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Truncating the '" + name + "' young stellar map ...")

            # Get the map
            young_map = self.young_maps[name]

            # Get the truncation mask
            mask = self.get_truncation_mask_for_map(young_map)

            # Mask
            young_map[mask] = 0.0

    # -----------------------------------------------------------------

    def truncate_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Truncating the '" + name + "' ionizing stellar map ...")

            # Get the map
            ionizing_map = self.ionizing_maps[name]

            # Get the truncation mask
            mask = self.get_truncation_mask_for_map(ionizing_map)

            # Mask
            ionizing_map[mask] = 0.0

    # -----------------------------------------------------------------

    def truncate_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Truncating the '" + name + "' dust map ...")

            # Get the map
            dust_map = self.dust_maps[name]

            # Get the truncation mask
            mask = self.get_truncation_mask_for_map(dust_map)

            # mask
            dust_map[mask] = 0.0

    # -----------------------------------------------------------------

    def clip_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the maps ...")

        # Old
        self.clip_old_maps()

        # Young
        self.clip_young_maps()

        # Ionizing
        self.clip_ionizing_maps()

        # Dust
        self.clip_dust_maps()

    # -----------------------------------------------------------------

    def make_clip_mask(self, origins):

        """
        This function ...
        :param origins:
        :return:
        """

        # Debugging
        log.debug("Making a clip mask for the filters: " + tostr(origins) + " ...")

        # Get frame list
        frames = self.dataset.get_framelist_for_filters(origins)

        # Get error map list
        errors = self.dataset.get_errormaplist_for_filters(origins)

        # Convolve
        frames.convolve_to_highest_fwhm()
        errors.convolve_to_highest_fwhm()

        # Rebin the frames to the same pixelgrid
        frames.rebin_to_highest_pixelscale()
        errors.rebin_to_highest_pixelscale()

        # Convert the frames to the same unit
        frames.convert_to_same_unit(unit="Jy")
        errors.convert_to_same_unit(unit="Jy")

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
            level = self.levels[frame.filter]
            mask = frame > level * errormap
            masks.append(mask)

        # Determine name for the mask
        #mask_name = "_".join(combination_names)

        # Create intersection mask
        mask = intersection(*masks)

        # Only keep central
        mask = mask.central()

        # Fill holes
        mask.fill_holes()

        # Invert
        mask.invert()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def clip_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Clipping the '" + name + "' old stellar map ...")

            # Get the origins
            origins = self.old_map_origins[name]

            # Create the clip mask
            mask = self.make_clip_mask(origins)

            # Set the mask
            self.old_masks[name] = mask

            # Clip
            self.old_maps[name][mask] = 0.0

    # -----------------------------------------------------------------

    def clip_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            origins = self.young_map_origins[name]

            # Create the clip mask
            mask = self.make_clip_mask(origins)

            # Set the mask
            self.young_masks[name] = mask

            # Clip
            self.young_maps[name][mask] = 0.0

    # -----------------------------------------------------------------

    def clip_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Create the clip mask
            mask = self.make_clip_mask(origins)

            # Set the mask
            self.ionizing_masks[name] = mask

            # Clip
            self.ionizing_maps[name][mask] = 0.0

    # -----------------------------------------------------------------

    def clip_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            origins = self.dust_map_origins[name]

            # Create the clip mask
            mask = self.make_clip_mask(origins)

            # Set the mask
            self.dust_masks[name] = mask

            # Clip
            self.dust_maps[name][mask] = 0.0

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.softening_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_radius(self):

        """
        This function ...
        :return:
        """

        return numbers.geometric_mean(self.config.softening_start, 1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(self.config.softening_start / self.softening_radius, 1. / self.softening_radius)

    # -----------------------------------------------------------------

    def soften_edges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the map edges ...")

        # Old
        self.soften_edges_old()

        # Young
        self.soften_edges_young()

        # Ionizing
        self.soften_edges_ionizing()

        # Dust
        self.soften_edges_dust()

    # -----------------------------------------------------------------

    def soften_edges_old(self):

        """
        Thisn function ...
        :return:
        """

        # Infomr the user
        log.info("Softening the edges of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Softening the edges of the '" + name + "' old stellar map ...")

            # Get ellipse
            ellipse = self.softening_ellipse.to_pixel(self.old_maps[name].wcs)

            # Soften edges
            self.old_maps[name].soften_edges(ellipse, self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Softening the edges of the '" + name + "' young stellar map ...")

            # Get ellipse
            ellipse = self.softening_ellipse.to_pixel(self.young_maps[name].wcs)

            # Soften edges
            self.young_maps[name].soften_edges(ellipse, self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Softening the edges of the '" + name + "' ionizing stellar map ...")

            # Get ellipse
            ellipse = self.softening_ellipse.to_pixel(self.ionizing_maps[name].wcs)

            # Soften edges
            self.ionizing_maps[name].soften_edges(ellipse, self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Softening the edges of the '" + name + "' dust map ...")

            # Get ellipse
            ellipse = self.softening_ellipse.to_pixel(self.dust_maps[name].wcs)

            # Soften edges
            self.dust_maps[name].soften_edges(ellipse, self.softening_range)

    # -----------------------------------------------------------------

    def deproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the maps ...")

        # Old
        self.deproject_old()

        # Young
        self.deproject_young()

        # Ionizing
        self.deproject_ionizing()

        # Dust
        self.deproject_dust()

    # -----------------------------------------------------------------

    def deproject_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the old stellar maps ...")

    # -----------------------------------------------------------------

    def deproject_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the young stellar maps ...")

    # -----------------------------------------------------------------

    def deproject_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the ionizing stellar maps ...")

    # -----------------------------------------------------------------

    def deproject_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the dust maps ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the clip masks
        self.write_masks()

        # Write the deprojected maps
        self.write_deprojected()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Old
        self.write_old_maps()

        # Young
        self.write_young_maps()

        # Ionizing
        self.write_ionizing_maps()

        # Dust
        self.write_dust_maps()

    # -----------------------------------------------------------------

    def write_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar map ...")

            # Determine the path
            path = fs.join(self.old_component_maps_path, name + ".fits")

            # Save the map
            self.old_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar map ...")

            # Determine the path
            path = fs.join(self.young_component_maps_path, name + ".fits")

            # Save the map
            self.young_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar map ...")

            # Determine the path
            path = fs.join(self.ionizing_component_maps_path, name + ".fits")

            # Save the map
            self.ionizing_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Writing the '" + name + "' dust map ...")

            # Determine the path
            path = fs.join(self.dust_component_maps_path, name + ".fits")

            # Save the map
            self.dust_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks ...")

        # Old
        self.write_old_masks()

        # Young
        self.write_young_masks()

        # Ionizing
        self.write_ionizing_masks()

        # Dust
        self.write_dust_masks()

    # -----------------------------------------------------------------

    def write_old_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the old stellar maps ...")

        # Loop over the masks
        for name in self.old_masks:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar mask ...")

            # Determine the path
            path = fs.join(self.old_masks_path, name + ".fits")

            # Write
            self.old_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the young stellar maps ...")

        # Loop over the masks
        for name in self.young_masks:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar mask ...")

            # Determine the path
            path = fs.join(self.young_masks_path, name + ".fits")

            # Write
            self.young_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the ionizing stellar maps ...")

        # Loop over the masks
        for name in self.ionizing_masks:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar mask ...")

            # Determine the path
            path = fs.join(self.ionizing_masks_path, name + ".fits")

            # Write
            self.ionizing_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the dust maps ...")

        # Loop over the masks
        for name in self.dust_masks:

            # Debugging
            log.debug("Writing the '" + name + "' dust mask ...")

            # Determine the path
            path = fs.join(self.dust_masks_path, name + ".fits")

            # Write
            self.dust_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected maps ...")

        # Old
        self.write_old_deprojected()

        # Young
        self.write_young_deprojected()

        # Ionizing
        self.write_ionizing_deprojected()

        # Dust
        self.write_dust_deprojected()

    # -----------------------------------------------------------------

    def write_old_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected old stellar maps ...")

        # Loop over the maps
        for name in self.old_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected old stellar map ...")

            # Determine path
            path = fs.join(self.old_deprojection_path, name + ".fits")

            # Write
            self.old_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected young stellar maps ...")

        # Loop over the maps
        for name in self.young_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected young stellar map ...")

            # Determine path
            path = fs.join(self.young_deprojection_path, name + ".fits")

            # Write
            self.young_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected ionizing stellar map ...")

            # Determine the path
            path = fs.join(self.ionizing_deprojection_path, name + ".fits")

            # Write
            self.ionizing_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected dust maps ...")

        # Loop over the maps
        for name in self.dust_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected dust map ...")

            # Determine the path
            path = fs.join(self.dust_deprojection_path, name + ".fits")

            # Write
            self.dust_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

# -----------------------------------------------------------------
