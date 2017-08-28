#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html.clip Contains the ClipMapsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc
from collections import defaultdict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts
from ....core.tools import browser
from ....core.tools.utils import lazyproperty
from ....core.tools import filesystem as fs
from ....core.tools import sequences
from ....magic.core.frame import Frame
from ..selectioncomponent import MapsSelectionComponent
from ....core.tools.stringify import tostr
from ....core.remote.remote import Remote
from ....core.basics.containers import hashdict
from ....magic.core.mask import Mask
from ....magic.core.alpha import AlphaMask

# -----------------------------------------------------------------

clipped_name = "clipped"
ncolumns = 2
colour_map = "jet"

# -----------------------------------------------------------------

class ClipMapsPageGenerator(MapsSelectionComponent):

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
        super(ClipMapsPageGenerator, self).__init__(*args, **kwargs)

        # Plot paths
        self.old_plot_path = None
        self.young_plot_path = None
        self.ionizing_plot_path = None
        self.dust_plot_path = None

        # Directory paths
        self.old_plot_map_paths = dict()
        self.young_plot_map_paths = dict()
        self.ionizing_plot_map_paths = dict()
        self.dust_plot_map_paths = dict()

        # The plots
        self.old_plots = dict()
        self.young_plots = dict()
        self.ionizing_plots = dict()
        self.dust_plots = dict()

        # The plot paths
        self.old_plot_paths = defaultdict(dict)
        self.young_plot_paths = defaultdict(dict)
        self.ionizing_plot_paths = defaultdict(dict)
        self.dust_plot_paths = defaultdict(dict)

        # The mask plot paths
        self.old_mask_plot_paths = defaultdict(dict)
        self.young_mask_plot_paths = defaultdict(dict)
        self.ionizing_mask_plot_paths = defaultdict(dict)
        self.dust_mask_plot_paths = defaultdict(dict)

        # The sliders
        self.sliders = None

        # The masks
        self.old_masks = dict()
        self.young_masks = dict()
        self.ionizing_masks = dict()
        self.dust_masks = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the maps
        self.load_maps()

        # 3. Set the paths
        self.set_paths()

        # 4. Process the maps
        self.process_maps()

        # 5. Make plots
        self.make_plots()

        # 6. Make sliders
        self.make_sliders()

        # 7. Generate the page
        self.generate_page()

        # 8. Writing
        self.write()

        # 9. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ClipMapsPageGenerator, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

        # Check
        if self.config.default_sigma_level not in self.config.sigma_levels: raise ValueError("Default sigma level must be one of the chosen sigma levels")

        # Make directory to contain the plots
        self.clipped_plots_path = fs.join(self.maps_html_path, clipped_name)
        if fs.is_directory(self.clipped_plots_path):
            if self.config.replot: fs.clear_directory(self.clipped_plots_path)
        else: fs.create_directory(self.clipped_plots_path)

        # Create subdirectories
        if self.config.add_old: self.old_plot_path = fs.create_directory_in(self.clipped_plots_path, "old", clear=self.config.replot_old)
        if self.config.add_young: self.young_plot_path = fs.create_directory_in(self.clipped_plots_path, "young", clear=self.config.replot_young)
        if self.config.add_ionizing: self.ionizing_plot_path = fs.create_directory_in(self.clipped_plots_path, "ionizing", clear=self.config.replot_ionizing)
        if self.config.add_dust: self.dust_plot_path = fs.create_directory_in(self.clipped_plots_path, "dust", clear=self.config.replot_dust)

        # Set random
        if self.config.random: self.config.random_old = self.config.random_young = self.config.random_ionizing = self.config.random_dust = self.config.random

        # Set all
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # Make selections
        self.old_selection = sequences.make_selection(self.old_map_names, self.config.old, self.config.not_old, nrandom=self.config.random_old, all=self.config.all_old, none=not self.config.add_old, indices=self.config.old_indices, not_indices=self.config.not_old_indices)
        self.young_selection = sequences.make_selection(self.young_map_names, self.config.young, self.config.not_young, nrandom=self.config.random_young, all=self.config.all_young, none=not self.config.add_young, indices=self.config.young_indices, not_indices=self.config.not_young_indices)
        self.ionizing_selection = sequences.make_selection(self.ionizing_map_names, self.config.ionizing, self.config.not_ionizing, nrandom=self.config.random_ionizing, all=self.config.all_ionizing, none=not self.config.add_ionizing, indices=self.config.ionizing_indices, not_indices=self.config.not_ionizing_indices)
        self.dust_selection = sequences.make_selection(self.dust_map_names, self.config.dust, self.config.not_dust, nrandom=self.config.random_dust, all=self.config.all_dust, none=not self.config.add_dust, indices=self.config.dust_indices, not_indices=self.config.not_dust_indices)

        # Prompt for user input (if selections not specified)
        self.prompt()

        # Create plot directories for each image
        self.create_plot_directories()

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

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        if self.config.remote is None: return None
        else: return Remote(host_id=self.config.remote)

    # -----------------------------------------------------------------

    def create_plot_directories(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating plot directories for the selected maps ...")

        # Old
        self.create_plot_directories_old()

        # Young
        self.create_plot_directories_young()

        # Ionizing
        self.create_plot_directories_ionizing()

        # Dust
        self.create_plot_directories_dust()

    # -----------------------------------------------------------------

    def create_plot_directories_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stellar maps
        for name in self.old_selection:

            # Determine path for this map
            dirpath = fs.join(self.old_plot_path, name)
            if not fs.is_directory(dirpath): fs.create_directory(dirpath)

            # Set
            self.old_plot_map_paths[name] = dirpath

    # -----------------------------------------------------------------

    def create_plot_directories_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the young stlelar maps
        for name in self.young_selection:

            # Determine the path for this map
            dirpath = fs.join(self.young_plot_path, name)
            if not fs.is_directory(dirpath): fs.create_directory(dirpath)

            # Set
            self.young_plot_map_paths[name] = dirpath

    # -----------------------------------------------------------------

    def create_plot_directories_ionizing(self):

        """
        Thisf unction ...
        :return:
        """

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Determine the path for this map
            dirpath = fs.join(self.ionizing_plot_path, name)
            if not fs.is_directory(dirpath): fs.create_directory(dirpath)

            # Set
            self.ionizing_plot_map_paths[name] = dirpath

    # -----------------------------------------------------------------

    def create_plot_directories_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust maps
        for name in self.dust_selection:

            # Determine the path for this map
            dirpath = fs.join(self.dust_plot_path, name)
            if not fs.is_directory(dirpath): fs.create_directory(dirpath)

            # Set
            self.dust_plot_map_paths[name] = dirpath

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Clip maps"

    # -----------------------------------------------------------------

    @property
    def image_width(self):

        """
        This fucntion ...
        :return:
        """

        return self.config.image_width

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        This function ...
        :return:
        """

        return self.config.image_height

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected maps ...")

        # Load the old stellar maps
        if self.config.add_old: self.load_old()

        # Load the young stellar maps
        if self.config.add_young: self.load_young()

        # Load the ionizing stellar maps
        if self.config.add_ionizing: self.load_ionizing()

        # Load the dust maps
        if self.config.add_dust: self.load_dust()

    # -----------------------------------------------------------------

    def load_old(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the old stellar maps ...")

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
        log.info("Loading the young stellar maps ...")

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
        log.info("Loading the ionizing stellar maps ...")

        # Load
        for name in self.ionizing_selection:

            # Debugging
            log.debug("Loading the '" + name + "' ionizing stellar map ...")

            # Load
            self.ionizing_maps[name] = Frame.from_file(self.ionizing_map_paths[name])

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust maps ...")

        # Load
        for name in self.dust_selection:

            # Debugging
            log.debug("Loading the '" + name + "' dust map ...")

            # Load
            self.dust_maps[name] = Frame.from_file(self.dust_map_paths[name])

    # -----------------------------------------------------------------

    def set_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the plot paths")

        # Old
        self.set_old_paths()

        # Young
        self.set_young_paths()

        # Ionizing
        self.set_ionizing_paths()

        # Dust
        self.set_dust_paths()

    # -----------------------------------------------------------------

    def set_old_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting plot paths for the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Get the origins
            origins = self.old_map_origins[name]

            # Get path
            dirpath = self.old_plot_map_paths[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

            # Make level combinations
            for sigma_levels in sequences.lists_combinations(*sigma_levels):

                # Create dictionary that says which sigma level was used for which frame
                levels_dict = hashdict({name: level for name, level in zip(origins, sigma_levels)})

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + ".png")

                # Set the path
                self.old_plot_paths[name][levels_dict] = filepath

                # Determine the mask filepath
                mask_filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + "_mask.png")

                # Set the path
                self.old_mask_plot_paths[name][levels_dict] = mask_filepath

    # -----------------------------------------------------------------

    def set_young_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting plot paths for the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get the origins
            origins = self.young_map_origins[name]

            # Get path
            dirpath = self.young_plot_map_paths[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

            # Make level combinations
            for sigma_levels in sequences.lists_combinations(*sigma_levels):

                # Create dictionary that says which sigma level was used for which frame
                levels_dict = hashdict({name: level for name, level in zip(origins, sigma_levels)})

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + ".png")

                # Set the path
                self.young_plot_paths[name][levels_dict] = filepath

                # Determine the mask filepath
                mask_filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + "_mask.png")

                # Set the path
                self.young_mask_plot_paths[name][levels_dict] = mask_filepath

    # -----------------------------------------------------------------

    def set_ionizing_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting plot paths for the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Get path
            dirpath = self.ionizing_plot_map_paths[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

            # Make level combinations
            for sigma_levels in sequences.lists_combinations(*sigma_levels):

                # Create dictionary that says which sigma level was used for which frame
                levels_dict = hashdict({name: level for name, level in zip(origins, sigma_levels)})

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + ".png")

                # Set the path
                self.ionizing_plot_paths[name][levels_dict] = filepath

                # Determine the mask filepath
                mask_filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + "_mask.png")

                # Set the path
                self.ionizing_mask_plot_paths[name][levels_dict] = mask_filepath

    # -----------------------------------------------------------------

    def set_dust_paths(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Setting plot paths for the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Get the origins
            origins = self.dust_map_origins[name]

            # Get path
            dirpath = self.dust_plot_map_paths[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

            # Make level combinations
            for sigma_levels in sequences.lists_combinations(*sigma_levels):

                # Create dictionary that says which sigma level was used for which frame
                levels_dict = hashdict({name: level for name, level in zip(origins, sigma_levels)})

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + ".png")

                # Set the path
                self.dust_plot_paths[name][levels_dict] = filepath

                # Determine the mask filepath
                mask_filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + "_mask.png")

                # Set the path
                self.dust_mask_plot_paths[name][levels_dict] = mask_filepath

    # -----------------------------------------------------------------

    def process_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the maps ...")

        # 1. Correct
        self.correct_maps()

        # 2. Crop
        self.crop_maps()

        # 3. Clip
        self.clip_maps()

    # -----------------------------------------------------------------

    def correct_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the maps ...")

        # Old
        if self.config.add_old: self.correct_old_maps()

        # Young
        if self.config.add_young: self.correct_young_maps()

        # Ionizing
        if self.config.add_ionizing: self.correct_ionizing_maps()

        # Dust
        if self.config.add_dust: self.correct_dust_maps()

    # -----------------------------------------------------------------

    def correct_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Correcting the '" + name + "' old stellar map ...")

            # Correct the map
            self.correct_map(self.old_maps[name])

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def correct_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Correcting the '" + name + "' young stellar map ...")

            # Correct the map
            self.correct_map(self.young_maps[name])

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def correct_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Correcting the '" + name + "' ionizing stellar map ...")

            # Correct the map
            self.correct_map(self.ionizing_maps[name])

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def correct_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Correcting the '" + name + "' dust map ...")

            # Correct the map
            self.correct_map(self.dust_maps[name])

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def crop_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the maps ...")

        # Old
        if self.config.add_old: self.crop_old_maps()

        # Young
        if self.config.add_young: self.crop_young_maps()

        # Ionizing
        if self.config.add_ionizing: self.crop_ionizing_maps()

        # Dust
        if self.config.add_dust: self.crop_dust_maps()

    # -----------------------------------------------------------------

    def crop_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Cropping the '" + name + "' old stellar map ...")

            # Crop the map
            self.crop_map(self.old_maps[name], self.config.cropping_factor)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def crop_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Cropping the '" + name + "' young stellar map ...")

            # Crop the map
            self.crop_map(self.young_maps[name], self.config.cropping_factor)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def crop_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Cropping the '" + name + "' ionizing stellar map ...")

            # Crop the map
            self.crop_map(self.ionizing_maps[name], self.config.cropping_factor)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def crop_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Cropping the '" + name + "' dust map ...")

            # Crop the map
            self.crop_map(self.dust_maps[name], self.config.cropping_factor)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def clip_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the maps ...")

        # Old
        if self.config.add_old: self.clip_old_maps()

        # Young
        if self.config.add_young: self.clip_young_maps()

        # Ionizing
        if self.config.add_ionizing: self.clip_ionizing_maps()

        # Dust
        if self.config.add_dust: self.clip_dust_maps()

    # -----------------------------------------------------------------

    def present_old_plots_level_combinations_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        combinations = []

        # Loop over the level combinations
        for levels in self.old_plot_paths[name]:

            # Get the filepath
            filepath = self.old_plot_paths[name][levels]

            # Check
            if fs.is_file(filepath):
                log.success("The plot of the '" + name + "' old stellar map for significance levels [" + self.levels_to_string(levels) + "'] is already present")
                combinations.append(levels)

        # Return the level combinations
        return combinations

    # -----------------------------------------------------------------

    def has_all_old_plots_for(self, name):

        """
        This function ...
        :param name:
        :param origins:
        :param sigma_levels:
        :return:
        """

        # Loop over the level combinations
        for levels in self.old_plot_paths[name]:

            # Get the filepath
            filepath = self.old_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def present_young_plots_level_combinations_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        combinations = []

        # Loop over the level combinations
        for levels in self.young_plot_paths[name]:

            # Get the filepath
            filepath = self.young_plot_paths[name][levels]

            # Check
            if fs.is_file(filepath):
                log.success("The plot of the '" + name + "' young stellar map for significance levels [" + self.levels_to_string(levels) + "'] is already present")
                combinations.append(levels)

        # Return the level combinations
        return combinations

    # -----------------------------------------------------------------

    def has_all_young_plots_for(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.young_plot_paths[name]:

            # Get the filepath
            filepath = self.young_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def present_ionizing_plots_level_combinations_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        combinations = []

        # Loop over the level combinations
        for levels in self.ionizing_plot_paths[name]:

            # Get the filepath
            filepath = self.ionizing_plot_paths[name][levels]

            # Check
            if fs.is_file(filepath):
                log.success("The plot of the '" + name + "' ionizing stellar map for significance levels [" + self.levels_to_string(levels) + "'] is already present")
                combinations.append(levels)

        # Return the level combinations
        return combinations

    # -----------------------------------------------------------------

    def has_all_ionizing_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.ionizing_plot_paths[name]:

            # Get the filepath
            filepath = self.ionizing_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def present_dust_plots_level_combinations_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        combinations = []

        # Loop over the level combinations
        for levels in self.dust_plot_paths[name]:

            # Get the filepath
            filepath = self.dust_plot_paths[name][levels]

            # Check
            if fs.is_file(filepath):
                log.success("The plot of the '" + name + "' dust map for significance levels [" + self.levels_to_string(levels) + "'] is already present")
                combinations.append(levels)

        # Return the level combinations
        return combinations

    # -----------------------------------------------------------------

    def has_all_dust_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.dust_plot_paths[name]:

            # Get the filepath
            filepath = self.dust_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): return False

        # All checks passed
        return True

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

            # Check
            if self.has_all_old_plots_for(name):
                log.success("All plots for the '" + name + "' image at the requested sigma levels are already present")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' old stellar map ...")

            # Get the origins
            origins = self.old_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.old_maps[name], origins, sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity,
                                          present=self.present_old_plots_level_combinations_for(name),
                                          fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                          fuzziness_offset=self.config.fuzzy_min_significance_offset, return_masks=True)

            # Replace by a dictionary of maps
            self.old_maps[name] = maps

            # Set the masks
            self.old_masks[name] = masks

            # Cleanup
            gc.collect()

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

            # Check
            if self.has_all_young_plots_for(name):
                log.success("All plots for the '" + name + "' image at requested sigma levels are already present")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            origins = self.young_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.young_maps[name], origins, sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity,
                                          present=self.present_young_plots_level_combinations_for(name),
                                          fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                          fuzziness_offset=self.config.fuzzy_min_significance_offset, return_masks=True)

            # Replace by a dictionary of maps
            self.young_maps[name] = maps

            # Set the masks
            self.young_masks[name] = masks

            # Cleanup
            gc.collect()

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

            # Check
            if self.has_all_ionizing_plots_for(name):
                log.success("All plots for the '" + name + "' image at requested sigma levels are already present")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.ionizing_maps[name], origins, sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity,
                                          present=self.present_ionizing_plots_level_combinations_for(name),
                                          fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                          fuzziness_offset=self.config.fuzzy_min_significance_offset, return_masks=True)

            # Replace by a dictionary of maps
            self.ionizing_maps[name] = maps

            # Set the masks
            self.ionizing_masks[name] = masks

            # Cleanup
            gc.collect()

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

            # Check
            if self.has_all_dust_plots_for(name):
                log.success("All plots for the '" + name + "' image at requested sigma levels are already present")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            origins = self.dust_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.dust_maps[name], origins, sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity,
                                          present=self.present_dust_plots_level_combinations_for(name),
                                          fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                          fuzziness_offset=self.config.fuzzy_min_significance_offset, return_masks=True)

            # Replace by a dictionary of maps
            self.dust_maps[name] = maps

            # Set the masks
            self.dust_masks[name] = masks

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    @lazyproperty
    def old_maps_image_names(self):

        """
        This function ...
        :return:
        """

        # names = set()
        # for map_name in self.old_maps:
        #     names.update(self.old_maps[map_name].keys())
        # return list(names)

        # WORKS BEFORE MAPS HAVE BEEN ADDED
        names = set()
        for map_name in self.old_selection:
            origins = self.old_map_origins[map_name]
            names.update(origins)
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_maps_image_names(self):

        """
        This function ...
        :return:
        """

        # names = set()
        # for map_name in self.young_maps:
        #     names.update(self.young_maps[map_name].keys())
        # return list(names)

        # WORKS BEFORE MAPS HAVE BEEN ADDED
        names = set()
        for map_name in self.young_selection:
            origins = self.young_map_origins[map_name]
            names.update(origins)
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_maps_image_names(self):

        """
        This function ...
        :return:
        """

        # names = set()
        # for map_name in self.ionizing_maps:
        #     names.update(self.ionizing_maps[map_name].keys())
        # return list(names)

        # WORKS BEFORE MAPS HAVE BEEN ADDED
        names = set()
        for map_name in self.ionizing_selection:
            origins = self.ionizing_map_origins[map_name]
            names.update(origins)
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_maps_image_names(self):

        """
        This function ...
        :return:
        """

        # names = set()
        # for map_name in self.dust_maps:
        #     names.update(self.dust_maps[map_name].keys())
        # return list(names)

        # WORKS BEFORE MAPS HAVE BEEN ADDED
        names = set()
        for map_name in self.dust_selection:
            origins = self.dust_map_origins[map_name]
            names.update(origins)
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        if self.config.add_old: names.update(self.old_maps_image_names)
        if self.config.add_young: names.update(self.young_maps_image_names)
        if self.config.add_ionizing: names.update(self.ionizing_maps_image_names)
        if self.config.add_dust: names.update(self.dust_maps_image_names)
        return names

    # -----------------------------------------------------------------

    def levels_to_string(self, levels):

        """
        This function ...
        :param levels:
        :return:
        """

        return tostr(levels, identity_symbol="_", delimiter="__", replace_spaces_keys="-", quote_key=False, quote_value=False)

    # -----------------------------------------------------------------

    def make_rgba_plot(self, the_map, filepath):

        """
        This function ...
        :param the_map:
        :param filepath:
        :return:
        """

        # Plot as RGB
        the_map.saveto_png(filepath, interval=self.config.interval, scale=self.config.scale,
                           alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha,
                           colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_mask_plot(self, mask, filepath):

        """
        This function ...
        :param mask:
        :param filepath:
        :return:
        """

        if isinstance(mask, AlphaMask): mask.saveto_png(filepath, colour="black", alpha=True) # alpha to recognize alpha masks
        elif isinstance(mask, Mask): mask.saveto_png(filepath, colour="black", alpha=False) # no alpha to recognize regular masks
        else: raise ValueError("Not a valid mask")

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Old
        if self.config.add_old: self.make_old_plots()

        # Young
        if self.config.add_young: self.make_young_plots()

        # Ionizing
        if self.config.add_ionizing: self.make_ionizing_plots()

        # Dust
        if self.config.add_dust: self.make_dust_plots()

    # -----------------------------------------------------------------

    def make_old_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the old stellar maps ...")

        # Maps
        self.make_old_map_plots()

        # Masks
        self.make_old_mask_plots()

    # -----------------------------------------------------------------

    def make_old_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' old stellar map ...")

            # Loop over the levels dicts
            for levels in self.old_maps[name]:

                # Determine the filepath
                filepath = self.old_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.old_maps[name][levels], filepath)

            # Clear
            del self.old_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_old_mask_plots(self):

        """
        This function ...
        :return: 
        """

        # Loop over the maps
        for name in self.old_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' old stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.old_maps[name]:

                # Determine the filepath
                filepath = self.old_mask_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as PNG file
                self.make_mask_plot(self.old_masks[name][levels], filepath)

            # Clear
            del self.old_masks[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_young_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the young stellar maps ...")

        # Maps
        self.make_young_map_plots()

        # Masks
        self.make_young_mask_plots()

    # -----------------------------------------------------------------

    def make_young_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' young stellar map ...")

            # Loop over the levels dicts
            for levels in self.young_maps[name]:

                # Determine the filepath
                filepath = self.young_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.young_maps[name][levels], filepath)

            # Clear
            del self.young_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_young_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the masp
        for name in self.young_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' young stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.young_maps[name]:

                # Determine the filepath
                filepath = self.young_mask_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as PNG
                self.make_mask_plot(self.young_masks[name][levels], filepath)

            # Clear
            del self.young_masks[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_ionizing_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the ionizing stellar maps ...")

        # Maps
        self.make_ionizing_map_plots()

        # Masks
        self.make_ionizing_mask_plots()

    # -----------------------------------------------------------------

    def make_ionizing_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' ionizing stellar map ...")

            # Loop over the levels dicts
            for levels in self.ionizing_maps[name]:

                # Determine the filepath
                filepath = self.ionizing_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.ionizing_maps[name][levels], filepath)

            # Clear
            del self.ionizing_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_ionizing_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' ionizing stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.ionizing_maps[name]:

                # Determine the filepath
                filepath = self.ionizing_mask_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as PNG
                self.make_mask_plot(self.ionizing_masks[name][levels], filepath)

            # Clear
            del self.ionizing_masks[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_dust_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the dust maps ...")

        # Maps
        self.make_dust_map_plots()

        # Masks
        self.make_dust_mask_plots()

    # -----------------------------------------------------------------

    def make_dust_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' dust map ...")

            # Loop over the levels dicts
            for levels in self.dust_maps[name]:

                # Determine the filepath
                filepath = self.dust_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.dust_maps[name][levels], filepath)

            # Clear
            del self.dust_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def make_dust_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' dust map masks ...")

            # Loop over the levels dicts
            for levels in self.dust_maps[name]:

                # Determine the filepath
                filepath = self.dust_mask_plot_paths[name][levels]

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as PNG
                self.make_mask_plot(self.dust_masks[name][levels], filepath)

            # Clear
            del self.dust_masks[name]
            gc.collect()

    # -----------------------------------------------------------------

    def get_sigma_levels_for_origins(self, origins, as_list=False):

        """
        This function ...
        :param origins:
        :param as_list:
        :return:
        """

        if as_list:

            levels = []
            #print(self.sigma_levels_for_images)
            for origin in origins: levels.append(self.sigma_levels_for_images[origin])
            return levels

        else:

            levels = dict()
            for origin in origins: levels[origin] = self.sigma_levels_for_images[origin]
            return levels

    # -----------------------------------------------------------------

    @lazyproperty
    def sigma_levels_for_images(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        levels = dict()

        # Loop over the images
        for name in self.image_names:

            # Get the level
            level = self.significance_levels[name]

            # Get the different sigma levels
            sigma_levels = [level * factor for factor in self.config.sigma_levels]

            # Set the levels
            levels[name] = sigma_levels

        # Return the levels
        return levels

    # -----------------------------------------------------------------

    @lazyproperty
    def default_sigma_levels_for_images(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        levels = dict()

        # Loop over the images
        for name in self.image_names:

            # Get the level
            level = self.significance_levels[name]

            # Set the default level
            default_level = self.config.default_sigma_level * level

            # Set the level
            levels[name] = default_level

        # Return the dictionary of default sigma levels
        return levels

    # -----------------------------------------------------------------

    def make_sliders(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sliders ...")

        # Group map names
        names = dict()
        names["old"] = self.old_selection
        names["young"] = self.young_selection
        names["ionizing"] = self.ionizing_selection
        names["dust"] = self.dust_selection

        # Group paths
        paths = dict()
        paths["old"] = self.old_plot_paths
        paths["young"] = self.young_plot_paths
        paths["ionizing"] = self.ionizing_plot_paths
        paths["dust"] = self.dust_plot_paths

        # Make the slider
        self.sliders = html.make_multi_image_sliders(names, paths, self.image_names, self.sigma_levels_for_images, self.default_sigma_levels_for_images, width=None, height=None, basic=True, img_class=None)

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = HTMLPage(self.title, style=page_style, css_path=css_paths, javascript_path=javascripts, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes))

        self.page += html.newline

        # Add the sliders
        self.page += self.sliders

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the page
        self.write_page()

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Save
        self.page.saveto(self.clip_maps_html_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.clip_maps_html_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return None

# -----------------------------------------------------------------
