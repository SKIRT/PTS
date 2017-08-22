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

        # The sliders
        self.sliders = None

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

        # 3. Process the maps
        self.process_maps()

        # 4. Make plots
        self.make_plots()

        # 5. Make sliders
        self.make_sliders()

        # 5. Generate the page
        self.generate_page()

        # 6. Writing
        self.write()

        # Show
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

        # Make directory to contain the plots
        self.clipped_plots_path = fs.join(self.maps_html_path, clipped_name)
        if fs.is_directory(self.clipped_plots_path):
            if self.config.replot: fs.clear_directory(self.clipped_plots_path)
        else: fs.create_directory(self.clipped_plots_path)

        # Create subdirectories
        if self.config.add_old: self.old_plot_path = fs.create_directory_in(self.clipped_plots_path, "old")
        if self.config.add_young: self.young_plot_path = fs.create_directory_in(self.clipped_plots_path, "young")
        if self.config.add_ionizing: self.ionizing_plot_path = fs.create_directory_in(self.clipped_plots_path, "ionizing")
        if self.config.add_dust: self.dust_plot_path = fs.create_directory_in(self.clipped_plots_path, "dust")

        # Set random
        if self.config.random: self.config.random_old = self.config.random_young = self.config.random_ionizing = self.config.random_dust = self.config.random

        # Set all
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # Make selections
        self.old_selection = sequences.make_selection(self.old_map_names, self.config.old, self.config.not_old, nrandom=self.config.random_old, all=self.config.all_old, none=not self.config.add_old)
        self.young_selection = sequences.make_selection(self.young_map_names, self.config.young, self.config.not_young, nrandom=self.config.random_young, all=self.config.all_young, none=not self.config.add_young)
        self.ionizing_selection = sequences.make_selection(self.ionizing_map_names, self.config.ionizing, self.config.not_ionizing, nrandom=self.config.random_ionizing, all=self.config.all_ionizing, none=not self.config.add_ionizing)
        self.dust_selection = sequences.make_selection(self.dust_map_names, self.config.dust, self.config.not_dust, nrandom=self.config.random_dust, all=self.config.all_dust, none=not self.config.add_dust)

        # Create plot directories for each image
        self.create_plot_directories()

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

    # @lazyproperty
    # def maps_filters(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return self.get_all_filters()

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

            # Clip the map
            maps = self.make_clipped_maps(self.old_maps[name], origins, self.config.sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # Replace by a dictionary of maps
            self.old_maps[name] = maps

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

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            origins = self.young_map_origins[name]

            # Clip the map
            maps = self.make_clipped_maps(self.young_maps[name], origins, self.config.sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # Replace by a dictionary of maps
            self.young_maps[name] = maps

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

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Clip the map
            maps = self.make_clipped_maps(self.ionizing_maps[name], origins, self.config.sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # Replace by a dictionary of maps
            self.ionizing_maps[name] = maps

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

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            origins = self.dust_map_origins[name]

            # Clip the map
            maps = self.make_clipped_maps(self.dust_maps[name], origins, self.config.sigma_levels,
                                          convolve=self.config.convolve, remote=self.remote,
                                          rebin_remote_threshold=self.config.rebin_remote_threshold,
                                          npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # Replace by a dictionary of maps
            self.dust_maps[name] = maps

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    @lazyproperty
    def old_maps_image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        for map_name in self.old_maps:
            names.update(self.old_maps[map_name].keys())
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_maps_image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        for map_name in self.young_maps:
            names.update(self.young_maps[map_name].keys())
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_maps_image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        for map_name in self.ionizing_maps:
            names.update(self.ionizing_maps[map_name].keys())
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_maps_image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        for map_name in self.dust_maps:
            names.update(self.dust_maps[map_name].keys())
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

        # Loop over the maps
        for name in self.old_maps:

            # Get path
            dirpath = self.old_plot_map_paths[name]

            # Debugging
            log.debug("Making plots of the '" + name + "' old stellar map ...")

            # Loop over the levels dicts
            for levels in self.old_maps[name]:

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels) + ".png")

                # Set the path
                self.old_plot_paths[name][levels] = filepath

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.old_maps[name][levels], filepath)

    # -----------------------------------------------------------------

    def make_young_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get path
            dirpath = self.young_plot_map_paths[name]

            # Debugging
            log.debug("Making plots of the '" + name + "' young stellar map ...")

            # Loop over the levels dicts
            for levels in self.young_maps[name]:

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels) + ".png")

                # Set the path
                self.young_plot_paths[name][levels] = filepath

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.young_maps[name][levels], filepath)

    # -----------------------------------------------------------------

    def make_ionizing_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Determine the path
            dirpath = self.ionizing_plot_map_paths[name]

            # Debugging
            log.debug("Making plots of the '" + name + "' ionizing stellar map ...")

            # Loop over the levels dicts
            for levels in self.ionizing_maps[name]:

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels) + ".png")

                # Set the path
                self.ionizing_plot_paths[name][levels] = filepath

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.ionizing_maps[name][levels], filepath)

    # -----------------------------------------------------------------

    def make_dust_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Determine the path
            dirpath = self.dust_plot_map_paths[name]

            # Debugging
            log.debug("Making plots of the '" + name + "' dust map ...")

            # Loop over the levels dicts
            for levels in self.dust_maps[name]:

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels) + ".png")

                # Set the path
                self.dust_plot_paths[name][levels] = filepath

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.dust_maps[name][levels], filepath)

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
        self.sliders = html.make_multi_image_sliders(names, paths, self.image_names, self.config.sigma_levels, self.config.default_sigma_level, width=None, height=None, basic=True, img_class=None)

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
