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
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts
from ....core.tools import browser
from ....core.tools.utils import lazyproperty
from ....core.tools import filesystem as fs
from ....core.tools import sequences
from ....magic.core.frame import Frame
from ..selectioncomponent import MapsSelectionComponent
from ....core.tools.stringify import tostr, stringify_list_fancy
from ....core.remote.remote import Remote
from ....core.basics.containers import hashdict
from ....magic.core.mask import Mask
from ....magic.core.alpha import AlphaMask, load_mask_or_alpha_mask
from ....core.tools.utils import memoize_method
from ....core.basics import containers
from ....core.tools.serialization import load_dict
from ....core.tools.stringify import stringify_dict

# -----------------------------------------------------------------

clipped_name = "clipped"
clipped_data_name = "clipped_data"
image_masks_name = "image_masks"
ncolumns = 2

# -----------------------------------------------------------------

old_name = "old"
young_name = "young"
ionizing_name = "ionizing"
dust_name = "dust"

# -----------------------------------------------------------------

page_width = 1000

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

        # Paths
        self.clipped_plots_path = None
        self.image_mask_plots_path = None

        # Plot paths for image masks
        self.image_mask_plot_paths = dict()
        self.image_mask_plot_paths_images = defaultdict(dict)

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

        # The clipped maps
        self.old_clipped_maps = dict()
        self.young_clipped_maps = dict()
        self.ionizing_clipped_maps = dict()
        self.dust_clipped_maps = dict()

        # The masks
        self.old_masks = dict()
        self.young_masks = dict()
        self.ionizing_masks = dict()
        self.dust_masks = dict()

        # Has all data (clipped maps and masks)?
        self.has_all_data_old = dict()
        self.has_all_data_young = dict()
        self.has_all_data_ionizing = dict()
        self.has_all_data_dust = dict()

        # The page
        self.page = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Set the paths
        self.set_paths()

        # 3. Load the maps
        if not self.has_all_plots: self.load_maps()

        # 4. Check present maps
        if not self.config.replot: self.check_present()

        # Check data
        if not self.config.clear_data: self.check_data()

        # 5. Process the maps
        if not self.has_all_plots: self.process_maps()

        # 6. Make plots
        if not self.has_all_plots: self.make_plots()

        # 7. Make image mask plots
        if not self.has_all_plots: self.plot_image_masks()

        # 8. Make sliders
        self.make_sliders()

        # 9. Generate the page
        self.generate_page()

        # 10. Writing
        self.write()

        # 11. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    @lazyproperty
    def do_old(self):

        """
        This function ...
        :return:
        """

        if self.config.only is not None: return "old" in self.config.only
        else: return self.config.add_old

    # -----------------------------------------------------------------

    @lazyproperty
    def do_young(self):

        """
        This function ...
        :return:
        """

        if self.config.only is not None: return "young" in self.config.only
        else: return self.config.add_young

    # -----------------------------------------------------------------

    @lazyproperty
    def do_ionizing(self):

        """
        This function ...
        :return:
        """

        if self.config.only is not None: return "ionizing" in self.config.only
        else: return self.config.add_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def do_dust(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.only is not None: return "dust" in self.config.only
        else: return self.config.add_dust

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
        if self.config.add_old: self.old_plot_path = fs.create_directory_in(self.clipped_plots_path, old_name, clear=self.config.replot_old)
        if self.config.add_young: self.young_plot_path = fs.create_directory_in(self.clipped_plots_path, young_name, clear=self.config.replot_young)
        if self.config.add_ionizing: self.ionizing_plot_path = fs.create_directory_in(self.clipped_plots_path, ionizing_name, clear=self.config.replot_ionizing)
        if self.config.add_dust: self.dust_plot_path = fs.create_directory_in(self.clipped_plots_path, dust_name, clear=self.config.replot_dust)

        # Create directory to contain the image mask plots
        self.image_mask_plots_path = fs.join(self.maps_html_path, image_masks_name)
        if fs.is_directory(self.image_mask_plots_path):
            if self.config.replot or self.config.replot_image_masks: fs.clear_directory(self.image_mask_plots_path)
        else: fs.create_directory(self.image_mask_plots_path)

        # Clear data?
        if self.config.clear_data and fs.is_directory(self.clipped_data_path): fs.clear_directory(self.clipped_data_path)
        if self.config.clear_old_data and fs.is_directory(self.clipped_data_old_path): fs.remove_directory(self.clipped_data_old_path)
        if self.config.clear_young_data and fs.is_directory(self.clipped_data_young_path): fs.remove_directory(self.clipped_data_young_path)
        if self.config.clear_ionizing_data and fs.is_directory(self.clipped_data_ionizing_path): fs.remove_directory(self.clipped_data_ionizing_path)
        if self.config.clear_dust_data and fs.is_directory(self.clipped_data_dust_path): fs.remove_directory(self.clipped_data_dust_path)

        # Copy data?
        if self.config.data_from is not None:
            path = fs.absolute_or_in(self.config.data_from, self.maps_html_path)
            old_path = fs.join(path, "old")
            young_path = fs.join(path, "young")
            ionizing_path = fs.join(path, "ionizing")
            dust_path = fs.join(path, "dust")
            if self.do_old: fs.clear_and_copy_from_directory(old_path, self.clipped_data_old_path)
            if self.do_young: fs.clear_and_copy_from_directory(young_path, self.clipped_data_young_path)
            if self.do_ionizing: fs.clear_and_copy_from_directory(ionizing_path, self.clipped_data_ionizing_path)
            if self.do_dust: fs.clear_and_copy_from_directory(dust_path, self.clipped_data_dust_path)

        # Load the selection
        self.load_selection()

        # Create plot directories for each image
        self.create_plot_directories()

        # Create plot directories for image masks
        self.create_image_mask_plot_directories()

    # -----------------------------------------------------------------

    @lazyproperty
    def clipped_data_path(self):

        """
        Thisf unction ...
        :return:
        """

        path = fs.join(self.maps_html_path, clipped_data_name)
        if self.config.write_data and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def clipped_data_old_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.clipped_data_path, old_name)
        if self.config.write_data and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def clipped_data_young_path(self):

        """
        Thisn function ...
        :return:
        """

        path = fs.join(self.clipped_data_path, young_name)
        if self.config.write_data and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def clipped_data_ionizing_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.clipped_data_path, ionizing_name)
        if self.config.write_data and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def clipped_data_dust_path(self):

        """
        Thisn function ...
        :return:
        """

        path = fs.join(self.clipped_data_path, dust_name)
        if self.config.write_data and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def clipped_data_old_path_for_map(self, name, create=True):

        """
        Thisf unction ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.clipped_data_old_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def clipped_data_young_path_for_map(self, name, create=True):

        """
        Thisn function ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.clipped_data_young_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def clipped_data_ionizing_path_for_map(self, name, create=True):

        """
        Thisn function ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.clipped_data_ionizing_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def clipped_data_dust_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        path = fs.join(self.clipped_data_dust_path, name)
        if create and not fs.is_directory(path): fs.create_directory(path)
        return path

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

    def load_selection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.debug("Loading maps selections ...")

        # Determine the path
        path = fs.join(self.maps_components_path, "selection_" + str(self.config.selection) + ".dat")

        # Load the selection file
        selection = load_dict(path)

        # Set the selections
        self.old_selection = selection["old"]
        self.young_selection = selection["young"]
        self.ionizing_selection = selection["ionizing"]
        self.dust_selection = selection["dust"]

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

    def create_image_mask_plot_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the plot directories for the image masks ...")

        # Loop over the images
        for name in self.image_names:

            # Determine the path
            path = fs.join(self.image_mask_plots_path, name)

            # Set
            self.image_mask_plot_paths[name] = path

            # Create the directory if necessary
            if not fs.is_directory(path): fs.create_directory(path)

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

    def set_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the plot paths ...")

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
        for name in self.old_selection:

            # Get path
            dirpath = self.old_plot_map_paths[name]

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_old_map(name):

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
        for name in self.young_selection:

            # Get path
            dirpath = self.young_plot_map_paths[name]

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_young_map(name):

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
        for name in self.ionizing_selection:

            # Get path
            dirpath = self.ionizing_plot_map_paths[name]

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_ionizing_map(name):

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
        for name in self.dust_selection:

            # Get path
            dirpath = self.dust_plot_map_paths[name]

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_dust_map(name):

                # Determine the filepath
                filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + ".png")

                # Set the path
                self.dust_plot_paths[name][levels_dict] = filepath

                # Determine the mask filepath
                mask_filepath = fs.join(dirpath, self.levels_to_string(levels_dict) + "_mask.png")

                # Set the path
                self.dust_mask_plot_paths[name][levels_dict] = mask_filepath

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected maps ...")

        # Load the old stellar maps
        if self.do_old: self.load_old()

        # Load the young stellar maps
        if self.do_young: self.load_young()

        # Load the ionizing stellar maps
        if self.do_ionizing: self.load_ionizing()

        # Load the dust maps
        if self.do_dust: self.load_dust()

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

    def check_present(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking which maps are already present ...")

        # Old stellar maps
        if not self.config.replot_old: self.check_old()

        # Young stellar maps
        if not self.config.replot_young: self.check_young()

        # Ionizing stellar maps
        if not self.config.replot_ionizing: self.check_ionizing()

        # Dust maps
        if not self.config.replot_dust: self.check_dust()

    # -----------------------------------------------------------------

    def check_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the old stellar maps ...")

        # Loop over the maps
        for name in self.old_selection:

            # Check
            if self.has_all_old_plots_for(name):
                log.success("All plots for the '" + name + "' old stellar map at the requested sigma levels are already present")
                if name in self.old_maps: del self.old_maps[name]

    # -----------------------------------------------------------------

    def check_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the young stellar maps ...")

        # Loop over the maps
        for name in self.young_selection:

            # Check
            if self.has_all_young_plots_for(name):
                log.success("All plots for the '" + name + "' young stellar map at the requested sigma levels are already present")
                if name in self.young_maps: del self.young_maps[name]

    # -----------------------------------------------------------------

    def check_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_selection:

            # Check
            if self.has_all_ionizing_plots_for(name):
                log.success("All plots for the '" + name + "' ionizing stellar map at the requested sigma levels are already present")
                if name in self.ionizing_maps: del self.ionizing_maps[name]

    # -----------------------------------------------------------------

    def check_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the dust maps ...")

        # Loop over the maps
        for name in self.dust_selection:

            # Check
            if self.has_all_dust_plots_for(name):
                log.success("All plots for the '" + name + "' dust map at the requested sigma levels are already present")
                if name in self.dust_maps: del self.dust_maps[name]

    # -----------------------------------------------------------------

    @lazyproperty
    def process_old_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_maps.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def process_young_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_maps.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def process_ionizing_map_names(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_maps.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def process_dust_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps.keys()

    # -----------------------------------------------------------------

    def get_map_and_masks_filenames_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Create string for the levels
        #levels_string = stringify_dict(levels_dict, quote_key=False, quote_value=False, identity_symbol="_", delimiter="__", value_delimiter='-')[1] # key_delimiter is not yet implemented, value_delimiter is not right
        levels_string = stringify_dict(levels_dict, quote_key=False, quote_value=False, identity_symbol="_", delimiter="__", replace_spaces_keys="-")[1]

        # Determine path for the map and mask
        map_filename = name + "___" + levels_string + ".fits"
        mask_filename = name + "___" + levels_string + "_mask.fits"

        # Return the names
        return map_filename, mask_filename

    # -----------------------------------------------------------------

    def check_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking existing data ...")

        # Old
        if self.do_old: self.check_data_old()

        # Young
        if self.do_young: self.check_data_young()

        # Ionizing
        if self.do_ionizing: self.check_data_ionizing()

        # Dust
        if self.do_dust: self.check_data_dust()

    # -----------------------------------------------------------------

    def get_sigma_levels_combinations_for_old_map(self, name):

        """
        Thisnf unction ...
        :param name:
        :return:
        """

        # Get the origins
        #origins = self.old_map_origins[name]
        origins = self.filtered_old_map_origins[name]

        # Get image names
        names = [str(fltr) for fltr in origins]

        # Get the corresponding sigma levels
        sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

        # Make the combinations
        return names, sequences.lists_combinations(*sigma_levels)

    # -----------------------------------------------------------------

    def get_levels_dicts_for_old_map(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Geneerate the combinations
        names, sigma_levels_combinations = self.get_sigma_levels_combinations_for_old_map(name)

        # Loop over the combinations
        for combination in sigma_levels_combinations:

            # Create dictionary that says which sigma level was used for which frame
            levels_dict = hashdict({name: level for name, level in zip(names, combination)})

            # Yield
            yield levels_dict

    # -----------------------------------------------------------------

    def check_data_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking existing data for old stellar maps ...")

        # Loop over the maps
        for map_name in self.process_old_map_names:

            # ALl plots already
            if self.has_all_old_map_and_mask_plots(map_name): continue

            # Debugging
            log.debug("Checking data for the '" + map_name + "' old stellar map ...")

            # Has all?
            has_all = True
            has_any = False

            # Initialize list of required filepaths
            required_filepaths = []

            # Loop over the different levels dictionaries
            for levels_dict in self.get_levels_dicts_for_old_map(map_name):

                # Already plotted
                if self.has_old_map_plot_for_levels(map_name, levels_dict) and self.has_old_mask_plot_for_levels(map_name, levels_dict): continue

                # Debugging
                levels_string = self.levels_to_string(levels_dict)
                log.debug("Checking data for the combination '" + levels_string + "' ...")

                # Debugging
                map_filepath, mask_filepath = self.get_old_map_and_mask_paths_for_levels(map_name, levels_dict, create=False)
                log.debug("Looking for files: '" + map_filepath + "' and '" + mask_filepath + "' ...")

                # Add filepaths
                required_filepaths.append(map_filepath)
                required_filepaths.append(mask_filepath)

                # Check whether file is present
                if self.has_old_map_and_mask_for_levels(map_name, levels_dict):

                    has_any = True

                    # Success
                    log.success("Data for combination '" + levels_string + "' is present")

                    # Load map and mask
                    the_map, mask = self.load_old_map_and_mask_for_levels(map_name, levels_dict)

                    # Set flag
                    the_map.metadata["finished"] = False

                    # Initialize
                    if map_name not in self.old_clipped_maps: self.old_clipped_maps[map_name] = dict()
                    if map_name not in self.old_masks: self.old_masks[map_name] = dict()

                    # Replace by a dictionary of maps
                    self.old_clipped_maps[map_name][levels_dict] = the_map

                    # Set the masks
                    self.old_masks[map_name][levels_dict] = mask

                # Doesn't have map and mask for each levels dict
                else: has_all = False

            # Remove other?
            if self.config.remove_other_data:
                log.debug("Removing the following files:")
                log.debug("")
                path = self.clipped_data_old_path_for_map(map_name, create=True)
                filenames = [fs.strip_extension(fs.name(filepath)) for filepath in required_filepaths]
                other_filepaths, other_filenames = fs.files_in_path(path, exact_not_name=filenames, returns=["path", "name"], unpack=True)
                for line in stringify_list_fancy(other_filenames, lines_prefix="  ")[1].split("\n"): log.debug(line)
                fs.remove_files(other_filepaths)
                log.debug("")

            # Set has all flag
            self.has_all_data_old[map_name] = has_all

            # If has all: remove the map
            if has_all and not self.config.reclip_from_masks: del self.old_maps[map_name]

            # Succes
            if has_all: log.success("All data is already present for the '" + map_name + "' old stellar map")

            # For reclipping
            if self.config.reclip_from_masks and has_any:

                # Message
                log.success("But reclipping will be applied")

                # Debugging
                log.debug("Correcting the '" + map_name + "' old stellar map to be reclipped ...")

                # Correct
                self.correct_map(self.old_maps[map_name])

                # Debugging
                log.debug("Rebinning the '" + map_name + "' old stellar map to be reclipped ...")

                # Rebin (so to be cropped as the clipped maps)
                last_levels_dict = self.old_clipped_maps[map_name].keys()[-1]
                self.old_maps[map_name].rebin(self.old_clipped_maps[map_name][last_levels_dict].wcs)

    # -----------------------------------------------------------------

    def get_sigma_levels_combinations_for_young_map(self, name):

        """
        Thisnf unction ...
        :param name:
        :return:
        """

        # Get the origins
        #origins = self.young_map_origins[name]
        origins = self.filtered_young_map_origins[name]

        # Get image names
        names = [str(fltr) for fltr in origins]

        # Get the corresponding sigma levels
        sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

        # Make the combinations
        return names, sequences.lists_combinations(*sigma_levels)

    # -----------------------------------------------------------------

    def get_levels_dicts_for_young_map(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Loop over the different combinations
        names, sigma_levels_combinations = self.get_sigma_levels_combinations_for_young_map(name)

        # Loop over the different combinations
        for combination in sigma_levels_combinations:

            # Create dictionary that says which sigma level was used for which frame
            levels_dict = hashdict({name: level for name, level in zip(names, combination)})

            # Yield
            yield levels_dict

    # -----------------------------------------------------------------

    def check_data_young(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Checking existing data for young stellar maps ...")

        # Loop over the maps
        for map_name in self.process_young_map_names:

            # ALl plots already
            if self.has_all_young_map_and_mask_plots(map_name): continue

            # Debugging
            log.debug("Checking data for the '" + map_name + "' young stellar map ...")

            # Has all?
            has_all = True
            has_any = False

            # Initialize list of required filepaths
            required_filepaths = []

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_young_map(map_name):

                # Already plotted
                if self.has_young_map_plot_for_levels(map_name, levels_dict) and self.has_young_mask_plot_for_levels(map_name, levels_dict): continue

                # Debugging
                levels_string = self.levels_to_string(levels_dict)
                log.debug("Checking data for the combination '" + levels_string + "' ...")

                # Debugging
                map_filepath, mask_filepath = self.get_young_map_and_mask_paths_for_levels(map_name, levels_dict, create=False)
                log.debug("Looking for files: '" + map_filepath + "' and '" + mask_filepath + "' ...")

                # Add filepaths
                required_filepaths.append(map_filepath)
                required_filepaths.append(mask_filepath)

                # Check whether file is present
                if self.has_young_map_and_mask_for_levels(map_name, levels_dict):

                    has_any = True

                    # Success
                    log.success("Data for combination '" + levels_string + "' is present")

                    # Load map and mask
                    the_map, mask = self.load_young_map_and_mask_for_levels(map_name, levels_dict)

                    # Set flag
                    the_map.metadata["finished"] = False

                    # Initialize
                    if map_name not in self.young_clipped_maps: self.young_clipped_maps[map_name] = dict()
                    if map_name not in self.young_masks: self.young_masks[map_name] = dict()

                    # Replace by a dictionary of maps
                    self.young_clipped_maps[map_name][levels_dict] = the_map

                    # Set the masks
                    self.young_masks[map_name][levels_dict] = mask

                # Doesn't have map and mask for each levels dict
                else: has_all = False

            # Remove other?
            if self.config.remove_other_data:
                log.debug("Removing the following files:")
                log.debug("")
                path = self.clipped_data_young_path_for_map(map_name, create=True)
                filenames = [fs.strip_extension(fs.name(filepath)) for filepath in required_filepaths]
                other_filepaths, other_filenames = fs.files_in_path(path, exact_not_name=filenames, returns=["path", "name"], unpack=True)
                for line in stringify_list_fancy(other_filenames, lines_prefix="  ")[1].split("\n"): log.debug(line)
                fs.remove_files(other_filepaths)
                log.debug("")

            # Set has all flag
            self.has_all_data_young[map_name] = has_all

            # If has all: remove the map
            if has_all and not self.config.reclip_from_masks: del self.young_maps[map_name]

            # Succes
            if has_all: log.success("All data is already present for the '" + map_name + "' young stellar map")

            # For reclipping
            if self.config.reclip_from_masks and has_any:

                # Message
                log.success("But reclipping will be applied")

                # Debugging
                log.debug("Correcting the '" + map_name + "' young stellar map to be reclipped ...")

                # Correct
                self.correct_map(self.young_maps[map_name])

                # Debugging
                log.debug("Rebinning the '" + map_name + "' young stellar map to be reclipped ...")

                # Rebin (so to be cropped as the clipped maps)
                last_levels_dict = self.young_clipped_maps[map_name].keys()[-1]
                self.young_maps[map_name].rebin(self.young_clipped_maps[map_name][last_levels_dict].wcs)

    # -----------------------------------------------------------------

    def get_sigma_levels_combinations_for_ionizing_map(self, name):

        """
        Thisnf unction ...
        :param name:
        :return:
        """

        # Get the origins
        #origins = self.ionizing_map_origins[name]
        origins = self.filtered_ionizing_map_origins[name]

        # Get image names
        names = [str(fltr) for fltr in origins]

        # Get the corresponding sigma levels
        sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

        # Make the combinations
        return names, sequences.lists_combinations(*sigma_levels)

    # -----------------------------------------------------------------

    def get_levels_dicts_for_ionizing_map(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Loop over the different combinations
        names, sigma_levels_combinations = self.get_sigma_levels_combinations_for_ionizing_map(name)

        # Loop over the different combinations
        for combination in sigma_levels_combinations:

            # Create dictionary that says which sigma level was used for which frame
            levels_dict = hashdict({name: level for name, level in zip(names, combination)})

            # Yield
            yield levels_dict

    # -----------------------------------------------------------------

    def check_data_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking existing data for ionizing stellar maps ...")

        # Loop over the maps
        for map_name in self.process_ionizing_map_names:

            # ALl plots already
            if self.has_all_ionizing_map_and_mask_plots(map_name): continue

            # Debugging
            log.debug("Checking data for the '" + map_name + "' ionizing stellar map ...")

            # Has all?
            has_all = True
            has_any = False

            # Initialize list of required filepaths
            required_filepaths = []

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_ionizing_map(map_name):

                # Already plotted
                if self.has_ionizing_map_plot_for_levels(map_name, levels_dict) and self.has_ionizing_mask_plot_for_levels(map_name, levels_dict): continue

                # Debugging
                levels_string = self.levels_to_string(levels_dict)
                log.debug("Checking data for the combination '" + levels_string + "' ...")

                # Debugging
                map_filepath, mask_filepath = self.get_ionizing_map_and_mask_paths_for_levels(map_name, levels_dict, create=False)
                log.debug("Looking for files: '" + map_filepath + "' and '" + mask_filepath + "' ...")

                # Add filepaths
                required_filepaths.append(map_filepath)
                required_filepaths.append(mask_filepath)

                # Check whether file is present
                if self.has_ionizing_map_and_mask_for_levels(map_name, levels_dict):

                    has_any = True

                    # Success
                    log.success("Data for combination '" + levels_string + "' is present")

                    # Load map and mask
                    the_map, mask = self.load_ionizing_map_and_mask_for_levels(map_name, levels_dict)

                    # Set flag
                    the_map.metadata["finished"] = False

                    # Initialize
                    if map_name not in self.ionizing_clipped_maps: self.ionizing_clipped_maps[map_name] = dict()
                    if map_name not in self.ionizing_masks: self.ionizing_masks[map_name] = dict()

                    # Replace by a dictionary of maps
                    self.ionizing_clipped_maps[map_name][levels_dict] = the_map

                    # Set the masks
                    self.ionizing_masks[map_name][levels_dict] = mask

                # Doesn't have map and mask for each levels dict
                else: has_all = False

            # Remove other?
            if self.config.remove_other_data:
                log.debug("Removing the following files:")
                log.debug("")
                path = self.clipped_data_ionizing_path_for_map(map_name, create=True)
                filenames = [fs.strip_extension(fs.name(filepath)) for filepath in required_filepaths]
                other_filepaths, other_filenames = fs.files_in_path(path, exact_not_name=filenames, returns=["path", "name"], unpack=True)
                for line in stringify_list_fancy(other_filenames, lines_prefix="  ")[1].split("\n"): log.debug(line)
                fs.remove_files(other_filepaths)
                log.debug("")

            # Set has all flag
            self.has_all_data_ionizing[map_name] = has_all

            # If has all: remove the map
            if has_all and not self.config.reclip_from_masks: del self.ionizing_maps[map_name]

            # Succes
            if has_all: log.success("All data is already present for the '" + map_name + "' ionizing stellar map")

            # For reclipping
            if self.config.reclip_from_masks and has_any:

                # Message
                log.success("But reclipping will be applied")

                # Debugging
                log.debug("Correcting the '" + map_name + "' ionizing stellar map to be reclipped ...")

                # Correct
                self.correct_map(self.ionizing_maps[map_name])

                # Debugging
                log.debug("Rebinning the '" + map_name + "' ionizing stellar map to be reclipped ...")

                # Rebin (so to be cropped as the clipped maps)
                last_levels_dict = self.ionizing_clipped_maps[map_name].keys()[-1]
                self.ionizing_maps[map_name].rebin(self.ionizing_clipped_maps[map_name][last_levels_dict].wcs)

    # -----------------------------------------------------------------

    def get_sigma_levels_combinations_for_dust_map(self, name):

        """
        Thisnf unction ...
        :param name:
        :return:
        """

        # Get the origins
        #origins = self.dust_map_origins[name]
        origins = self.filtered_dust_map_origins[name]

        # Get image names
        names = [str(fltr) for fltr in origins]

        # Get the corresponding sigma levels
        sigma_levels = self.get_sigma_levels_for_origins(origins, as_list=True)

        # Make the combinations
        return names, sequences.lists_combinations(*sigma_levels)

    # -----------------------------------------------------------------

    def get_levels_dicts_for_dust_map(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Loop over the different combinations
        names, sigma_levels_combination = self.get_sigma_levels_combinations_for_dust_map(name)

        # Loop over the different combinations
        for combination in sigma_levels_combination:

            # Create dictionary that says which sigma level was used for which frame
            levels_dict = hashdict({name: level for name, level in zip(names, combination)})

            # Yield
            yield levels_dict

    # -----------------------------------------------------------------

    def check_data_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking existing data for dust maps ...")

        # Loop over the maps
        for map_name in self.process_dust_map_names:

            # ALl plots already
            if self.has_all_dust_map_and_mask_plots(map_name): continue

            # Debugging
            log.debug("Checking data for the '" + map_name + "' dust map ...")

            # Has all?
            has_all = True
            has_any = False

            # Initialize list of required filepaths
            required_filepaths = []

            # Loop over the levels dictionaries
            for levels_dict in self.get_levels_dicts_for_dust_map(map_name):

                # Already plotted
                if self.has_dust_map_plot_for_levels(map_name, levels_dict) and self.has_dust_mask_plot_for_levels(map_name, levels_dict): continue

                # Debugging
                levels_string = self.levels_to_string(levels_dict)
                log.debug("Checking data for the combination '" + levels_string + "' ...")

                # Debugging
                map_filepath, mask_filepath = self.get_dust_map_and_mask_paths_for_levels(map_name, levels_dict, create=False)
                log.debug("Looking for files: '" + map_filepath + "' and '" + mask_filepath + "' ...")

                # Add filepaths
                required_filepaths.append(map_filepath)
                required_filepaths.append(mask_filepath)

                # Check whether file is present
                if self.has_dust_map_and_mask_for_levels(map_name, levels_dict):

                    has_any = True

                    # Success
                    log.success("Data for combination '" + levels_string + "' is present")

                    # Load map and mask
                    the_map, mask = self.load_dust_map_and_mask_for_levels(map_name, levels_dict)

                    # Set flag
                    the_map.metadata["finished"] = False

                    # Initialize
                    if map_name not in self.dust_clipped_maps: self.dust_clipped_maps[map_name] = dict()
                    if map_name not in self.dust_masks: self.dust_masks[map_name] = dict()

                    # Replace by a dictionary of maps
                    self.dust_clipped_maps[map_name][levels_dict] = the_map

                    # Set the masks
                    self.dust_masks[map_name][levels_dict] = mask

                # Doesn't have map and mask for each levels dict
                else: has_all = False

            # Remove other?
            if self.config.remove_other_data:
                log.debug("Removing the following files:")
                log.debug("")
                path = self.clipped_data_dust_path_for_map(map_name, create=True)
                filenames = [fs.strip_extension(fs.name(filepath)) for filepath in required_filepaths]
                other_filepaths, other_filenames = fs.files_in_path(path, exact_not_name=filenames, returns=["path", "name"], unpack=True)
                for line in stringify_list_fancy(other_filenames, lines_prefix="  ")[1].split("\n"): log.debug(line)
                fs.remove_files(other_filepaths)
                log.debug("")

            # Set has all flag
            self.has_all_data_dust[map_name] = has_all

            # If has all: remove the map
            if has_all and not self.config.reclip_from_masks: del self.dust_maps[map_name]

            # Succes
            if has_all: log.success("All data is already present for the '" + map_name + "' dust stellar map")

            # For reclipping
            if self.config.reclip_from_masks and has_any:

                # Message
                log.success("But reclipping will be applied")

                # Debugging
                log.debug("Correcting the '" + map_name + "' dust map to be reclipped ...")

                # Correct
                self.correct_map(self.dust_maps[map_name])

                # Debugging
                log.debug("Rebinning the '" + map_name + "' dust map to be reclipped ...")

                # Rebin (so to be cropped as the clipped maps)
                last_levels_dict = self.dust_clipped_maps[map_name].keys()[-1]
                self.dust_maps[map_name].rebin(self.dust_clipped_maps[map_name][last_levels_dict].wcs)

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
        if self.do_old: self.correct_old_maps()

        # Young
        if self.do_young: self.correct_young_maps()

        # Ionizing
        if self.do_ionizing: self.correct_ionizing_maps()

        # Dust
        if self.do_dust: self.correct_dust_maps()

    # -----------------------------------------------------------------

    def correct_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the old stellar maps ...")

        # Loop over the maps
        for name in self.process_old_map_names:

            # Needed?
            if name not in self.old_maps:
                if not self.has_all_data_old[name]: raise ValueError("Something went wrong")
                if not len(self.old_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_old[name]: continue

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
        for name in self.process_young_map_names:

            # Needed?
            if name not in self.young_maps:
                if not self.has_all_data_young[name]: raise ValueError("Something went wrong")
                if not len(self.young_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_young[name]: continue

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
        for name in self.process_ionizing_map_names:

            # Needed?
            if name not in self.ionizing_maps:
                if not self.has_all_data_ionizing[name]: raise ValueError("Something went wrong")
                if not len(self.ionizing_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already corrected?
            if self.has_all_data_ionizing[name]: continue

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
        for name in self.process_dust_map_names:

            # Needed?
            if name not in self.dust_maps:
                if not self.has_all_data_dust[name]: raise ValueError("Something went wrong")
                if not len(self.dust_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already corrected?
            if self.has_all_data_dust[name]: continue

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
        if self.do_old: self.crop_old_maps()

        # Young
        if self.do_young: self.crop_young_maps()

        # Ionizing
        if self.do_ionizing: self.crop_ionizing_maps()

        # Dust
        if self.do_dust: self.crop_dust_maps()

    # -----------------------------------------------------------------

    def crop_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the old stellar maps ...")

        # Loop over the maps
        for name in self.process_old_map_names:

            # Needed?
            if name not in self.old_maps:
                if not self.has_all_data_old[name]: raise ValueError("Something went wrong")
                if not len(self.old_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_old[name]: continue

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
        for name in self.process_young_map_names:

            # Needed?
            if name not in self.young_maps:
                if not self.has_all_data_young[name]: raise ValueError("Something went wrong")
                if not len(self.young_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_young[name]: continue

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
        for name in self.process_ionizing_map_names:

            # Needed?
            if name not in self.ionizing_maps:
                if not self.has_all_data_ionizing[name]: raise ValueError("Something went wrong")
                if not len(self.ionizing_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_ionizing[name]: continue

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
        for name in self.process_dust_map_names:

            # Needed?
            if name not in self.dust_maps:
                if not self.has_all_data_dust[name]: raise ValueError("Something went wrong")
                if not len(self.dust_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Not needed because the map is already cropped?
            if self.has_all_data_dust[name]: continue

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
        if self.do_old: self.clip_old_maps()

        # Young
        if self.do_young: self.clip_young_maps()

        # Ionizing
        if self.do_ionizing: self.clip_ionizing_maps()

        # Dust
        if self.do_dust: self.clip_dust_maps()

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

    @memoize_method
    def has_all_old_plots_for(self, name):

        """
        This function ...
        :param name:
        :param origins:
        :param sigma_levels:
        :return:
        """

        return self.has_all_old_map_plots_for(name) and self.has_all_old_mask_plots_for(name)

    # -----------------------------------------------------------------

    def has_all_old_map_plots_for(self, name):

        """
        This function ...
        :param name:
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

    def has_all_old_mask_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.old_mask_plot_paths[name]:

            # Get the filepath
            filepath = self.old_mask_plot_paths[name][levels]

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

    @memoize_method
    def has_all_young_plots_for(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return self.has_all_young_map_plots_for(name) and self.has_all_young_mask_plots_for(name)

    # -----------------------------------------------------------------

    def has_all_young_map_plots_for(self, name):

        """
        This function ...
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

    def has_all_young_mask_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.young_mask_plot_paths[name]:

            # Get the filepath
            filepath = self.young_mask_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): return False

        # ALl checks passed
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

    @memoize_method
    def has_all_ionizing_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_all_ionizing_map_plots_for(name) and self.has_all_ionizing_mask_plots_for(name)

    # -----------------------------------------------------------------

    def has_all_ionizing_map_plots_for(self, name):

        """
        Thisf unction ...
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

    def has_all_ionizing_mask_plots_for(self, name):

        """
        This fucntion ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.ionizing_mask_plot_paths[name]:

            # Get the filepath
            filepath = self.ionizing_mask_plot_paths[name][levels]

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

    @memoize_method
    def has_all_dust_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_all_dust_map_plots_for(name) and self.has_all_dust_mask_plots_for(name)

    # -----------------------------------------------------------------

    def has_all_dust_map_plots_for(self, name):

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

    def has_all_dust_mask_plots_for(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the level combinations
        for levels in self.dust_mask_plot_paths[name]:

            # Get the filepath
            filepath = self.dust_mask_plot_paths[name][levels]

            # Check
            if not fs.is_file(filepath): continue

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def old_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.old_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_settings_old(self, name, current=None, current_masks=None):

        """
        This function ...
        :param name:
        :param current:
        :param current_masks:
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["present"] = self.present_old_plots_level_combinations_for(name)
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["return_masks"] = True
        settings["current"] = current
        settings["current_masks"] = current_masks
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["resoften_current_masks"] = self.config.resoften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_old
        settings["boundary"] = self.old_maps_boundary
        settings["plot"] = self.config.plot_clipping_old
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_old

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the old stellar maps ...")

        # Loop over the maps
        for name in self.process_old_map_names:

            # Needed?
            if name not in self.old_maps:
                if not self.has_all_data_old[name]: raise ValueError("Something went wrong")
                if not len(self.old_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' old stellar map ...")

            # Get the origins
            #origins = self.old_map_origins[name]
            origins = self.filtered_old_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Get current maps to reuse
            if self.config.reclip_from_masks: current = None
            else: current = self.old_clipped_maps[name] if name in self.old_clipped_maps else None

            # Get current masks
            current_masks = self.old_masks[name] if name in self.old_masks else None

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.old_maps[name], origins, sigma_levels, **self.get_clip_settings_old(name, current, current_masks))

            # Replace by a dictionary of maps
            self.old_clipped_maps[name] = maps

            # Set the masks
            self.old_masks[name] = masks

            # Write the data
            if self.config.write_data: self.write_old_maps_and_masks(name, maps, masks)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def get_old_map_and_mask_paths_for_levels(self, name, levels_dict, create=True):

        """
        Thisf unction ....
        :param name:
        :param levels_dict:
        :param create:
        :return:
        """

        # Get filenames
        map_filename, mask_filename = self.get_map_and_masks_filenames_for_levels(name, levels_dict)

        # Determine filepaths
        path = self.clipped_data_old_path_for_map(name, create=create)
        map_filepath = fs.join(path, map_filename)
        mask_filepath = fs.join(path, mask_filename)

        # Return the filepaths
        return map_filepath, mask_filepath

    # -----------------------------------------------------------------

    def has_old_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get the filepaths
        map_filepath, mask_filepath = self.get_old_map_and_mask_paths_for_levels(name, levels_dict, create=False)

        # Return
        return fs.is_file(map_filepath) and fs.is_file(mask_filepath)

    # -----------------------------------------------------------------

    def load_old_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get the filepaths
        map_filepath, mask_filepath = self.get_old_map_and_mask_paths_for_levels(name, levels_dict)

        # Return
        return Frame.from_file(map_filepath, no_filter=True), load_mask_or_alpha_mask(mask_filepath, no_wcs=True)

    # -----------------------------------------------------------------

    def write_old_maps_and_masks(self, name, maps, masks):

        """
        This function ...
        :param name:
        :param maps:
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Writing the clipped old stellar maps and clip masks for the '" + name + "' map ...")

        # Loop over the level combinations
        for levels_dict in maps:

            # Get filepaths
            map_filepath, mask_filepath = self.get_old_map_and_mask_paths_for_levels(name, levels_dict)

            # Save
            maps[levels_dict].saveto(map_filepath)
            masks[levels_dict].saveto(mask_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.young_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_settings_young(self, name, current=None, current_masks=None):

        """
        This function ...
        :param name:
        :param current:
        :param current_masks:
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["present"] = self.present_young_plots_level_combinations_for(name)
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["return_masks"] = True
        settings["current"] = current
        settings["current_masks"] = current_masks
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["resoften_current_masks"] = self.config.resoften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_young
        settings["boundary"] = self.young_maps_boundary
        settings["plot"] = self.config.plot_clipping_young
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_young

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the young stellar maps ...")

        # Loop over the maps
        for name in self.process_young_map_names:

            # Needed?
            if name not in self.young_maps:
                if not self.has_all_data_young[name]: raise ValueError("Something went wrong")
                if not len(self.young_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            #origins = self.young_map_origins[name]
            origins = self.filtered_young_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Get current maps to reuse
            if self.config.reclip_from_masks: current = None
            else: current = self.young_clipped_maps[name] if name in self.young_clipped_maps else None

            # Get current masks
            current_masks = self.young_masks[name] if name in self.young_masks else None

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.young_maps[name], origins, sigma_levels, **self.get_clip_settings_young(name, current, current_masks))

            # Replace by a dictionary of maps
            self.young_clipped_maps[name] = maps

            # Set the masks
            self.young_masks[name] = masks

            # Write the data
            if self.config.write_data: self.write_young_maps_and_masks(name, maps, masks)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def get_young_map_and_mask_paths_for_levels(self, name, levels_dict, create=True):

        """
        Thisf unction ....
        :param name:
        :param levels_dict:
        :param create:
        :return:
        """

        # Get filenames
        map_filename, mask_filename = self.get_map_and_masks_filenames_for_levels(name, levels_dict)

        # Determine filepaths
        path = self.clipped_data_young_path_for_map(name, create=create)
        map_filepath = fs.join(path, map_filename)
        mask_filepath = fs.join(path, mask_filename)

        # Return the filepaths
        return map_filepath, mask_filepath

    # -----------------------------------------------------------------

    def has_young_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_young_map_and_mask_paths_for_levels(name, levels_dict, create=False)

        # Check and return
        return fs.is_file(map_filepath) and fs.is_file(mask_filepath)

    # -----------------------------------------------------------------

    def load_young_map_and_mask_for_levels(self, name, levels_dict):

        """
        Thisnf unction ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_young_map_and_mask_paths_for_levels(name, levels_dict)

        # Return
        return Frame.from_file(map_filepath, no_filter=True), load_mask_or_alpha_mask(mask_filepath, no_wcs=True)

    # -----------------------------------------------------------------

    def write_young_maps_and_masks(self, name, maps, masks):

        """
        This function ...
        :param name:
        :param maps:
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Writing the clipped young stellar maps and clip masks for the '" + name + "' map ...")

        # Loop over the level combinations
        for levels_dict in maps:

            # Get filepaths
            map_filepath, mask_filepath = self.get_young_map_and_mask_paths_for_levels(name, levels_dict)

            # Save
            maps[levels_dict].saveto(map_filepath)
            masks[levels_dict].saveto(mask_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.ionizing_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_settings_ionizing(self, name, current=None, current_masks=None):

        """
        This function ...
        :param name:
        :param current:
        :param current_masks:
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["present"] = self.present_ionizing_plots_level_combinations_for(name)
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["return_masks"] = True
        settings["current"] = current
        settings["current_masks"] = current_masks
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks,
        settings["resoften_current_masks"] = self.config.resoften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_ionizing
        settings["boundary"] = self.config.ionizing_maps_boundary
        settings["plot"] = self.config.plot_clipping_ionizing
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_ionizing

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.process_ionizing_map_names:

            # Needed?
            if name not in self.ionizing_maps:
                if not self.has_all_data_ionizing[name]: raise ValueError("Something went wrong")
                if not len(self.ionizing_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            #origins = self.ionizing_map_origins[name]
            origins = self.filtered_ionizing_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Get current maps to reuse
            if self.config.reclip_from_masks: current = None
            else: current = self.ionizing_clipped_maps[name] if name in self.ionizing_clipped_maps else None

            # Get current masks
            current_masks = self.ionizing_masks[name] if name in self.ionizing_masks else None

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.ionizing_maps[name], origins, sigma_levels, **self.get_clip_settings_ionizing(name, current, current_masks))

            # Replace by a dictionary of maps
            self.ionizing_clipped_maps[name] = maps

            # Set the masks
            self.ionizing_masks[name] = masks

            # Write the data
            if self.config.write_data: self.write_ionizing_maps_and_masks(name, maps, masks)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def get_ionizing_map_and_mask_paths_for_levels(self, name, levels_dict, create=True):

        """
        Thisf unction ....
        :param name:
        :param levels_dict:
        :param create:
        :return:
        """

        # Get filenames
        map_filename, mask_filename = self.get_map_and_masks_filenames_for_levels(name, levels_dict)

        # Determine filepaths
        path = self.clipped_data_ionizing_path_for_map(name, create=create)
        map_filepath = fs.join(path, map_filename)
        mask_filepath = fs.join(path, mask_filename)

        # Return the filepaths
        return map_filepath, mask_filepath

    # -----------------------------------------------------------------

    def has_ionizing_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_ionizing_map_and_mask_paths_for_levels(name, levels_dict, create=False)

        # Check and return
        return fs.is_file(map_filepath) and fs.is_file(mask_filepath)

    # -----------------------------------------------------------------

    def load_ionizing_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_ionizing_map_and_mask_paths_for_levels(name, levels_dict)

        # Return
        return Frame.from_file(map_filepath, no_filter=True), load_mask_or_alpha_mask(mask_filepath, no_wcs=True)

    # -----------------------------------------------------------------

    def write_ionizing_maps_and_masks(self, name, maps, masks):

        """
        This function ...
        :param name:
        :param maps:
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Writing the clipped ionizing stellar maps and clip masks for the '" + name + "' map ...")

        # Loop over the level combinations
        for levels_dict in maps:

            # Get filepaths
            map_filepath, mask_filepath = self.get_ionizing_map_and_mask_paths_for_levels(name, levels_dict)

            # Save
            maps[levels_dict].saveto(map_filepath)
            masks[levels_dict].saveto(mask_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.dust_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_settings_dust(self, name, current=None, current_masks=None):

        """
        This function ...
        :param name:
        :param current:
        :param current_masks:
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["present"] = self.present_dust_plots_level_combinations_for(name)
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["return_masks"] = True
        settings["current"] = current
        settings["current_masks"] = current_masks
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["resoften_current_masks"] = self.config.resoften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_dust
        settings["boundary"] = self.dust_maps_boundary
        settings["plot"] = self.config.plot_clipping_dust
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_dust

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the dust maps ...")

        # Loop over the maps
        for name in self.process_dust_map_names:

            # Needed?
            if name not in self.dust_maps:
                if not self.has_all_data_dust[name]: raise ValueError("Something went wrong")
                if not len(self.dust_clipped_maps[name]) > 0: raise ValueError("Something went wrong")
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            #origins = self.dust_map_origins[name]
            origins = self.filtered_dust_map_origins[name]

            # Get the sigma levels
            sigma_levels = self.get_sigma_levels_for_origins(origins)

            # Debugging
            log.debug("The relevant sigma levels are: " + tostr(sigma_levels))

            # Get current maps to reuse
            if self.config.reclip_from_masks: current = None
            else: current = self.dust_clipped_maps[name] if name in self.dust_clipped_maps else None

            # Get current masks to use
            current_masks = self.dust_masks[name] if name in self.dust_masks else None

            # Clip the map
            maps, masks = self.make_clipped_maps(name, self.dust_maps[name], origins, sigma_levels, **self.get_clip_settings_dust(name, current, current_masks))

            # Replace by a dictionary of maps
            self.dust_clipped_maps[name] = maps

            # Set the masks
            self.dust_masks[name] = masks

            # Write the data
            if self.config.write_data: self.write_dust_maps_and_masks(name, maps, masks)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def get_dust_map_and_mask_paths_for_levels(self, name, levels_dict, create=True):

        """
        Thisf unction ....
        :param name:
        :param levels_dict:
        :param create:
        :return:
        """

        # Get filenames
        map_filename, mask_filename = self.get_map_and_masks_filenames_for_levels(name, levels_dict)

        # Determine filepaths
        path = self.clipped_data_dust_path_for_map(name, create=create)
        map_filepath = fs.join(path, map_filename)
        mask_filepath = fs.join(path, mask_filename)

        # Return the filepaths
        return map_filepath, mask_filepath

    # -----------------------------------------------------------------

    def has_dust_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_dust_map_and_mask_paths_for_levels(name, levels_dict, create=False)

        # Return
        return fs.is_file(map_filepath) and fs.is_file(mask_filepath)

    # -----------------------------------------------------------------

    def load_dust_map_and_mask_for_levels(self, name, levels_dict):

        """
        This function ...
        :param name:
        :param levels_dict:
        :return:
        """

        # Get filepaths
        map_filepath, mask_filepath = self.get_dust_map_and_mask_paths_for_levels(name, levels_dict)

        # Return
        return Frame.from_file(map_filepath, no_filter=True), load_mask_or_alpha_mask(mask_filepath, no_wcs=True)

    # -----------------------------------------------------------------

    def write_dust_maps_and_masks(self, name, maps, masks):

        """
        This function ...
        :param name:
        :param maps:
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Writing the clipped dust maps and clip masks for the '" + name + "' map ...")

        # Loop over the level combinations
        for levels_dict in maps:

            # Get filepaths
            map_filepath, mask_filepath = self.get_dust_map_and_mask_paths_for_levels(name, levels_dict)

            # Save
            maps[levels_dict].saveto(map_filepath)
            masks[levels_dict].saveto(mask_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def filtered_old_map_origins(self):

        """
        This function ...
        :return:
        """

        # No filters to be ignored
        if self.config.ignore_filters is None: return self.old_map_origins

        # Some filters should be ignored
        else:

            origins = dict()
            for name in self.old_map_origins:
                origins_name = sequences.removed(self.old_map_origins[name], self.config.ignore_filters)
                origins[name] = origins_name
            return origins

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

            origins = self.filtered_old_map_origins[map_name]
            filter_names = [str(origin) for origin in origins]
            names.update(filter_names)

        # Return the list of image names
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def filtered_young_map_origins(self):

        """
        This function ...
        :return:
        """

        # No filters to be ignored
        if self.config.ignore_filters is None: return self.young_map_origins

        # Some filters should be ignored
        else:

            origins = dict()
            for name in self.young_map_origins:
                origins_name = sequences.removed(self.young_map_origins[name], self.config.ignore_filters)
                origins[name] = origins_name
            return origins

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

            origins = self.filtered_young_map_origins[map_name]
            filter_names = [str(origin) for origin in origins]
            names.update(filter_names)

        # Return the list of image names
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def filtered_ionizing_map_origins(self):

        """
        Thisn function ...
        :return:
        """

        # No filters to be ignored
        if self.config.ignore_filters is None: return self.ionizing_map_origins

        # Some filters should be ignored
        else:

            origins = dict()
            for name in self.ionizing_map_origins:
                origins_name = sequences.removed(self.ionizing_map_origins[name], self.config.ignore_filters)
                origins[name] = origins_name
            return origins

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

            origins = self.filtered_ionizing_map_origins[map_name]
            filter_names = [str(origin) for origin in origins]
            names.update(filter_names)

        # Return the list of image names
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def filtered_dust_map_origins(self):

        """
        This function ...
        :return:
        """

        # No filters to be ignored
        if self.config.ignore_filters is None: return self.dust_map_origins

        # Some filters should be ignored
        else:

            origins = dict()
            for name in self.dust_map_origins:
                origins_name = sequences.removed(self.dust_map_origins[name], self.config.ignore_filters)
                origins[name] = origins_name
            return origins

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

            origins = self.filtered_dust_map_origins[map_name]
            filter_names = [str(origin) for origin in origins]
            names.update(filter_names)

        # Return the list of image names
        return list(names)

    # -----------------------------------------------------------------

    @lazyproperty
    def image_names(self):

        """
        This function ...
        :return:
        """

        names = set()
        #if self.config.add_old: names.update(self.old_maps_image_names)
        #if self.config.add_young: names.update(self.young_maps_image_names)
        #if self.config.add_ionizing: names.update(self.ionizing_maps_image_names)
        #if self.config.add_dust: names.update(self.dust_maps_image_names)
        if self.do_old: names.update(self.old_maps_image_names)
        if self.do_young: names.update(self.young_maps_image_names)
        if self.do_ionizing: names.update(self.ionizing_maps_image_names)
        if self.do_dust: names.update(self.dust_maps_image_names)
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

    def make_rgba_plot(self, the_map, filepath, colours="red", scale="log"):

        """
        This function ...
        :param the_map:
        :param filepath:
        :param colours:
        :param scale:
        :return:
        """

        # Plot as RGB
        the_map.saveto_png(filepath, interval=self.config.interval, scale=scale,
                           alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha,
                           colours=colours)

    # -----------------------------------------------------------------

    def make_mask_plot(self, mask, filepath):

        """
        This function ...
        :param mask:
        :param filepath:
        :return:
        """

        if isinstance(mask, AlphaMask): mask.saveto_png(filepath, colour=self.config.mask_colour, alpha=True) # alpha to recognize alpha masks
        elif isinstance(mask, Mask): mask.saveto_png(filepath, colour=self.config.mask_colour, alpha=False) # no alpha to recognize regular masks
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
        if self.do_old: self.make_old_plots()

        # Young
        if self.do_young: self.make_young_plots()

        # Ionizing
        if self.do_ionizing: self.make_ionizing_plots()

        # Dust
        if self.do_dust: self.make_dust_plots()

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

    def get_old_map_filepath(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.old_plot_paths[name][levels]
        except KeyError:

            if len(self.old_plot_paths[name]) == 0: raise RuntimeError("No levels in the subdirectory")
            else:

                log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
                log.warning("and dictionary keys:")
                for key in self.old_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

                log.warning("Looping over the keys and checking for equivalency ...")

                filepath = None
                for key in self.old_plot_paths[name]:
                    if containers.close_dicts(levels, key): filepath = self.old_plot_paths[name][key]
                    break

                else:
                    for key in self.old_plot_paths[name]:
                        differences = containers.dicts_differences(key, levels)
                        log.error(" - {" + tostr(differences) + "}")
                    raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def has_old_map_plot_for_levels(self, name, levels):

        """
        Thisf unction ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_old_map_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_old_mask_plot_for_levels(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_old_mask_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_all_old_map_plots(self, name):


        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_old_map(name):

            # Check
            if not self.has_old_map_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_old_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_old_map(name):

            # Check
            if not self.has_old_mask_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_old_map_and_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_all_old_map_plots(name) and self.has_all_old_mask_plots(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_all_old_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stelar maps
        for name in self.old_selection:

            # Check map
            if not self.has_all_old_map_plots(name): return False

            # Check mask
            if not self.has_all_old_mask_plots(name): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def make_old_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_old_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' old stellar map ...")

            # Loop over the levels dicts
            for levels in self.old_clipped_maps[name]:

                # Get the filepath
                filepath = self.get_old_map_filepath(name, levels)

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.old_clipped_maps[name][levels], filepath, colours=self.old_color, scale=self.old_scale)

            # Clear
            del self.old_clipped_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def get_old_mask_filepath(self, name, levels):

        """
        This fucntion ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.old_mask_plot_paths[name][levels]
        except KeyError:

            log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
            log.warning("and dictionary keys:")
            for key in self.old_mask_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

            log.warning("Looping over the keys and checking for equivalency ...")

            filepath = None
            for key in self.old_mask_plot_paths[name]:
                if containers.close_dicts(levels, key): filepath = self.old_mask_plot_paths[name][key]
                break

            else:
                for key in self.old_mask_plot_paths[name]:
                    differences = containers.dicts_differences(key, levels)
                    log.error(" - {" + tostr(differences) + "}")
                raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def make_old_mask_plots(self):

        """
        This function ...
        :return: 
        """

        # Loop over the maps
        for name in self.process_old_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' old stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.old_masks[name]:

                # Get the filepath
                filepath = self.get_old_mask_filepath(name, levels)

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

    def get_young_map_filepath(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        # Determine the filepath
        try: filepath = self.young_plot_paths[name][levels]
        except KeyError:

            if len(self.young_plot_paths[name]) == 0: raise RuntimeError("No levels in the subdirectory")
            else:

                log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
                log.warning("and dictionary keys:")
                for key in self.young_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

                log.warning("Looping over the keys and checking for equivalency ...")

                filepath = None
                for key in self.young_plot_paths[name]:
                    if containers.close_dicts(levels, key):
                        filepath = self.young_plot_paths[name][key]
                        break
                else:
                    for key in self.young_plot_paths[name]:
                        differences = containers.dicts_differences(key, levels)
                        log.error(" - {" + tostr(differences) + "}")
                    raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def has_young_map_plot_for_levels(self, name, levels):

        """
        Thisf unction ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_young_map_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_young_mask_plot_for_levels(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_young_mask_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_all_young_map_plots(self, name):


        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_young_map(name):

            # Check
            if not self.has_young_map_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_young_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_young_map(name):

            # Check
            if not self.has_young_mask_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_young_map_and_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_all_young_map_plots(name) and self.has_all_young_mask_plots(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_all_young_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stelar maps
        for name in self.young_selection:

            # Check map
            if not self.has_all_young_map_plots(name): return False

            # Check mask
            if not self.has_all_young_mask_plots(name): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def make_young_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_young_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' young stellar map ...")

            # Loop over the levels dicts
            #for levels in self.young_maps[name]:
            for levels in self.young_clipped_maps[name]:

                # Get filepath
                filepath = self.get_young_map_filepath(name, levels)

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.young_clipped_maps[name][levels], filepath, colours=self.young_color, scale=self.young_scale)

            # Clear
            del self.young_clipped_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def get_young_mask_filepath(self, name, levels):

        """
        Thisf unction ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.young_mask_plot_paths[name][levels]
        except KeyError:

            log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
            log.warning("and dictionary keys:")
            for key in self.young_mask_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

            log.warning("Looping over the keys and checking for equivalency ...")

            filepath = None
            for key in self.young_mask_plot_paths[name]:
                if containers.close_dicts(levels, key):
                    filepath = self.young_mask_plot_paths[name][key]
                    break

            else:
                for key in self.young_mask_plot_paths[name]:
                    differences = containers.dicts_differences(key, levels)
                    log.error(" - {" + tostr(differences) + "}")
                raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def make_young_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the masp
        for name in self.process_young_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' young stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.young_masks[name]:

                # Get the filepath
                filepath = self.get_young_mask_filepath(name, levels)

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

    def get_ionizing_map_filepath(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.ionizing_plot_paths[name][levels]
        except KeyError:

            if len(self.ionizing_plot_paths[name]) == 0: raise RuntimeError("No levels in subdirectory")
            else:

                log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
                log.warning("and dictionary keys:")
                for key in self.ionizing_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

                log.warning("Looping over the keys and checking for equivalency ...")

                filepath = None
                for key in self.ionizing_plot_paths[name]:
                    if containers.close_dicts(levels, key): filepath = self.ionizing_plot_paths[name][key]
                    break

                else:
                    for key in self.ionizing_plot_paths[name]:
                        differences = containers.dicts_differences(key, levels)
                        log.error(" - {" + tostr(differences) + "}")
                    raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def has_ionizing_map_plot_for_levels(self, name, levels):

        """
        Thisf unction ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_ionizing_map_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_ionizing_mask_plot_for_levels(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_ionizing_mask_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_all_ionizing_map_plots(self, name):


        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_ionizing_map(name):

            # Check
            if not self.has_ionizing_map_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_ionizing_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_ionizing_map(name):

            # Check
            if not self.has_ionizing_mask_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_ionizing_map_and_mask_plots(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return self.has_all_ionizing_map_plots(name) and self.has_all_ionizing_mask_plots(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_all_ionizing_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stelar maps
        for name in self.ionizing_selection:

            # Check map
            if not self.has_all_ionizing_map_plots(name): return False

            # Check mask
            if not self.has_all_ionizing_mask_plots(name): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def make_ionizing_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_ionizing_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' ionizing stellar map ...")

            # Loop over the levels dicts
            #for levels in self.ionizing_maps[name]:
            for levels in self.ionizing_clipped_maps[name]:

                # Get the filepath
                filepath = self.get_ionizing_map_filepath(name, levels)

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.ionizing_clipped_maps[name][levels], filepath, colours=self.ionizing_color, scale=self.ionizing_scale)

            # Clear
            del self.ionizing_clipped_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def get_ionizing_mask_filepath(self, name, levels):

        """
        Thisn function ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.ionizing_mask_plot_paths[name][levels]
        except KeyError:

            log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
            log.warning("and dictionary keys:")
            for key in self.ionizing_mask_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

            log.warning("Looping over the keys and checking for equivalency ...")

            filepath = None
            for key in self.ionizing_mask_plot_paths[name]:
                if containers.close_dicts(levels, key): filepath = self.ionizing_mask_plot_paths[name][key]
                break

            else:
                for key in self.ionizing_mask_plot_paths[name]:
                    differences = containers.dicts_differences(key, levels)
                    log.error(" - {" + tostr(differences) + "}")
                raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def make_ionizing_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_ionizing_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' ionizing stellar map masks ...")

            # Loop over the levels dicts
            for levels in self.ionizing_masks[name]:

                # Get the filepath
                filepath = self.get_ionizing_mask_filepath(name, levels)

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

    def get_dust_map_filepath(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.dust_plot_paths[name][levels]
        except KeyError:

            if len(self.dust_plot_paths[name]) == 0: raise RuntimeError("No levels in subdictionary")
            else:

                log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
                log.warning("and dictionary keys:")
                for key in self.dust_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

                log.warning("Looping over the keys and checking for equivalency ...")

                filepath = None
                for key in self.dust_plot_paths[name]:
                    if containers.close_dicts(levels, key): filepath = self.dust_plot_paths[name][key]
                    break

                else:
                    for key in self.dust_plot_paths[name]:
                        differences = containers.dicts_differences(key, levels)
                        log.error(" - {" + tostr(differences) + "}")
                    raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def has_dust_map_plot_for_levels(self, name, levels):

        """
        Thisf unction ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_dust_map_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_dust_mask_plot_for_levels(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        # Get the filepath
        filepath = self.get_dust_mask_filepath(name, levels)

        # Check
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_all_dust_map_plots(self, name):


        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_dust_map(name):

            # Check
            if not self.has_dust_map_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_dust_mask_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the levels dicts
        for levels_dict in self.get_levels_dicts_for_dust_map(name):

            # Check
            if not self.has_dust_mask_plot_for_levels(name, levels_dict): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def has_all_dust_map_and_mask_plots(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return self.has_all_dust_map_plots(name) and self.has_all_dust_mask_plots(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_all_dust_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stelar maps
        for name in self.dust_selection:

            # Check map
            if not self.has_all_dust_map_plots(name): return False

            # Check mask
            if not self.has_all_dust_mask_plots(name): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def make_dust_map_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_dust_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' dust map ...")

            # Loop over the levels dicts
            #for levels in self.dust_maps[name]:
            for levels in self.dust_clipped_maps[name]:

                # Get the filepath
                filepath = self.get_dust_map_filepath(name, levels)

                # Check
                if fs.is_file(filepath): continue

                # Debugging
                log.debug(" - Making plots for the levels {" + tostr(levels) + "} ...")

                # Save as RGBA
                self.make_rgba_plot(self.dust_clipped_maps[name][levels], filepath, colours=self.dust_color, scale=self.dust_scale)

            # Clear
            del self.dust_clipped_maps[name]
            gc.collect()

    # -----------------------------------------------------------------

    def get_dust_mask_filepath(self, name, levels):

        """
        This function ...
        :param name:
        :param levels:
        :return:
        """

        try: filepath = self.dust_mask_plot_paths[name][levels]
        except KeyError:

            log.warning("Mismatch between key = {" + tostr(levels) + "} (" + str(type(levels)) + ")")
            log.warning("and dictionary keys:")
            for key in self.dust_mask_plot_paths[name]: log.warning(" - {" + tostr(key) + "} (" + str(type(key)) + ")")

            log.warning("Looping over the keys and checking for equivalency ...")

            filepath = None
            for key in self.dust_mask_plot_paths[name]:
                if containers.close_dicts(levels, key): filepath = self.dust_mask_plot_paths[name][key]
                break

            else:
                for key in self.dust_mask_plot_paths[name]:
                    differences = containers.dicts_differences(key, levels)
                    log.error(" - {" + tostr(differences) + "}")
                raise RuntimeError("Could not find an equivalent key")

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def make_dust_mask_plots(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.process_dust_map_names:

            # Debugging
            log.debug("Making plots of the '" + name + "' dust map masks ...")

            # Loop over the levels dicts
            for levels in self.dust_masks[name]:

                # Get the filepath
                filepath = self.get_dust_mask_filepath(name, levels)

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

    @lazyproperty
    def has_all_plots(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.has_all_old_plots and self.has_all_young_plots and self.has_all_ionizing_plots and self.has_all_dust_plots

    # -----------------------------------------------------------------

    def get_sigma_levels_for_origins(self, origins, as_list=False):

        """
        This function ...
        :param origins:
        :param as_list:
        :return:
        """

        # Return as a list?
        if as_list:

            levels = []
            #print(self.sigma_levels_for_images)
            for origin in origins: levels.append(self.sigma_levels_for_images[origin])
            return levels

        # Return as a dictionary
        else:

            levels = dict()
            for origin in origins: levels[origin] = self.sigma_levels_for_images[origin]
            return levels

    # -----------------------------------------------------------------

    @property
    def has_additional_levels(self):

        """
        This function ...
        :return:
        """

        return self.config.additional_levels is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def additional_levels_image_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.config.additional_levels]

    # -----------------------------------------------------------------

    @lazyproperty
    def additional_levels_images(self):

        """
        This function ...
        :return:
        """

        levels = dict()
        for fltr in self.config.additional_levels:
            name = str(fltr)
            levels[name] = self.config.additional_levels[fltr]
        return levels

    # -----------------------------------------------------------------

    def has_additional_levels_for_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_additional_levels and name in self.additional_levels_image_names

    # -----------------------------------------------------------------

    def get_additional_levels_for_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.additional_levels_images[name]

    # -----------------------------------------------------------------

    def get_relative_levels_for_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if not self.has_additional_levels_for_image(name):  relative_levels = self.config.sigma_levels
        else:

            all_levels = self.config.sigma_levels + self.get_additional_levels_for_image(name)
            relative_levels = list(sorted(set(all_levels)))

        # Check whether default is in
        if self.has_default_levels_from_file:

            # Get relative default level
            relative_default_level = self.get_relative_default_sigma_level_for_image_from_file(name)
            relative_levels.append(relative_default_level)

            # Put it in
            relative_levels = list(sorted(set(relative_levels)))

        # Return
        return relative_levels

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

            # Get the factors
            factors = self.get_relative_levels_for_image(name)

            # Get the different sigma levels
            sigma_levels = [level * factor for factor in factors]

            # Set the levels
            levels[name] = sigma_levels

        # Return the levels
        return levels

    # -----------------------------------------------------------------

    @property
    def has_default_levels_from_file(self):

        """
        This function ...
        :return:
        """

        return self.config.default_levels_from is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def default_levels_from_file(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(self.maps_components_path, "levels_" + str(self.config.default_levels_from) + ".dat")

        # Load the levels file
        levels = load_dict(path)

        # Return the levels
        return levels

    # -----------------------------------------------------------------

    @lazyproperty
    def default_levels_from_file_images(self):

        """
        Thisf unction ...
        :return:
        """

        levels = dict()
        for fltr in self.default_levels_from_file:
            name = str(fltr)
            levels[name] = self.default_levels_from_file[fltr]
        return levels

    # -----------------------------------------------------------------

    def get_default_sigma_level_for_image_from_file(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.default_levels_from_file_images[name]

    # -----------------------------------------------------------------

    def get_relative_default_sigma_level_for_image_from_file(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the level
        level = self.significance_levels[name]

        # Get the absolute default level
        default_level = self.get_default_sigma_level_for_image_from_file(name)

        # Return relative
        return default_level / level

    # -----------------------------------------------------------------

    def get_default_sigma_level_for_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Has default levels from file
        if self.has_default_levels_from_file: relative = self.get_relative_default_sigma_level_for_image_from_file(name)

        # Default level from factor
        else: relative = self.config.default_sigma_level

        # Get the level
        level = self.significance_levels[name]

        # Set the default level
        default_level = relative * level

        # Return
        return default_level

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

            # Get the default level
            default_level = self.get_default_sigma_level_for_image(name)

            # Set the level
            levels[name] = default_level

        # Return the dictionary of default sigma levels
        return levels

    # -----------------------------------------------------------------

    def has_all_mask_plots_for_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the plot path for this image
        plot_path = self.image_mask_plot_paths[name]

        # Loop over the levels
        for level in self.sigma_levels_for_images[name]:

            # Determine path
            path = fs.join(plot_path, str(level) + ".png")

            # Determine mask path
            #mask_path = fs.join(plot_path, str(level) + "_mask.png")

            # Check
            # if fs.is_file(path) and fs.is_file(mask_path): pass
            # else: has_all = False
            #if not fs.is_file(path) or not fs.is_file(mask_path): return False
            if not fs.is_file(path): return False

        # Return
        return True

    # -----------------------------------------------------------------

    def plot_image_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the image masks ...")

        # Loop over the images
        for name in self.image_names:

            # Get the sigma levels
            levels = self.sigma_levels_for_images[name]

            # Get plot path
            plot_path = self.image_mask_plot_paths[name]

            # Check whether not all plots are already present
            if self.has_all_mask_plots_for_image(name):
                log.success("All mask plots for the '" + name + "' image are already present")
                for level in levels:
                    path = fs.join(plot_path, str(level) + ".png")
                    self.image_mask_plot_paths_images[name][level] = path
                continue

            # Debugging
            log.debug("Plotting the masks for the '" + name + "' image ...")

            # Get frame
            frame = self.dataset.get_frame(name)

            # Crop the frame
            x_min, x_max, y_min, y_max = frame.crop_to(self.truncation_box, factor=self.config.cropping_factor)

            # Get error map and crop as well
            errormap = self.dataset.get_errormap(name)
            errormap.crop(x_min, x_max, y_min, y_max)

            # Create the significance map
            significance = frame / errormap

            # Create the plots for the different levels
            for level in levels:
                path = self.plot_image_mask_for_level(name, significance, plot_path, level)
                self.image_mask_plot_paths_images[name][level] = path

    # -----------------------------------------------------------------

    def plot_image_mask_for_level(self, name, significance, plot_path, level):

        """
        This function ...
        :param name:
        :param significance:
        :param plot_path:
        :param level:
        :return:
        """

        # Determine path
        path = fs.join(plot_path, str(level) + ".png")

        # Check
        if fs.is_file(path):
            log.success("The image mask plot for the '" + name + "' image at a sigma level of '" + str(level) + "' is already present")
            return path

        # Debugging
        log.debug("Making plot for the '" + name + "' image mask at a sigma level of '" + str(level) + "' ...")

        # Create the mask
        if self.config.fuzzy_mask: mask = self.create_fuzzy_mask_for_level(significance, level, self.config.fuzziness, offset=self.config.fuzzy_min_significance_offset)
        else: mask = self.create_mask_for_level(significance, level)

        # FINISH MASK: keep largest and fill holes
        mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)
        mask.fill_holes()

        # Save the mask
        if self.config.fuzzy_mask: alpha = True
        else: alpha = False

        # Save the mask plot
        mask.saveto_png(path, colour=self.config.mask_colour, alpha=alpha)

        # Return the path
        return path

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
        self.sliders = html.make_multi_image_sliders(names, paths, self.image_names, self.sigma_levels_for_images,
                                                     self.default_sigma_levels_for_images, width=self.image_width,
                                                     height=self.image_height, basic=True, img_class=None,
                                                     label_images=self.image_mask_plot_paths_images, table_class=table_class)

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

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths, javascript_path=javascripts, footing=updated_footing())

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
        with browser.serve_local_host: browser.open_path(self.clip_maps_html_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return None

# -----------------------------------------------------------------
