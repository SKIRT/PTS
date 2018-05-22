#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.view.multi View multiple images with JS9.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .base import ImageViewer
from ...core.tools.html import HTMLPage, SimpleTable, newline, updated_footing, center, make_script_button, dictionary, make_usable
from ...core.basics.log import log
from .html import body_settings, javascripts, css_scripts, JS9Loader, JS9Spawner, JS9Window, make_load_regions, JS9Preloader, load_region_or_regions
from .html import make_replace_infs_by_nans_multiple, make_replace_negatives_by_nans_multiple
from .html import make_spawn_code, add_to_div
from ...core.tools import filesystem as fs
from ..tools.info import get_image_info_from_header_file
from ...core.tools import browser
from ..region.list import load_as_pixel_region_list, RegionList, PixelRegionList, SkyRegionList
from ..region.region import PixelRegion, SkyRegion, Region
from ...core.tools import types
from ..basics.coordinatesystem import CoordinateSystem
from ...core.tools import strings

# -----------------------------------------------------------------

ncolumns = 3
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

base_url = "http://users.ugent.be/~sjversto"
stylesheet_filename = "stylesheet.css"
stylesheet_url = fs.join(base_url, stylesheet_filename)

# -----------------------------------------------------------------

style = "ugentstyle"

# -----------------------------------------------------------------

class MultiImageViewer(ImageViewer):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MultiImageViewer, self).__init__(*args, **kwargs)

        # The image paths
        self.paths = None

        # The regions
        self.regions = None

        # The preloader
        self.preloader = None

        # The loaders
        self.loaders = dict()

        # The region loaders
        self.region_loaders = dict()

        # The windows
        self.windows = dict()

        # All images loader
        self.all_loader = None

        # The table
        self.table = None

        # Info
        self.info = dict()

        # The display IDs
        self.display_ids = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make info
        if self.config.info: self.get_info()

        # 3. Make the views
        self.make_views()

        # 4. Make table
        self.make_table()

        # 5. Generate the page
        if self.config.page: self.generate_page()

        # 6. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(MultiImageViewer, self).setup(**kwargs)

        # Get paths
        if kwargs.get("paths", None) is not None: self.paths = kwargs.pop("paths")
        elif self.config.image_paths is not None: self.set_paths(self.config.image_paths)
        else: self.load_paths()

        # Check: are there images?
        if not self.has_names: raise ValueError("There are no images")

        # Get regions
        if kwargs.get("regions", None) is not None:

            regions = kwargs.pop("regions")
            if isinstance(regions, RegionList): self.set_regions(regions)
            elif types.is_dictionary(regions): self.regions = regions
            else: raise ValueError("Invalid type for 'regions': must be region list or dictionary of region lists")

        # Single region is passed
        elif kwargs.get("region", None) is not None: self.set_region(kwargs.pop("region"))

        # Load regions from file
        else: self.load_regions()

        # Initialize preloader if necessary
        if self.preload_any: self.preloader = JS9Preloader()

    # -----------------------------------------------------------------

    def set_paths(self, image_paths):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Setting the image paths ...")

        # Initialize the dictionary
        self.paths = dict()

        # Loop over the paths
        for path in image_paths:

            # Get name
            name = fs.strip_extension(fs.name(path))

            # Debugging
            log.debug("Processing '" + name + "' image ...")

            # Set
            self.paths[name] = path

    # -----------------------------------------------------------------

    def load_paths(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Loading the image paths ...")

        # Get the files
        self.paths = fs.files_in_cwd(extension=self.config.extensions, recursive=self.config.recursive,
                                          contains=self.config.contains, not_contains=self.config.not_contains,
                                          exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                          returns="dict")

    # -----------------------------------------------------------------

    def set_regions(self, region_list):

        """
        This function ...
        :param region_list:
        :return:
        """

        # Inform the user
        log.info("Setting the regions ...")

        # Initialize
        self.regions = dict()

        # Loop over the images
        for name in self.names:

            # Debugging
            log.debug("Processing '" + name + "' image regions ...")

            # Check type
            if isinstance(region_list, PixelRegionList): regions = region_list
            elif isinstance(region_list, SkyRegionList):
                wcs = self.get_coordinate_system(name)
                regions = region_list.to_pixel(wcs)
            else: raise ValueError("Invalid argument")

            # Set the regions
            self.regions[name] = regions

    # -----------------------------------------------------------------

    def set_region(self, region):

        """
        This function ...
        :param region:
        :return:
        """

        # Inform the user
        log.info("Setting the region ...")

        # Initialize
        self.regions = dict()

        # Loop over the images
        for name in self.names:

            # Debugging
            log.debug("Processing '" + name + "' image region ...")

            # Check type
            if isinstance(region, PixelRegion): pixel_region = region
            elif isinstance(region, SkyRegion):
                wcs = self.get_coordinate_system(name)
                pixel_region = region.to_pixel(wcs)
            else: raise ValueError("Invalid argument")

            # Create pixel region list
            regions = PixelRegionList.single(pixel_region)

            # Set
            self.regions[name] = regions

    # -----------------------------------------------------------------

    def load_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions ...")

        # Initialize
        self.regions = dict()

        # Loop over the images
        for name in self.names:

            # Get the path
            path = self.paths[name]
            directory_path = fs.directory_of(path)

            # Regions filename is defined
            if self.config.regions_name is not None:

                # Determine path
                filepath = fs.absolute_or_in(self.config.regions_name, directory_path)

                # Check whether extension is present
                if not fs.has_extension(filepath): filepath += "." + self.config.regions_extension

            # No filename is determined: assume name of the image
            else: filepath = fs.join(directory_path, self.config.regions_prefix + name + self.config.regions_suffix + "." + self.config.regions_extension)

            # Check if file exists
            if not fs.is_file(filepath): continue

            # Load the regions
            regions = load_as_pixel_region_list(filepath, path)

            # Set
            self.regions[name] = regions

        # Reset to None
        if len(self.regions) == 0: self.regions = None

    # -----------------------------------------------------------------

    def get_coordinate_system(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return CoordinateSystem.from_file(self.paths[name])

    # -----------------------------------------------------------------

    @property
    def has_regions(self):

        """
        Thisf unction ...
        :return:
        """

        return self.regions is not None

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.paths.keys()

    # -----------------------------------------------------------------

    @property
    def nnames(self):

        """
        This function ...
        :return:
        """

        return len(self.names)

    # -----------------------------------------------------------------

    @property
    def has_names(self):

        """
        This function ...
        :return:
        """

        return self.nnames > 0

    # -----------------------------------------------------------------

    def get_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info for the images ...")

        # Loop over the images
        for name in self.names:

            # Get the path
            path = self.paths[name]

            # Get the image info
            info = get_image_info_from_header_file(name, path, path=False, name=False)

            # Make list
            code = dictionary(info, key_color=key_color)

            # Add
            self.info[name] = code

    # -----------------------------------------------------------------

    @property
    def base_zoom(self):

        """
        This function ...
        :return:
        """

        return self.config.zoom.split(";")[0]

    # -----------------------------------------------------------------

    @property
    def next_zoom(self):

        """
        This function ...
        :return:
        """

        if ";" in self.config.zoom: return self.config.zoom.split(";")[1]
        else: return None

    # -----------------------------------------------------------------

    @property
    def preload_any(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.preload_all: return True
        else:
            for name in self.names:
                if self.config.preload is not None and name in self.config.preload: return True
            return False

    # -----------------------------------------------------------------

    def make_views(self):

        """
        Thisj function ...
        :return:
        """

        load_info = dict()
        images = dict()
        region_loads = dict()
        placeholders = dict()
        views = dict()

        # Loop over all images
        for name in self.names:

            # Get name
            display_name = make_usable(name)

            # Get path
            path = self.paths[name]

            settings = dict()
            settings["scale"] = self.config.scale
            settings["colormap"] = self.config.colormap
            settings["zoom"] = self.base_zoom
            # settings["fits2png"] = "true"

            # Get regions
            region = self.regions[name] if self.has_regions and name in self.regions else None
            regions_for_loader = region if self.config.load_regions else None

            # Add preload
            if self.config.preload_all or (self.config.preload is not None and name in self.config.preload):

                # Add to preloader
                image = self.preloader.add_path(name, path, settings=settings, display=display_name,
                                                regions=regions_for_loader, zoom=self.next_zoom)

                # Create window
                self.windows[name] = JS9Window(display_name, width=self.config.width, height=self.config.height,
                                               background_color=background_color, menubar=self.config.menubar,
                                               colorbar=self.config.colorbar, resize=self.config.resize)
                display_id = display_name

                # Set load info
                views[display_id] = self.windows[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = image

            # Add dynamic load
            elif self.config.dynamic:

                # Create the loader
                self.loaders[name] = JS9Spawner.from_path("Load image", name, path, settings=settings, button=True,
                                                          menubar=self.config.menubar, colorbar=self.config.colorbar,
                                                          regions=regions_for_loader, add_placeholder=False,
                                                          background_color=background_color, zoom=self.next_zoom)
                display_id = self.loaders[name].display_id

                self.windows[name] = self.loaders[name].placeholder

                # Set load info
                views[display_id] = self.loaders[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = self.loaders[name].image
                placeholders[display_id] = self.loaders[name].spawn_div_name

            # Regular button load in a pre-existing viewer
            else:

                # Create loader
                self.loaders[name] = JS9Loader.from_path("Load image", name, path, display=display_name,
                                                         settings=settings, button=True, regions=regions_for_loader,
                                                         zoom=self.next_zoom)

                # Create window
                self.windows[name] = JS9Window(display_name, width=self.config.width, height=self.config.height,
                                               background_color=background_color, menubar=self.config.menubar,
                                               colorbar=self.config.colorbar, resize=self.config.resize)

                display_id = display_name

                # Set load info
                views[display_id] = self.windows[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = self.loaders[name].image

            # Set display ID
            self.display_ids[name] = display_id

            # Add region
            if region is not None:

                # Regions button
                region_button_id = display_name + "regionsbutton"
                load_region_function_name = "load_regions_" + display_name

                # Make load region code
                load_region = load_region_or_regions(region, display=display_id, changeable=self.config.regions.changeable,
                                                     movable=self.config.regions.movable, rotatable=self.config.regions.rotatable,
                                                     removable=self.config.regions.removable, resizable=self.config.resizable,
                                                     color=self.config.regions.color, quote_character="'")
                region_loads[display_id] = load_region

                # Create region loader
                self.region_loaders[name] = make_script_button(region_button_id, "Load regions", load_region,
                                                               load_region_function_name)


        all_loader_name = "allimagesloaderbutton"
        all_loader_text = "Load all images"

        load_script = ""

        # Load over the images (displays)
        for display_id in load_info:

            # name, path, regions_for_loader = load_info[display_id]
            image = images[display_id]
            view = views[display_id]

            load_image = image.load()
            if display_id in region_loads: load_region = region_loads[display_id]
            else: load_region = None

            if display_id in placeholders:

                spawn_div_name = placeholders[display_id]

                # Make spawn code
                spawn_code = make_spawn_code(view, image, menubar=self.config.menubar, colorbar=self.config.colorbar,
                                             width=self.config.width, background_color=background_color)

                # Add html code to DIV
                load_script += add_to_div(spawn_div_name, spawn_code)

                # Add DIV to JS9
                load_script += "JS9.AddDivs('" + display_id + "');\n"

            # Load image and region code
            load_script += load_image
            load_script += "\n"
            if load_region is not None: load_script += load_region + "\n"
            load_script += "\n"

        # Make all images loader
        function_name = "loadAllImages"
        self.all_loader = make_script_button(all_loader_name, all_loader_text, load_script, function_name)

    # -----------------------------------------------------------------

    @property
    def display_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in self.names: names.append(self.display_ids[name])
        return names

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        if self.config.dynamic:
            if self.config.info: return 2
            else: return 1
        else:
            if self.config.info: return 2
            else: return ncolumns

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "realtable"

    # -----------------------------------------------------------------

    def make_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making table ...")

        cells = []

        # Loop over the images
        for name in self.names:

            string = ""

            # Split the name of the image over multiple lines
            if len(name) > self.config.max_ncharacters_title:

                lines = strings.split_in_lines(name, length=self.config.max_ncharacters_title, as_list=True)
                for line in lines: string += center(line)
                string += newline

            # Add name of the image
            else: string += center(name) + newline

            # Add loader if necessary
            if name in self.loaders:
                string += center(str(self.loaders[name]))

            # Add region loader if necessary
            if name in self.region_loaders:
                string += center(str(self.region_loaders[name]))

            # Get display name
            display_name = self.display_ids[name]

            # Infs
            string += center(str(self.make_infs_button(name, display_name)))

            # Ngeatives
            string += center(str(self.make_negatives_button(name, display_name)))

            # Add the window
            if name in self.windows: string += str(self.windows[name])

            # Add the cell
            cells.append(string)

            # Add info
            if self.config.info: cells.append(self.info[name])

        # Make the table
        self.table = SimpleTable.rasterize(cells, ncolumns=self.ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "test"

    # -----------------------------------------------------------------

    @property
    def all_infs_button(self):

        """
        This function ...
        :return:
        """

        # Create nan/infs replacer button
        button_id = "allnansinfs"
        replace_function_name = "replace_infs_nans_all"
        replace_nans_infs = make_replace_infs_by_nans_multiple(self.display_names)

        # Create the button
        return make_script_button(button_id, "Replace infs", replace_nans_infs, replace_function_name)

    # -----------------------------------------------------------------

    @property
    def all_negatives_button(self):

        """
        Thisf unction ...
        :return:
        """

        # Create nan/infs replacer button
        button_id = "allnegatives"
        replace_function_name = "replace_negatives_nans_all"
        replace_nans_negatives = make_replace_negatives_by_nans_multiple(self.display_names)

        # Create the button
        return make_script_button(button_id, "Replace negatives", replace_nans_negatives, replace_function_name)

    # -----------------------------------------------------------------

    def _initialize_page(self):

        """
        This function ...
        :return:
        """

        settings_body = body_settings if self.preloader is not None and self.preloader.has_images else None

        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = HTMLPage(self.title, body_settings=settings_body, javascript_path=javascripts,
                             css_path=css_paths, style=style, footing=updated_footing())

    # -----------------------------------------------------------------

    def _add_theme_button(self):

        """
        This function ...
        :return:
        """

        # Make theme button
        self.page += center(self.theme_button) + newline

    # -----------------------------------------------------------------

    def _add_infs_button(self):

        """
        This function ...
        :return:
        """

        # Infs
        self.page += center(self.all_infs_button) + newline

    # -----------------------------------------------------------------

    def _add_negatives_button(self):

        """
        Thisf unction ...
        :return:
        """

        # Ngeatives
        self.page += center(self.all_negatives_button) + newline

    # -----------------------------------------------------------------

    def _add_all_loader(self):

        """
        Thisfunction ...
        :return:
        """

        # Add all images loader
        if self.all_loader is not None:

            self.page += center(str(self.all_loader))
            self.page += newline

    # -----------------------------------------------------------------

    def _add_preloader(self):

        """
        This function ...
        :return:
        """

        # Add preloader
        if self.preloader is not None and self.preloader.has_images: self.page += self.preloader

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Initialize
        self._initialize_page()

        # Add theme button
        self._add_theme_button()

        # Infs and negatives
        self._add_infs_button()
        self._add_negatives_button()

        # All images loader
        self._add_all_loader()

        # Add the table
        self.page += self.table

        # Add preloader
        self._add_preloader()

    # -----------------------------------------------------------------

    def show(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Open the page
        browser.open_page(self.page)

# -----------------------------------------------------------------
