#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.html Contains the TruncationPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc

# Import the relevant PTS classes and modules
from .component import TruncationComponent
from ...core.tools.html import HTMLPage, SimpleTable, newline, updated_footing, make_theme_button, center, sleep_function, other_sleep_function, make_script_button
from ...core.tools.logging import log
from ...magic.view.html import JS9Viewer, JS9Preloader, body_settings, javascripts, css_scripts, JS9Menubar, JS9Loader, JS9Spawner, JS9Colorbar, JS9Window, make_load_region_function, make_load_region
from ...core.tools import filesystem as fs
from ...core.tools.utils import lazyproperty
from ...core.filter.filter import parse_filter
from .analytics import mask_names
from ...core.tools import browser

# -----------------------------------------------------------------

ncolumns = 3
image_width = 300
image_height = 300

# -----------------------------------------------------------------

base_url = "http://users.ugent.be/~sjversto"
stylesheet_filename = "stylesheet.css"
stylesheet_url = fs.join(base_url, stylesheet_filename)

style = "ugentstyle"

# -----------------------------------------------------------------

class TruncationPageGenerator(TruncationComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(TruncationPageGenerator, self).__init__(*args, **kwargs)

        # --- Attributes ---

        # The preloader
        self.preloader = None

        # The loaders
        self.loaders = dict()

        # The region loaders
        self.region_loaders = dict()

        # The windows
        self.windows = dict()

        # The page
        self.page = None

        # The plots directory
        self.plots_path = None

        # Plot paths
        self.plots_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the masks
        if self.config.mask: self.load_masks()

        # Make plots
        self.make_plots()

        # Make the views
        self.make_views()

        # Make table
        self.make_table()

        # Generate the page
        self.generate_page()

        # 6. Writing
        self.write()

        # 7. Showing
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TruncationPageGenerator, self).setup(**kwargs)

        # Set the preloader
        #if self.config.preload: self.preloader = JS9Preloader()
        self.preloader = JS9Preloader()

        # Make directory to contain the plots
        self.plots_path = fs.join(self.truncation_html_path, "plots")
        if fs.is_directory(self.plots_path): fs.clear_directory(self.plots_path)
        else: fs.create_directory(self.plots_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        # Sorted on wavelength!
        if self.config.filters is not None: return sorted(self.config.filters, key=lambda fltr: fltr.wavelength.to("micron").value)
        #else: return sorted(self.dataset.filters, key=lambda fltr: fltr.wavelength) # TOO SLOW
        else: return sorted((parse_filter(name) for name in self.dataset.names), key=lambda fltr: fltr.wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Truncation"

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        return [self.dataset.get_name_for_filter(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Loop over the images
        #for fltr in self.filters:
        for name, fltr in zip(self.names, self.filters):

            # Get the mask
            mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Add the mask
            self.masks[name] = mask

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make image plots
        if self.config.png: self.make_image_plots()

    # -----------------------------------------------------------------

    def make_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making image plots ...")

        # Loop over the images
        #for name in self.dataset.names:
        for name, fltr in zip(self.names, self.filters):

            # Get name
            #name = self.dataset.get_name_for_filter(fltr)

            # Get frame
            frame = self.dataset.get_frame(name)

            # Determine path
            filepath = fs.join(self.plots_path, name + ".png")

            # Set the path
            self.plots_paths[name] = filepath

            # Make plot
            # Save as PNG image
            frame.saveto_png(filepath, colours=self.config.colormap, absolute_alpha=True)

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def make_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views ...")

        # Get the frames
        #self.frames = self.dataset.get_framelist()

        # Get the error maps
        #self.errormaps = self.dataset.get_errormaplist()

        # Loop over all prepared images, get the images
        #self.masks = dict()
        #for name in self.dataset.names:
        for name, fltr in zip(self.names, self.filters):

            # Get the mask
            #mask_names = ["padded", "bad"]
            #mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Set the mask
            #if mask is None: continue
            #self.masks[name] = mask

            # Get name
            #name = self.dataset.get_name_for_filter(fltr)

            display_name = name.replace(" ", "") # CANNOT CONTAIN SPACES!!

            # Get path
            if self.config.png: path = self.plots_paths[name]
            else: path  = self.dataset.get_frame_path(name)

            settings = dict()
            settings["scale"] = self.config.scale
            settings["colormap"] = self.config.colormap
            #settings["fits2png"] = "true"

            #regions = "ellipse"

            # Get the angle
            center = self.disk_ellipse.center  # in sky coordinates
            semimajor = self.disk_ellipse.semimajor
            semiminor = self.disk_ellipse.semiminor
            angle = self.disk_ellipse.angle

            # Determine the ratio of semimajor and semiminor
            #ratio = semiminor / semimajor

            regions_for_loader = region if self.config.load_regions else None

            # Add preload
            if self.config.preload_all or (self.config.preload is not None and fltr in self.config.preload):

                # Add to preloader
                self.preloader.add_path(name, path, settings=settings, display=display_name, regions=regions_for_loader)

                # Create window
                self.windows[name] = JS9Window(display_name, width=image_width, height=image_height, background_color="white", menubar=self.config.menubar, colorbar=self.config.colorbar, resize=self.config.resize)

                display_id = display_name

            # Add dynamic load
            elif self.config.dynamic:

                self.loaders[name] = JS9Spawner.from_path("Load image", name, path, settings=settings, button=True,
                                                          menubar=self.config.menubar, colorbar=self.config.colorbar, regions=regions_for_loader, add_placeholder=False)
                display_id = self.loaders[name].display_id

                self.windows[name] = self.loaders[name].placeholder

            # Regular button load in a pre-existing viewer
            else:

                # Create loader
                self.loaders[name] = JS9Loader.from_path("Load image", name, path, display=display_name,
                                                         settings=settings, button=True, regions=regions_for_loader)

                # Create window
                self.windows[name] = JS9Window(display_name, width=image_width, height=image_height, background_color="white", menubar=self.config.menubar, colorbar=self.config.colorbar, resize=self.config.resize)

                display_id = display_name

            # Regions button
            region_button_id = display_name + "regionsbutton"
            load_region_function_name = "load_regions_" + display_name
            # load_region = make_load_region_function(load_region_function_name, regions, display=None)
            load_region = make_load_region(regions, display=display_id)

            # Create region loader
            self.region_loaders[name] = make_script_button(region_button_id, "Load regions", load_region,
                                                           load_region_function_name)

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        if self.config.dynamic: return 1
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
        #for name in self.dataset:
        #for fltr in self.filters:
        for name, fltr in zip(self.names, self.filters):

            # Get name
            #name = self.dataset.get_name_for_filter(fltr)

            # Add name of the image
            string = center(name + newline)

            # Add loader if necessary
            if name in self.loaders:
                string += center(str(self.loaders[name]))
                string += newline

            # Add region loader if necessary
            if name in self.region_loaders:
                string += center(str(self.region_loaders[name]))
                string += newline

            # Add the window
            if name in self.windows: string += str(self.windows[name])

            # Add the cell
            cells.append(string)

        # Make the table
        self.table = SimpleTable.rasterize(cells, ncolumns=self.ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        #settings_body = body_settings if self.config.preload else None
        #settings_body = body_settings
        settings_body = body_settings if self.preloader.has_images else None

        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = HTMLPage(self.title, body_settings=settings_body, javascript_path=javascripts, css_path=css_paths, style=style, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"

        self.page += "<script>" + other_sleep_function + "</script>"

        self.page += center(make_theme_button(classes=classes))

        # Add the table
        self.page += self.table

        # Add preloader
        #if self.config.preload: self.page += self.preloader
        if self.preloader.has_images: self.page += self.preloader

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
        self.page.saveto(self.truncation_html_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.truncation_html_page_path)

# -----------------------------------------------------------------
