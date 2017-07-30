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
from ...core.tools.html import HTMLPage, SimpleTable, newline, updated_footing, make_theme_button, center
from ...core.tools.logging import log
from ...magic.view.html import JS9Viewer, JS9Preloader, body_settings, javascripts, css_scripts, JS9Menubar, JS9Loader, JS9Spawner
from ...core.tools import filesystem as fs
from ...core.tools.utils import lazyproperty

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

        # The menubars
        self.menus = dict()

        # The views
        self.views = dict()

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
        if self.config.filters is not None: return sorted(self.config.filters, key=lambda fltr: fltr.wavelength)
        else: sorted(self.dataset.filters, key=lambda fltr: fltr.wavelength)

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
        for fltr in self.filters:

            # Get name
            name = self.dataset.get_name_for_filter(fltr)

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
        for fltr in self.filters:

            # Get the mask
            #mask_names = ["padded", "bad"]
            #mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Set the mask
            #if mask is None: continue
            #self.masks[name] = mask

            # Get name
            name = self.dataset.get_name_for_filter(fltr)

            display_name = name.replace(" ", "") # CANNOT CONTAIN SPACES!!

            # Get path
            if self.config.png: path = self.plots_paths[name]
            else: path  = self.dataset.get_frame_path(name)

            settings = dict()
            settings["scale"] = self.config.scale
            settings["colormap"] = self.config.colormap
            #settings["fits2png"] = "true"

            # Set menu name
            # menu_name = display_name + "_menu"
            menu_name = display_name + "Menubar"
            # displays = [display_name]
            displays = None

            # Add preload
            if self.config.preload_all or (self.config.preload is not None and fltr in self.config.preload):

                # Add to preloader
                self.preloader.add_path(name, path, settings=settings, display=display_name)

                # Create view
                self.views[name] = JS9Viewer(display_name, width=image_width, height=image_height)

                # Add menu bar
                if self.config.menubar: self.menus[name] = JS9Menubar(menu_name, displays=displays, width=image_width, background_color="white")

            # Add dynamic load
            elif self.config.dynamic: self.loaders[name] = JS9Spawner.from_path("Load image", name, path, settings=settings, button=True, menubar=self.config.menubar)

            # Regular button load in a pre-existing viewer
            else:

                # Create loader
                self.loaders[name] = JS9Loader.from_path("Load image", name, path, display=display_name, settings=settings, button=True)

                # Create view
                self.views[name] = JS9Viewer(display_name, width=image_width, height=image_height)

                # Add menu bar
                if self.config.menubar: self.menus[name] = JS9Menubar(menu_name, displays=displays, width=image_width, background_color="white")

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
        for fltr in self.filters:

            # Get name
            name = self.dataset.get_name_for_filter(fltr)

            # Add name of the image
            string = center(name + newline)

            # Add loader if necessary
            if name in self.loaders: string += center(str(self.loaders[name]))

            string += newline

            # Add menu bar (if requested)
            if name in self.menus: string += str(self.menus[name])

            string += newline

            # Add view (if not dynamic)
            if name in self.views: string += str(self.views[name])

            # Add the cell
            cells.append(string)

        # Make the table
        self.table = SimpleTable.rasterize(cells, ncolumns=ncolumns)

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
        fs.open_in_browser(self.truncation_html_page_path)

# -----------------------------------------------------------------
