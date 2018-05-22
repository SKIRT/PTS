#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.view.single View an single image with JS9.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from pts.magic.view.html import body_settings, make_replace_infs_by_nans, make_load_regions, make_replace_negatives_by_nans
from pts.core.tools import html
from pts.core.tools import filesystem as fs
from pts.core.tools import browser
from pts.magic.region.list import load_region_list
from pts.magic.view.html import javascripts, css_scripts, JS9Preloader, JS9Window
from ...core.tools.utils import lazyproperty
from .base import ImageViewer

# -----------------------------------------------------------------

stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"
background_color = "white"
page_style = "ugentstyle"

# -----------------------------------------------------------------

class SingleImageViewer(ImageViewer):

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
        super(SingleImageViewer, self).__init__(*args, **kwargs)

        # The window
        self.window = None

        # The preloader
        self.preloader = None

        # The regions button
        self.regions_button = None

        # The page
        self.page = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. View window
        self.make_view()

        # 3. Regions
        if self.has_regions: self.make_regions()

        # 4. Page
        if self.config.page: self.generate_page()

        # 5. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    @property
    def has_regions(self):

        """
        This function ...
        :return:
        """

        return self.config.regions is not None

    # -----------------------------------------------------------------

    @property
    def filename(self):

        """
        This function ...
        :return:
        """

        return fs.strip_extension(fs.name(self.config.image))

    # -----------------------------------------------------------------

    @property
    def filepath(self):

        """
        This function ...
        :return:
        """

        return self.config.image

    # -----------------------------------------------------------------

    @lazyproperty
    def image_name(self):

        """
        This function ...
        :return:
        """

        return html.make_usable(self.filename)

    # -----------------------------------------------------------------

    @property
    def display_name(self):

        """
        Thisn function ...
        :return:
        """

        return "JS9" + self.image_name

    # -----------------------------------------------------------------

    def make_view(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform theb user
        log.info("Making the view ...")

        # Make settings
        settings = dict()
        settings["scale"] = self.config.scale
        settings["colormap"] = self.config.colormap
        settings["zoom"] = self.config.zoom

        # Create preloader
        self.preloader = JS9Preloader()
        image = self.preloader.add_path(self.filename, self.filepath, settings=settings, display=self.display_name)

        # Create window
        self.window = JS9Window(self.display_name, width=self.config.width, height=self.config.height,
                           background_color=background_color, menubar=self.config.menubar,
                           colorbar=self.config.colorbar, resize=self.config.resize, scrolling=self.config.scrolling,
                           panner=self.config.panner, magnifier=self.config.magnifier,
                           combine_panner_and_magnifier=self.config.combine_panner_and_magnifier)

    # -----------------------------------------------------------------

    @property
    def view(self):

        """
        This function ...
        :return:
        """

        return self.window.view

    # -----------------------------------------------------------------

    def make_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the regions ...")

        # Load the regions
        regions = load_region_list(self.config.regions)

        # Regions button
        region_button_id = self.image_name + "regionsbutton"
        load_region_function_name = "load_regions_" + self.image_name
        load_region = make_load_regions(regions, display=self.display_name, movable=self.config.movable,
                                        rotatable=self.config.rotatable, removable=self.config.removable,
                                        resizable=self.config.resizable, quote_character="'")

        # Create region loader
        self.regions_button = html.make_script_button(region_button_id, "Load regions", load_region,
                                                      load_region_function_name)

    # -----------------------------------------------------------------

    @property
    def infs_button(self):

        """
        This function ...
        :return:
        """

        return self.make_infs_button(self.image_name, self.display_name)

    # -----------------------------------------------------------------

    @property
    def negatives_button(self):

        """
        This function ...
        :return:
        """

        return self.make_negatives_button(self.image_name, self.display_name)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return self.filename

    # -----------------------------------------------------------------

    def _initialize_page(self):

        """
        Thisf unction ...
        :return:
        """

        # Create CSS for the page width
        css = html.make_page_width(self.config.page_width)

        # Create list of css scripts
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = html.HTMLPage(self.title, css=css, body_settings=body_settings, style=page_style,
                                  css_path=css_paths,
                                  javascript_path=javascripts, footing=html.updated_footing())

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

        # Make theme button
        self.page += html.center(self.theme_button) + html.newline

        # Infs
        self.page += html.center(self.infs_button) + html.newline

        # Negatives
        self.page += html.center(self.negatives_button) + html.newline

        # Regions
        if self.has_regions: self.page += html.center(self.regions_button) + html.newline

        # Add the window
        self.page += self.window

        # Add the preloader
        self.page += self.preloader

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
