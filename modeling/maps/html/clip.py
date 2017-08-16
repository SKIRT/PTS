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

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import MapsComponent
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts
from ....core.tools import browser
from ....core.tools.stringify import tostr
from ....core.tools.utils import lazyproperty
from ....magic.core.rgb import RGBImage
from ....core.tools import filesystem as fs
from ....magic.core.mask import Mask

# -----------------------------------------------------------------

clipped_name = "clipped"
ncolumns = 2
colour_map = "jet"

# -----------------------------------------------------------------

class ClipMapsPageGenerator(MapsComponent):

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

        # Plot paths for each filter
        self.filter_plot_paths = dict()

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

        # Generate the page
        self.generate_page()

        # 5. Writing
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

        # Create directories for each filter
        for fltr in self.maps_filters:
            fltr_path = fs.join(self.clipped_plots_path, str(fltr))
            self.filter_plot_paths[fltr] = fltr_path
            if not fs.is_directory(fltr_path): fs.create_directory(fltr_path)

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

        #return 150
        return None

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        This function ...
        :return:
        """

        return 300

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters()

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Loop over the filters
        for fltr in self.maps_filters:

            # Get plot path
            plot_path = self.filter_plot_paths[fltr]

            # Get frame
            frame = self.dataset.get_frame_for_filter(fltr)

            # Get error map list
            errormap = self.dataset.get_errormap_for_filter(fltr)

            # Create the significance map
            significance = frame / errormap

            # Create the plots
            for level in self.config.sigma_levels:

                # Create the mask
                mask = Mask(significance > level)

                # Fill holes
                mask.fill_holes()

                # Invert
                #mask.invert()

                # Create RGB image
                image = RGBImage.from_mask(mask)

                # Determine path
                path = fs.join(plot_path, str(level) + ".png")

                # Save the image
                image.saveto(path)

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



    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return None

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
