#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.html.significance Contains the SignificanceLevelsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import TruncationComponent
from ...html.component import stylesheet_url, page_style
from ...html.component import slider_stylesheet_url, slider_url
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts
from ....core.tools import browser
from ....core.tools.stringify import tostr
from ....core.tools.utils import lazyproperty
from ....magic.core.rgb import RGBImage
from ....core.tools import filesystem as fs
from ....magic.core.mask import Mask
from ....core.tools import sequences
from ....core.filter.filter import parse_filter

# -----------------------------------------------------------------

significance_plots_name = "significance_plots"
ncolumns = 2
colour_map = "jet"

page_width = 600

# -----------------------------------------------------------------

class SignificanceLevelsPageGenerator(TruncationComponent):

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
        super(SignificanceLevelsPageGenerator, self).__init__(*args, **kwargs)

        # Plot paths for each filter
        self.filter_plot_paths = dict()

        # Paths of the images
        self.level_plot_paths = defaultdict(OrderedDict)

        # The sliders
        self.sliders = dict()

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

        # Make the sliders
        self.make_sliders()

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
        super(SignificanceLevelsPageGenerator, self).setup(**kwargs)

        # Make directory to contain the plots
        self.plots_path = fs.join(self.truncation_html_path, "significance")
        if fs.is_directory(self.plots_path):
            if self.config.replot: fs.clear_directory(self.plots_path)
        else: fs.create_directory(self.plots_path)

        # Create directories for each filter
        for fltr in self.filters:
            fltr_path = fs.join(self.plots_path, str(fltr))
            self.filter_plot_paths[fltr] = fltr_path
            if not fs.is_directory(fltr_path): fs.create_directory(fltr_path)

        # Set the default sigma level
        self.config.default_level = sequences.find_closest_value(self.config.sigma_levels, self.config.default_level)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        return self.dataset.filters

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        return self.dataset.names

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Significance maps"

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

    def has_all_plots(self, name):

        """
        This function ...
        :return:
        """

        # Get plot path
        plot_path = self.filter_plot_paths[parse_filter(name)]

        #print(fs.files_in_path(plot_path, returns="name", extension="png", convert=float), self.config.sigma_levels)

        # Loop over the levels
        for level in self.config.sigma_levels:

            # Determine path
            path = fs.join(plot_path, str(level) + ".png")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def check_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the plot path
        plot_path = self.filter_plot_paths[parse_filter(name)]

        has_all = True

        # Loop over the levels
        for level in self.config.sigma_levels:

            # Determine path
            path = fs.join(plot_path, str(level) + ".png")

            # Check
            if fs.is_file(path): self.level_plot_paths[name][level] = path
            else: has_all = False

        # Return
        return has_all

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Loop over the frames
        for name in self.names:

            # Check whether not all plots are already present
            #if self.has_all_plots(name): continue
            if self.check_plots(name): continue

            # Get the filter
            fltr = self.dataset.get_filter(name)

            # Get plot path
            plot_path = self.filter_plot_paths[fltr]

            # Get frame
            #frame = self.dataset.get_frame_for_filter(fltr)
            frame = self.dataset.get_frame(name)

            # Get error map list
            #errormap = self.dataset.get_errormap_for_filter(fltr)
            errormap = self.dataset.get_errormap(name)

            # Create the significance map
            significance = frame / errormap

            # Create the plots
            for level in self.config.sigma_levels:

                # Determine path
                path = fs.join(plot_path, str(level) + ".png")

                # Check
                if fs.is_file(path): continue

                # Create the mask
                mask = Mask(significance > level)

                # Fill holes
                mask.fill_holes()

                # Invert
                #mask.invert()

                # Create RGB image
                image = RGBImage.from_mask(mask)

                # Save the image
                image.saveto(path)

                # Set the path
                self.level_plot_paths[name][level] = path

    # -----------------------------------------------------------------

    def make_sliders(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the image sliders ...")

        # Loop over the images
        for name in self.names:

            # Get the labels and the urls
            labels = self.level_plot_paths[name].keys()
            paths = self.level_plot_paths[name].values()

            # Create image ID
            image_id = name.replace("_", "").replace(" ", "")

            # Create the slider
            slider = html.make_image_slider(image_id, paths, labels, self.config.default_level, width=self.image_width, height=self.image_height)

            # Set the slider
            self.sliders[name] = slider

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Create list of css scripts
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)
        css_paths.append(slider_stylesheet_url)

        # Create CSS for the page width
        css = html.make_page_width(page_width)

        # Make javascripts urls
        javascript_paths = javascripts[:]
        #javascript_paths.append(sortable_url)
        #javascript_paths.append(preview_url)
        javascript_paths.append(slider_url)

        # Create the page
        self.page = HTMLPage(self.title, style=page_style, css_path=css_paths,
                             javascript_path=javascripts, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes))

        self.page += html.newline

        self.page += html.line
        self.page += html.newline

        # Add the sliders
        for name in self.names:

            # Add the name
            self.page += name.upper()
            self.page += html.newline

            # Add the image slider
            self.page += self.sliders[name]

            # Add line
            self.page += html.line
            self.page += html.newline

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
        self.page.saveto(self.significance_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.significance_page_path)

# -----------------------------------------------------------------
