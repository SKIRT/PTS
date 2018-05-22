#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.photometry Contains the PhotometryPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, sortable_table_class
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ...core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

page_width = 600
thumbnail_title = "Thumbnail"

# -----------------------------------------------------------------

bokeh_stylesheet_url = "https://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.css"
bokeh_javascript_url = "https://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.js"

# HAS TO BE ADDED FOR BOKEH:
# <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.css" type="text/css" />
# <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.js"></script>

# -----------------------------------------------------------------

class PhotometryPageGenerator(HTMLPageComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PhotometryPageGenerator, self).__init__(*args, **kwargs)

        # The sed plot
        self.sed_plot = None
        self.sed_plot_script = None

        # The sed table
        self.sed_table = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make plots
        self.make_plots()

        # 3. Make tables
        self.make_tables()

        # 4. Generate the html
        self.generate()

        # 5. Write
        self.write()

        # 6. Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PhotometryPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # SED plot
        self.make_sed_plot()

        # Image plot
        self.make_image_plots()

    # -----------------------------------------------------------------

    def make_sed_plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making SED plot ...")

        # Create the SED plotter
        plotter = SEDPlotter()

        # Use bokeh
        plotter.config.library = "bokeh"

        # Don't show, but get the figure
        plotter.config.show = False

        # Add the SED
        plotter.add_sed(self.observed_sed, "Observation")
        plotter.add_sed(self.asymptotic_sed_path, "Asymptotic")
        plotter.add_sed(self.truncated_sed_path, "Truncated")

        # Run the plotter
        plotter.run()

        # Get the figure
        self.sed_plot_script, self.sed_plot = plotter.figure.to_html_components(wrap_script=False)

    # -----------------------------------------------------------------

    def make_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making image plots ...")

        # Loop over the images
        for name in self.photometry_image_names:

            # Make plots

            # Eventually make (JS9) views

            pass

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        Thisf ucntion ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Create the table
        self.sed_table = SimpleTable.from_table(self.observed_sed, css_class=sortable_table_class)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the page
        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Create list of css scripts
        #css_paths = css_scripts[:]

        css_paths = []
        css_paths.append(stylesheet_url)
        css_paths.append(bokeh_stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Make javascripts urls
        #javascript_paths = javascripts[:]

        javascript_paths = []
        javascript_paths.append(sortable_url)
        #javascript_paths.append(preview_url)
        javascript_paths.append(bokeh_javascript_url)

        javascript = self.sed_plot_script

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing(), javascript_header=javascript)

        #classes = dict()
        #classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        self.page += "Click on any column to sort"
        self.page += html.newline

        # Add the table
        self.page += self.sed_table
        self.page += html.newline + html.newline

        # Add the plot
        self.page += self.sed_plot
        self.page += html.newline

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write status page
        self.write_page()

    # -----------------------------------------------------------------

    @property
    def page_path(self):

        """
        This function ...
        :return:
        """

        return self.photometry_page_path

# -----------------------------------------------------------------
