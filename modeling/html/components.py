#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.components Contains the ComponentsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, hover_table_class
from ...core.tools import html
from ..plotting.model import plot_galaxy_components, generate_html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

widget_width = 400
widget_height = 500
widget_style = "minimal"

# -----------------------------------------------------------------

class ComponentsPageGenerator(HTMLPageComponent):

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
        super(ComponentsPageGenerator, self).__init__(*args, **kwargs)

        # Plots
        self.bulge_plot = None
        self.disk_plot = None
        self.model_plot = None

        # Tables
        self.bulge_table = None
        self.disk_table = None
        self.model_table = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make plots
        self.make_plots()

        # Make tables
        self.make_tables()

        # Generaet the html
        self.generate()

        # Write
        self.write()

        # Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentsPageGenerator, self).setup(**kwargs)

        # Set the scripts path
        self.scripts_path = fs.create_directory_in(self.html_path, "scripts_components")

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Plot the bulge
        self.plot_bulge()

        # Plot the disk
        self.plot_disk()

        # Plot the model
        self.plot_model()

    # -----------------------------------------------------------------

    @lazyproperty
    def make_plot_data_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["width"] = widget_width
        kwargs["height"] = widget_height
        kwargs["style"] = widget_style
        return kwargs

    # -----------------------------------------------------------------

    @lazyproperty
    def render_kwargs(self):

        """
        Thisf unction ...
        :return:
        """

        kwargs = dict()
        kwargs["only_body"] = True
        return kwargs

    # -----------------------------------------------------------------

    def plot_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the bulge component ...")

        # Plot, create HTML
        components = {"bulge": self.bulge_model}
        box = plot_galaxy_components(components, draw=True, show=False, **self.make_plot_data_kwargs)
        title = "Bulge model"
        self.bulge_plot = generate_html(box, title, self.scripts_path, **self.render_kwargs)

    # -----------------------------------------------------------------

    def plot_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the disk component ...")

        # Plot, create HTML
        components = {"disk": self.disk_model}
        box = plot_galaxy_components(components, draw=True, show=False, **self.make_plot_data_kwargs)
        title = "Disk model"
        self.disk_plot = generate_html(box, title, self.scripts_path, **self.render_kwargs)

    # -----------------------------------------------------------------

    def plot_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model (disk + bulge) ...")

        # Plot, create HTML
        components = {"disk": self.disk_model, "bulge": self.bulge_model}
        box = plot_galaxy_components(components, draw=True, show=False, **self.make_plot_data_kwargs)
        title = "Model"
        self.model_plot = generate_html(box, title, self.scripts_path, **self.render_kwargs)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make bulge table
        self.make_bulge_table()

        # Make disk table
        self.make_disk_table()

    # -----------------------------------------------------------------

    def make_bulge_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the bulge component table ...")

        # Make the table
        self.bulge_table = html.SimpleTable.from_composite(self.bulge_model, css_class=hover_table_class)

    # -----------------------------------------------------------------

    def make_disk_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the disk component table ...")

        # Make the table
        self.disk_table = html.SimpleTable.from_composite(self.disk_model, css_class=hover_table_class)

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

        # Create list of css scripts
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Make javascripts urls
        javascript_paths = javascripts[:]
        javascript_paths.append(sortable_url)
        javascript_paths.append(preview_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes, images=False))
        self.page += html.newline

        self.page += "BULGE"
        self.page += html.newline + html.newline
        self.page += self.bulge_table
        self.page += html.newline + html.newline
        self.page += self.bulge_plot
        self.page += html.newline + html.newline

        self.page += "DISK"
        self.page += html.newline + html.newline
        self.page += self.disk_table
        self.page += html.newline + html.newline
        self.page += self.disk_plot
        self.page += html.newline + html.newline

        self.page += "MODEL"
        self.page += html.newline + html.newline
        self.page += self.model_plot
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

        return self.components_page_path

# -----------------------------------------------------------------
