#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.status Contains the StatusPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class, hover_table_class
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

class StatusPageGenerator(HTMLPageComponent):

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
        super(StatusPageGenerator, self).__init__(*args, **kwargs)

        # Tables
        self.status_table = None
        self.history_table = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make tables
        self.make_tables()

        # Make plots
        self.make_plots()

        # Generate the html
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
        super(StatusPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make status table
        self.make_status_table()

        # Make the history table
        self.make_history_table()

    # -----------------------------------------------------------------

    def make_status_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the status table ...")

        # Set the background colours of the cells
        bgcolors = [(None, color) for color in self.status.colors]

        # Create the table
        self.status_table = html.SimpleTable(self.status, header=["Step", "Status"], bgcolors=bgcolors, css_class=hover_table_class)

    # -----------------------------------------------------------------

    def make_history_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the history table ...")

        # Create the table
        self.history_table = html.SimpleTable.from_table(self.history, css_class=hover_table_class)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the status
        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the status page ...")

        # Create list of css scripts
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Make javascripts urls
        javascript_paths = javascripts[:]
        javascript_paths.append(sortable_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        #classes = dict()
        #classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        #body = self.heading

        # Create titles
        title_status = html.underline_template.format(text="Modeling status")
        title_history = html.underline_template.format(text="Modeling history")

        self.page += html.line + html.newline
        self.page += title_status + html.newline + html.newline + str(self.status_table) + html.newline + html.newline
        self.page += html.line + html.newline
        self.page += title_history + html.newline + html.newline + str(self.history_table) + html.newline + html.newline

        #body += self.footing

        # Create the status page
        #self.make_page(body)

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

        return self.status_page_path

# -----------------------------------------------------------------
