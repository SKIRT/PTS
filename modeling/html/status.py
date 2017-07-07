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
from ...core.tools.logging import log
from .component import HTMLPageComponent, stylesheet_url
from ...core.tools import html
from ...core.tools import filesystem as fs

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
        super(HTMLPageComponent, self).__init__(*args, **kwargs)

        # Tables
        self.info_table = None
        self.properties_table = None
        self.status_table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

        # Make tables
        self.make_tables()

        # Make plots
        self.make_plots()

        # Generaet the html
        self.generate()

        # Write
        self.write()

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

        # Make info table
        self.make_info_table()

        # Make properties table
        self.make_properties_table()

        # Make status table
        self.make_status_table()

    # -----------------------------------------------------------------

    def make_info_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the info table ...")

        # Create the table
        self.info_table = html.SimpleTable(self.galaxy_info.items(), header_row=["Property", "Value"])

    # -----------------------------------------------------------------

    def make_properties_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the properties table ...")

        # Create the table
        self.properties_table = html.SimpleTable(self.galaxy_properties.as_tuples(), header_row=["Property", "Value"])

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
        self.status_table = html.SimpleTable(self.status, header_row=["Step", "Status"], bgcolors=bgcolors)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the status
        self.generate_status()

    # -----------------------------------------------------------------

    def generate_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the status page ...")

        # Create title
        title = html.fontsize_template.format(size=20, text=html.underline_template.format(text="Modeling of " + self.galaxy_name))

        # Create titles
        title_info = html.underline_template.format(text="GALAXY INFO")
        title_properties = html.underline_template.format(text="GALAXY PROPERTIES")
        title_status = html.underline_template.format(text="MODELING STATUS")

        body = title + html.newline + html.newline + title_info + html.newline + html.newline + str(self.info_table) + html.newline + html.newline
        body += title_properties + html.newline + html.newline + str(self.properties_table) + html.newline + html.newline
        body += title_status + html.newline + html.newline + str(self.status_table) + html.newline + html.newline

        # Create contents
        contents = dict()
        contents["title"] = "Modeling of " + self.galaxy_name
        contents["head"] = html.link_stylesheet_header_template.format(url=stylesheet_url)
        contents["body"] = body

        # Create the status page
        self.page = html.page_template.format(**contents)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write status page
        self.write_status_page()

    # -----------------------------------------------------------------

    def write_status_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing status page ...")

        # Write
        fs.write_text(self.status_page_path, self.page)

# -----------------------------------------------------------------
