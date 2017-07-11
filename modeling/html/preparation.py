#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.preparation Contains the StatusPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import HTMLPageComponent, table_class
from ...core.tools import html
from ..preparation.preparer import load_statistics
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

class PreparationPageGenerator(HTMLPageComponent):

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
        super(PreparationPageGenerator, self).__init__(*args, **kwargs)

        # Tables
        self.statistics_tables = dict()
        self.statistics_table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Make plots
        self.make_plots()

        # 3. Make tables
        self.make_tables()

        # 4. Generaet the html
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
        super(PreparationPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make statistics tables
        self.make_statistics_tables()

        # Make the table
        self.make_statistics_table()

    # -----------------------------------------------------------------

    def make_statistics_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the statistics tables ...")

        # Loop over the preparation names
        for prep_name in self.preparation_names:

            # Load the statistics
            statistics = load_statistics(self.config.path, prep_name)

            # Convert to table
            table = html.SimpleTable(statistics.as_tuples(), ["Property", "Value"], css_class=table_class, tostr_kwargs=self.tostr_kwargs)

            # Set the table
            self.statistics_tables[prep_name] = table

    # -----------------------------------------------------------------

    def make_statistics_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the statistics table ...")

        # Fill cells
        cells = []
        for name in self.preparation_names:

            # Check whether statistics are found


            # Set title
            text = ""
            text += html.center_template.format(text=html.bold_template.format(text=name))
            text += html.newline + html.newline
            text += str(self.statistics_tables[name])

            # Add to cells
            cells.append(text)

        # Create the table
        self.statistics_table = html.SimpleTable.rasterize(cells, 4, css_class=table_class, tostr_kwargs=self.tostr_kwargs)

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

        # Heading
        body = self.heading

        # Add table
        body += str(self.statistics_table)

        # Footing
        body += self.footing

        # Create the status page
        self.make_page(body)

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

        return self.preparation_page_path

# -----------------------------------------------------------------
