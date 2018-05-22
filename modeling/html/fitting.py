#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.fitting Contains the FittingPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class
from ...core.tools import html

# -----------------------------------------------------------------

class FittingPageGenerator(HTMLPageComponent):

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
        super(FittingPageGenerator, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

        # The statistics table
        self.statistics_table = None

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
        super(FittingPageGenerator, self).setup(**kwargs)

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

        # Making the statistics table
        self.make_statistics_table()

    # -----------------------------------------------------------------

    def make_statistics_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the statistics table ...")

        # Get the statistics table
        table = self.fitting_run.statistics

        # Generate HTML table
        self.statistics_table = html.SimpleTable(table.as_tuples(), table.column_names, css_class=table_class, tostr_kwargs=self.tostr_kwargs)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

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

        return self.fitting_page_path

# -----------------------------------------------------------------
