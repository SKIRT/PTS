#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.preparation Contains the PreparationPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import html
from ..preparation.preparer import load_statistics, has_statistics
from ...core.basics.table import SmartTable
from .component import HTMLPageComponent, sortable_table_class
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width

# -----------------------------------------------------------------

page_width = 600

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

        # The statistics for each image
        self.statistics = dict()

        # The complete statistics table
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

        # Load the statistics
        self.load_statistics()

        # Make the statistics table
        self.make_statistics_table()

        # Make statistics tables
        #self.make_statistics_tables()

        # Make the table
        #self.make_statistics_table()

    # -----------------------------------------------------------------

    def load_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the statistics ...")

        # Loop over the preparation names
        for prep_name in self.preparation_names:

            # Check
            if not has_statistics(self.config.path, prep_name):

                log.warning("Statistics not found for the '" + prep_name + "' image")
                continue

            # Load the statistics
            statistics = load_statistics(self.config.path, prep_name)

            # Add
            self.statistics[prep_name] = statistics

    # -----------------------------------------------------------------

    def make_statistics_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the statistics table ...")

        # Get the statistics in a list for each image
        statistics = []
        for prep_name in self.preparation_names:
            if prep_name in self.statistics: statistics.append(self.statistics[prep_name])
            else: statistics.append(None)

        # Create the table
        table = SmartTable.from_composites(*statistics, labels=self.preparation_names, label="Image")

        # Create background colours
        #for i in range(table.nrows):
        #for fltr in self.environment.filters:
        colours = []
        for colour in self.environment.plotting_colours:
            colours.append([colour] + [None] * (table.ncolumns - 1))

        # Generate HTML table
        self.statistics_table = html.SimpleTable(table.as_tuples(), table.column_names, css_class=sortable_table_class, tostr_kwargs=self.tostr_kwargs, bgcolors=colours)

    # -----------------------------------------------------------------

    # def make_statistics_tables_old(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Making the statistics tables ...")
    #
    #     # Loop over the preparation names
    #     for prep_name in self.preparation_names:
    #
    #         # Load the statistics
    #         statistics = load_statistics(self.config.path, prep_name)
    #
    #         # Convert to table
    #         table = html.SimpleTable(statistics.as_tuples(), ["Property", "Value"], css_class=table_class, tostr_kwargs=self.tostr_kwargs)
    #
    #         # Set the table
    #         self.statistics_tables[prep_name] = table
    #
    # # -----------------------------------------------------------------
    #
    # def make_statistics_table_old(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Making the statistics table ...")
    #
    #     # Fill cells
    #     cells = []
    #     for name in self.preparation_names:
    #
    #         # Check whether statistics are found
    #         if name not in self.statistics_tables:
    #             log.warning("Preparation statistics not found for the '" + name + "' image")
    #             text = ""
    #
    #         else:
    #
    #         # Set title
    #         text = ""
    #         text += html.center_template.format(text=html.bold_template.format(text=name))
    #         text += html.newline + html.newline
    #         text += str(self.statistics_tables[name])
    #
    #         # Add to cells
    #         cells.append(text)
    #
    #     # Create the table
    #     self.statistics_table = html.SimpleTable.rasterize(cells, 4, css_class=table_class, tostr_kwargs=self.tostr_kwargs)

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
        log.info("Generating the page ...")

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

        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline + html.newline

        # Title
        self.page += "STATISTICS:"
        self.page += html.newline + html.newline

        # Add the table
        self.page += self.statistics_table

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
