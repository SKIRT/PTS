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
from .component import HTMLPageComponent, stylesheet_url, table_class
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class MapsPageGenerator(HTMLPageComponent):

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
        super(MapsPageGenerator, self).__init__(*args, **kwargs)

        # The plots
        self.old_plot = None
        self.young_plot = None
        self.ionizing_plot = None
        self.dust_plot = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Make tables
        self.make_tables()

        # 3. Make plots
        self.make_plots()

        # 4. Generaet the html
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
        super(MapsPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_model(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_name(self):

        """
        This function ...
        :return:
        """



    # -----------------------------------------------------------------

    @lazyproperty
    def old_component_path(self):

        """
        This function ...
        :return:
        """

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        self.make_old_plot()
n
    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Genreate the page
        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

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

        return self.maps_page_path

# -----------------------------------------------------------------
