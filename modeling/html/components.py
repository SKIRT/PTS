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

from ipyvolume.embed import embed_html

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import HTMLPageComponent, table_class
from ...core.tools import html
from ...dustpedia.core.properties import DustPediaProperties
from ...dustpedia.core.database import get_account
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...core.filter.filter import parse_filter
from ...magic.core.remote import RemoteFrame
from ..plotting.model import plot_galaxy_components, generate_html

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

        self.bulge_table = None
        self.bulge_plot = None

        self.disk_table = None
        self.disk_plot = None

        self.model_table = None
        self.model_plot = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

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

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        self.plot_bulge()

        self.plot_disk()

        self.plot_model()

    # -----------------------------------------------------------------

    def plot_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the bulge component ...")

        # Plot, create HTML
        components = {"bulge": None}
        kwargs = {}
        box = plot_galaxy_components(components, draw=True, show=False, **kwargs)
        self.bulge_plot = generate_html(box, **kwargs)

    # -----------------------------------------------------------------

    def plot_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the disk component ...")

        # Plot, create HTML
        components = {"disk": None}
        kwargs = {}
        box = plot_galaxy_components(components, draw=True, show=False, **kwargs)
        self.disk_plot = generate_html(box, **kwargs)

    # -----------------------------------------------------------------

    def plot_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model (disk + bulge) ...")

        # Plot, create HTML
        components = {"disk": None, "bulge": None}
        kwargs = {}
        box = plot_galaxy_components(components, draw=True, show=False, **kwargs)
        self.model_plot = generate_html(box, **kwargs)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")



    # -----------------------------------------------------------------


    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        body = self.heading + html.newline



        body += self.footing

        # Make page
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

        return self.data_page_path

# -----------------------------------------------------------------
