#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.all Contains the AllPagesGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.html.status import StatusPageGenerator
from pts.modeling.html.data import DataPageGenerator
from pts.modeling.html.preparation import PreparationPageGenerator
from pts.modeling.html.components import ComponentsPageGenerator
from pts.modeling.html.maps import MapsPageGenerator
from pts.modeling.html.model import ModelPageGenerator
from pts.core.tools.logging import log
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from ..core.progression import create_modeling_progression

# -----------------------------------------------------------------

class AllPagesGenerator(GalaxyModelingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(AllPagesGenerator, self).__init__(*args, **kwargs)

        # The modeling progression to use to generate the pages
        self.progression = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

        # Generate the status page
        if self.history.finished("fetch_properties"): self.generate_status()

        # Generate the data page
        if self.history.finished("fetch_images"): self.generate_data()

        # Generate the preparation page
        if self.history.finished("prepare_data"): self.generate_preparation()

        # Generate the components page
        if self.history.finished("decompose"): self.generate_components()

        # Generate the maps page
        if self.history.finished("build_model"): self.generate_maps()

        # GEnerate the model page
        if self.history.finished("configure_fit"): self.generate_model()

        # Write
        self.write()

        # Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AllPagesGenerator, self).setup(**kwargs)

        # Create the progression
        if "progression" in kwargs: self.progression = kwargs.pop("progression")
        else: self.create_progression()

    # -----------------------------------------------------------------

    def create_progression(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling progression ...")

        # Create
        self.progression = create_modeling_progression(self.config.path)

    # -----------------------------------------------------------------

    def generate_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generate the status page ...")

        # Generate
        generator = StatusPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the data page ...")

        # Generate
        generator = DataPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the preparation page ...")

        # Generate
        generator = PreparationPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the components page ...")

        # Generate
        generator = ComponentsPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the maps page ...")

        # Generate
        generator = MapsPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the model page ...")

        # Generate the model page
        generator = ModelPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the pages ...")

        # Open
        fs.open_in_browser(self.environment.html_status_path)

# -----------------------------------------------------------------
