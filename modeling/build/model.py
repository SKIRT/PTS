#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.builder Contains the ModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import BuildComponent
from .dust import DustBuilder
from .stars import StarsBuilder
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..component.galaxy import GalaxyModelingComponent

# -----------------------------------------------------------------

class ModelBuilder(BuildComponent, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        #super(ModelBuilder, self).__init__(*args, **kwargs)
        ModelBuilder.__init__(self, *args, **kwargs)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The path for this model
        self.model_path = None

        # The path for the stellar components
        self.model_stellar_path = None

        # The path for the dust components
        self.model_dust_path = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Build stars
        self.build_stars()

        # 3. Build dust component
        self.build_dust()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(ModelBuilder, self).setup(**kwargs)
        ModelBuilder.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Set the model path and create it
        self.model_path = fs.create_directory_in(self.models_path, self.model_name)

        # Set the path of the directory for the stellar components
        self.model_stellar_path = fs.create_directory_in(self.model_path, "stellar")

        # Set the path of the directory for the dust components
        self.model_dust_path = fs.create_directory_in(self.model_path, "dust")

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def build_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the stellar components ...")

        # Create the builder
        builder = StarsBuilder(interactive=True)

        # Set the model name
        builder.config.name = self.model_name

        # Set the output path
        builder.config.output = self.model_stellar_path

        # Set SFR estimation
        builder.config.default_sfr = self.config.sfr

        # Run
        builder.run()

    # -----------------------------------------------------------------

    def build_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust components ...")

        # Create the builder
        builder = DustBuilder(interactive=True)

        # Set the model name
        builder.config.name = self.model_name

        # Set the output path
        builder.config.output = self.model_dust_path

        # Set dust mass estimation
        builder.config.default_dust_mass = self.config.dust_mass

        # Run
        builder.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write models table
        self.write_table()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model table ...")

        # Add the model
        table = self.models_table
        table.add_model(self.model_name, description, old_stars_map_path, young_stars_map_path, ionizing_stars_map_path, dust_map_path)

        # Save the table
        table.saveto(self.models_table_path)

# -----------------------------------------------------------------
