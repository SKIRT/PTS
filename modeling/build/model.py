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

# -----------------------------------------------------------------

class ModelBuilder(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ModelBuilder, self).__init__(config, interactive)

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
        super(ModelBuilder, self).setup(**kwargs)

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

        # Set the output path
        builder.config.output = self.model_stellar_path

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

        # Set the output path
        builder.config.output = self.model_dust_path

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

        # Write the maps
        #self.write_maps()

        # Write the deprojections
        #self.write_deprojections()

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
        table.add_model(self.model_name, None, None, None, None, None)

        # Save the table
        table.saveto(self.models_table_path)

# -----------------------------------------------------------------
