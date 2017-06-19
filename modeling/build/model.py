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

class ModelBuilderBase(BuildComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ModelBuilderBase, self).__init__(*args, **kwargs)

        # The path for this model
        self.model_path = None

        # The path for the stellar components
        self.model_stellar_path = None

        # The path for the dust components
        self.model_dust_path = None

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This functio n...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelBuilderBase, self).setup(**kwargs)

        # Set the model path and create it
        self.model_path = fs.create_directory_in(self.models_path, self.model_name)

        # Set the path of the directory for the stellar components
        self.model_stellar_path = fs.create_directory_in(self.model_path, "stellar")

        # Set the path of the directory for the dust components
        self.model_dust_path = fs.create_directory_in(self.model_path, "dust")

# -----------------------------------------------------------------

class ModelBuilder(ModelBuilderBase, GalaxyModelingComponent):
    
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
        ModelBuilderBase.__init__(self, *args, **kwargs)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The scaleheight of the old stars
        self.old_scaleheight = None

        # Map paths
        self.old_stars_map_path = None
        self.young_stars_map_path = None
        self.ionizing_stars_map_path = None
        self.dust_map_path = None

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
        ModelBuilderBase.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

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

        # Set the scaleheight of the old stars
        self.old_scaleheight = builder.old_scaleheight

        # Set map paths
        self.old_stars_map_path = builder.old_stars_map_path
        self.young_stars_map_path = builder.young_stars_map_path
        self.ionizing_stars_map_path = builder.ionizing_stars_map_path

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
        builder.run(old_scaleheight=self.old_scaleheight)

        # Set map paths
        self.dust_map_path = builder.dust_map_path

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
        log.info("Writing the models table ...")

        # Add the model
        table = self.models_table
        table.add_model(self.model_name, self.config.description, self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path, self.dust_map_path)

        # Save the table
        table.saveto(self.models_table_path)

# -----------------------------------------------------------------
