#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.setup Contains the ModelingSetupTool class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable
from ..core.tools.logging import log
from ..magic.tools.catalogs import get_ngc_name, get_hyperleda_name
from ..core.tools import filesystem as fs
from .core.component import get_config_file_path
from ..core.basics.configuration import Configuration

# -----------------------------------------------------------------

class ModelingSetupTool(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(ModelingSetupTool, self).__init__(config)

        # The path to the modeling directory
        self.modeling_path = None

        # The NGC name of the galaxy
        self.ngc_name = None

        # The HYPERLEDA name of the galaxy
        self.hyperleda_name = None

        # The modeling configuration
        self.modeling_config = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Resolve the name of the galaxy
        self.resolve_name()

        # 3. Create the modeling directory
        self.create_directory()

        # 4. Create the configuration
        self.create_config()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingSetupTool, self).setup()

        # Set the path to the modeling directory
        self.modeling_path = fs.join(self.config.path, self.config.galaxy_name)

    # -----------------------------------------------------------------

    def create_directory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling directory ...")

        # Check whether a directory with this name is already present
        if fs.is_directory(self.modeling_path): raise ValueError("A directory with the name '" + self.config.galaxy_name + "' already exists in the current working directory")

        # Create the directory
        fs.create_directory(self.modeling_path)

    # -----------------------------------------------------------------

    def resolve_name(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Resolving the galaxy name ...")

        # Get the NGC name of the galaxy
        self.ngc_name = get_ngc_name(self.config.galaxy_name)

        # Inform the user
        log.info("Galaxy NGC ID is '" + self.ngc_name + "'")

        # Get the name in the HYPERLEDA catalog
        self.hyperleda_name = get_hyperleda_name(self.config.galaxy_name)

        # Inform the user
        log.info("Galaxy HYPERLEDA ID is '" + self.hyperleda_name + "'")

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling configuration ...")

        # Create a Configuration object
        self.modeling_config = Configuration()

        # Set the configuration settings
        self.modeling_config.name = self.config.galaxy_name
        self.modeling_config.ngc_name = self.ngc_name
        self.modeling_config.hyperleda_name = self.hyperleda_name
        self.modeling_config.method = self.config.method
        self.modeling_config.host_ids = self.config.host_ids
        self.modeling_config.fitting_host_ids = self.config.fitting_host_ids

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the configuration
        self.write_config()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the modeling configuration ...")

        # Determine the path
        path = get_config_file_path(self.modeling_path)

        # Save the config
        self.modeling_config.saveto(path)

# -----------------------------------------------------------------
