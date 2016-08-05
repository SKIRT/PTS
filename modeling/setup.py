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
from ..magic.tools.catalogs import get_ngc_name
from ..core.tools import filesystem as fs

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

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Create the modeling directory
        self.create_directory()

        # 3. Create the meta file
        self.create_meta()

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

        if fs.is_directory(self.modeling_path): raise ValueError("A directory with the name '" + self.config.galaxy_name + "' already exists in the current working directory")

        # Create the directory
        fs.create_directory(self.modeling_path)

        # Inform the user
        log.info("Resolving the galaxy name ...")

    # -----------------------------------------------------------------

    def create_meta(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the meta file ...")

        # Get the NGC name of the galaxy
        ngc_name = get_ngc_name(self.config.galaxy_name)

        # Inform the user
        log.info("Galaxy NGC ID is '" + ngc_name + "'")

        # Determine the path to the meta file
        meta_path = fs.join(self.modeling_path, "modeling.meta")

        # Write the NGC name of the galaxy to the meta file
        with open(meta_path, 'w') as meta_file: print("Modeling directory for galaxy: " + ngc_name, file=meta_file)

# -----------------------------------------------------------------
