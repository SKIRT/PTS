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

# -----------------------------------------------------------------

class ModelBuilder(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelBuilder, self).__init__(config)

        # The maps
        self.old_stars = None
        self.young_stars = None
        self.ionizing_stars = None
        self.dust = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Build stars
        self.build_stars()

        # Build dust component
        self.build_dust()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelBuilder, self).setup()

    # -----------------------------------------------------------------

    def build_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the stellar components ...")

        # Create the builder
        builder = StarsBuilder()

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
        builder = DustBuilder()

        # Run
        builder.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """
    
        self.write_maps()

        self.write_deprojections()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_deprojections(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
