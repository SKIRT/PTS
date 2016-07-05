#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.inspector Contains the SourceInspector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import OldConfigurable
from ...core.tools.logging import log
from ..view import MagicViewer

# -----------------------------------------------------------------

class SourceInspector(OldConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SourceInspector, self).__init__(config, "magic")

        # -- Attributes --

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        """

        # Create a new SourceFinder instance
        if arguments.config is not None: inspector = cls(arguments.config)
        elif arguments.settings is not None: inspector = cls(arguments.settings)
        else: inspector = cls()

        # Return the new instance
        return inspector

    # -----------------------------------------------------------------

    def run(self, frame, galactic_catalog, stellar_catalog, special_region=None, ignore_region=None, bad_mask=None, animation=None):

        """
        This function ...
        :param frame:
        :param galactic_catalog:
        :param stellar_catalog:
        :param special_region:
        :param ignore_region:
        :param bad_mask:
        :param animation:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galactic_catalog, stellar_catalog, special_region, ignore_region, bad_mask, animation)

        # 2. Find the galaxies
        self.find_galaxies()
        
        # 3. Find the stars
        if self.config.find_stars: self.find_stars()

        # 4. Look for other sources
        if self.config.find_other_sources: self.find_other_sources()

        # 5. Build and update catalog
        self.build_and_synchronize_catalog()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Base class implementation removes the children
        super(SourceInspector, self).clear()

        # Set default values for all attributes

    # -----------------------------------------------------------------

    def setup(self, frame, galactic_catalog, stellar_catalog, special_region, ignore_region, bad_mask=None, animation=None):

        """
        This function ...
        :param frame:
        :param galactic_catalog:
        :param stellar_catalog:
        :param special_region:
        :param ignore_region:
        :param bad_mask:
        :param animation:
        :return:
        """

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(SourceInspector, self).setup()

# -----------------------------------------------------------------
