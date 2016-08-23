#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.inspector Contains the PreparationInspector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...magic.sources.finder import SourceFinder
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.misc.imageimporter import ImageImporter
from ...magic.catalog.importer import CatalogImporter
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...core.basics.animation import Animation
from ...core.tools import time
from ...core.tools import parsing

# -----------------------------------------------------------------

levels = [1.0, 3.0, 5.0, 10., 20.]

# -----------------------------------------------------------------

class PreparationInspector(PreparationComponent):
    
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
        super(PreparationInspector, self).__init__(config)

        # -- Attributes --

        self.inspect_path = None

        # Maps of the significance levels
        self.significance_maps = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Inspect errors
        self.inspect_significance()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PreparationInspector, self).setup()

        self.inspect_path = fs.join(self.config.path, "inspect")
        if not fs.is_directory(self.inspect_path): fs.create_directory(self.inspect_path)

    # -----------------------------------------------------------------

    def inspect_significance(self):

        """
        This function ...
        :return:
        """

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Calculating significance map for the " + name + " image ...")

            # Add the level map to the dictionary
            self.significance_maps[name] = self.dataset.get_significance(name, levels)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write significance maps
        self.write_significance_maps()

    # -----------------------------------------------------------------

    def write_significance_maps(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for name in self.significance_maps:

            # Determine the path
            path = fs.join(self.inspect_path, name + "_significance.fits")

            # Save the error level map
            self.significance_maps[name].save(path)

# -----------------------------------------------------------------
