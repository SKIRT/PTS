#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.unitfinder Contains the UnitFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..catalog.importer import CatalogImporter
from ..core.frame import Frame

# -----------------------------------------------------------------

class UnitFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(UnitFinder, self).__init__(*args, **kwargs)

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get catalogs
        self.get_catalogs()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(UnitFinder, self).setup(**kwargs)

        # Load the frame (from config or input kwargs)
        if "frame" in kwargs: self.frames = kwargs.pop("frame")
        else: self.load_frame()

    # -----------------------------------------------------------------

    def load_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the frame ...")

        # Load the frame
        self.frame = Frame.from_file(self.config.dataset)

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Create a CatalogImporter instance
        catalog_importer = CatalogImporter()

        # Get the coordinate box and minimum pixelscale
        coordinate_box = self.frame.bounding_box
        pixelscale = self.frame.pixelscale

        # Run the catalog importer
        catalog_importer.run(coordinate_box=coordinate_box, pixelscale=pixelscale)

        # Set the catalogs
        self.galactic_catalog = catalog_importer.galactic_catalog
        self.stellar_catalog = catalog_importer.stellar_catalog

# -----------------------------------------------------------------
