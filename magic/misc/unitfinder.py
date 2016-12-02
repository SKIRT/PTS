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

# -----------------------------------------------------------------

class UnitFinder(Configurable):

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
        super(UnitFinder, self).__init__(config)

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        super(UnitFinder, self).setup(**kwargs)

# -----------------------------------------------------------------
