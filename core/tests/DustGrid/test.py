#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log
from pts.modeling.tests.base import m81_data_path
from pts.core.tools import filesystem as fs
from pts.core.launch.launcher import SKIRTLauncher
from pts.core.prep.smile import SKIRTSmileSchema

# -----------------------------------------------------------------

description = "testing dust grid generation, writing and importing by SKIRT, optimization and plotting"

# -----------------------------------------------------------------

# Determine path to maps directory
maps_path = fs.join(m81_data_path, "maps")

# Determine the path to the dust map
dust_map_path = fs.join(maps_path, "dust.fits")

# -----------------------------------------------------------------

class DustGridTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DustGridTest, self).__init__(*args, **kwargs)

        # Smile
        self.smile = SKIRTSmileSchema()

        # Create the SKIRT launcher
        self.launcher = SKIRTLauncher()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create the ski file
        self.create_ski()

        # Launch
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustGridTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Create template ski file
        ski = self.smile.create_panchromatic_template()

        ski_path = fs.join(self.path, "template.ski")

        # Save
        ski.saveto(ski_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Launching ...")

# -----------------------------------------------------------------
