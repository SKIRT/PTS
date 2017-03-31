#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame
from pts.magic.core.mask import Mask

# -----------------------------------------------------------------

# https://photutils.readthedocs.io/en/stable/photutils/background.html

# -----------------------------------------------------------------

description = "Test the source extraction"

# -----------------------------------------------------------------

class SourcesTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(SourcesTest, self).__init__(config, interactive)

        # The frame
        self.frame = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # FInd
        self.find()

        # Extract
        self.extract()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SourcesTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------
