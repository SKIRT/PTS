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

# Import astronomical modules
from photutils.datasets import make_100gaussians_image
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame

# -----------------------------------------------------------------

# https://photutils.readthedocs.io/en/stable/photutils/background.html

# -----------------------------------------------------------------

class SkyTest(TestImplementation):

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
        super(SkyTest, self).__init__(config, interactive)

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

        # Generate the data
        self.generate_data()

        self.reference()

        self.write()

        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SkyTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def generate_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the data ...")

        # Create the image
        data = make_100gaussians_image()
        self.frame = Frame(data)

    # -----------------------------------------------------------------

    def add_background(self):

        """
        This function ...
        :return:
        """

        ny, nx = data.shape
        y, x = np.mgrid[:ny, :nx]
        gradient = x * y / 5000.
        data2 = data + gradient
        plt.imshow(data2, norm=norm, origin='lower', cmap='Greys_r')

    # -----------------------------------------------------------------

    def reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, norm=norm, origin='lower', cmap='Greys_r')

# -----------------------------------------------------------------