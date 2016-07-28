#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.animation.imageblink Contains the ImageBlinkAnimation class.

# -----------------------------------------------------------------

# Import standard modules
import io
import numpy as np
import copy
import imageio

# Import the relevant PTS classes and modules
from ...core.basics.animation import Animation
from ..tools import plotting

# -----------------------------------------------------------------

class ImageBlinkAnimation(Animation):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ImageBlinkAnimation, self).__init__()

        # Set the number of frames per second
        self.fps = 2

        # Properties for the plotting
        self.scale = "log"
        self.interval = "pts"
        self.cmap = "viridis"
        self.vmin = None
        self.vmax = None

    # -----------------------------------------------------------------

    def add_image(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Determine the interval (depends on whether this is the first frame that is added or not)
        if self.vmin is None and self.vmax is None: interval = self.interval
        else: interval = (self.vmin, self.vmax)

        # Make a plot of the image
        buf = io.BytesIO()
        self.vmin, self.vmax = plotting.plot_box(image, path=buf, format="png", interval=interval, scale=self.scale, cmap=self.cmap)
        buf.seek(0)
        im = imageio.imread(buf)
        buf.close()

        self.add_frame(im)

# -----------------------------------------------------------------
