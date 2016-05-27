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

        # Maximum value of the frame
        self.max_frame_value = None

        # Set the number of frames per second
        self.fps = 2

    # -----------------------------------------------------------------

    def add_image(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Create an animation to show the result of the source extraction step
        if self.max_frame_value is None: self.max_frame_value = np.nanmax(image)

        # Make a plot of the image
        buf = io.BytesIO()
        plotting.plot_box(image, path=buf, format="png", vmin=0.0, vmax=self.max_frame_value)
        buf.seek(0)
        im = imageio.imread(buf)
        buf.close()

        self.add_frame(im)

# -----------------------------------------------------------------
