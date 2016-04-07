#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.animation Contains the FitAnimator

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from PIL import Image, ImageSequence
import imageio

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import gif, filesystem

# -----------------------------------------------------------------

class FitAnimator(FittingComponent):
    
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
        super(FitAnimator, self).__init__(config)

        # -- Attributes --

        # The list of frames
        self.frames = []

        # The animation
        self.animation = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new FitAnimator instance
        animator = cls(arguments.config)

        # Set the modeling path
        animator.config.path = arguments.path

        # Return the new instance
        return animator

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the images
        self.load_images()

        # 3. Create the animation
        self.create_animation()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FitAnimator, self).setup()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SED plot files ...")

        # Find all PNG files within the the fit/plot directory
        for path in filesystem.files_in_path(self.fit_plot_path, extension="png", recursive=True):

            # Load the image
            #image = Image.open(path)

            image = imageio.imread(path)

            # Add the image to the list of frames
            self.frames.append(image)

    # -----------------------------------------------------------------

    def create_animation(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the GIF animation ...")

        # Determine the path to the animation file
        path = self.full_output_path("fitting.gif")

        # Create and write the GIF file
        #gif.writeGif(path, self.frames, duration=0.2)

        imageio.mimwrite(path, self.frames, format="gif")

# -----------------------------------------------------------------
