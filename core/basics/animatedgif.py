#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.animatedgif Contains the AnimatedGif class.

# -----------------------------------------------------------------

# Import standard modules
import imageio

# -----------------------------------------------------------------

class AnimatedGif(object):

    """
    This class ...
    """

    def __init__(self, frames=None):

        """
        The constructor ...
        :param frames:
        """

        # Set the frames
        self.frames = frames if frames is not None else []

        # The number of frames per second
        self.fps = 10

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        frames = imageio.mimread(path)
        return cls(frames)

    # -----------------------------------------------------------------

    def add_frame(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Add the frame
        self.frames.append(frame)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :param fps:
        :return:
        """

        # Create and write the GIF file
        imageio.mimwrite(path, self.frames, fps=self.fps)

# -----------------------------------------------------------------
