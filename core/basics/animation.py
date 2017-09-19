#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.animation Contains the Animation class.

# -----------------------------------------------------------------

# Import standard modules
import copy
import imageio

# Import the relevant PTS classes and modules
from .apng import APNG

# -----------------------------------------------------------------

class Animation(object):

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

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        frames = imageio.mimread(path)
        animation = cls(frames)
        animation.path = path
        return animation

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

    def add_frame_from_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        frame = imageio.imread(path)
        self.add_frame(frame)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check if the path is defined
        if self.path is None: raise RuntimeError("Path is not defined for this animation")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :param fps:
        :return:
        """

        # APNG: special
        if path.endswith(".apng"): write_apng(path, self.frames)

        # Use ImageIO
        else: imageio.mimwrite(path, self.frames, fps=self.fps)

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return len(self.frames)

# -----------------------------------------------------------------

def invert_colors(animation):

    """
    This function ...
    :param animation:
    :return:
    """

    for frame in animation.frames: frame[:, :, 0:3] = 255 - frame[:, :, 0:3]

# -----------------------------------------------------------------

def inverted_colors(animation):

    """
    This function ...
    :param animation:
    :return:
    """

    animation_copy = copy.deepcopy(animation)
    invert_colors(animation_copy)
    return animation_copy

# -----------------------------------------------------------------

def write_apng(path, frames):

    """
    This funtion ...
    :param path:
    :param frames:
    :return:
    """

    raise NotImplementedError("Not supported yet")
    #im = APNG()

# -----------------------------------------------------------------
