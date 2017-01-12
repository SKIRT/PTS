#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.animation.scatter Contains the ScatterAnimation class.

# -----------------------------------------------------------------

# Import standard modules
import io
import numpy as np
import copy
import imageio

# Import the relevant PTS classes and modules
from ...core.basics.animation import Animation
from ...core.plot.scatter import Scatter3DPlotter
from ...core.basics.scatter import Scatter3D

# -----------------------------------------------------------------

class ScatterAnimation(Animation):

    """
    This class ...
    """

    def __init__(self, x_limits, y_limits, z_limits):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ScatterAnimation, self).__init__()

        # Set the number of frames per second
        self.fps = 5

        # Properties
        self.x_limits = x_limits
        self.y_limits = y_limits
        self.z_limits = z_limits

        self.x_label = None
        self.y_label = None
        self.z_label = None

        self.density = True

        # The scatter data
        self.scatter = Scatter3D()

        # The plotter
        self._plotter = Scatter3DPlotter()

        # Add the scatter data
        self._plotter.add_data("Scatter data", self.scatter)

    # -----------------------------------------------------------------

    def add_point(self, x, y, z):

        """
        This function ...
        :return:
        """

        # Add a point to the scatter
        self.scatter.add_point(x, y, z)

        # Add a point to the plotter
        #self._plotter.add_point(x, y, z)

        buf = io.BytesIO()

        self._plotter.set_x_limits(self.x_limits[0], self.x_limits[1])
        self._plotter.set_y_limits(self.y_limits[0], self.y_limits[1])
        self._plotter.set_z_limits(self.z_limits[0], self.z_limits[1])
        if self.x_label is not None: self._plotter.set_x_label(self.x_label)
        if self.y_label is not None: self._plotter.set_y_label(self.y_label)
        if self.z_label is not None: self._plotter.set_z_label(self.z_label)
        self._plotter.format = "png"

        self._plotter.density = self.density

        # Run the scatter plotter
        self._plotter.run(buf)

        buf.seek(0)
        im = imageio.imread(buf)
        buf.close()
        self.add_frame(im)

        # Clear the scatter plotter
        self._plotter.clear_figure()

# -----------------------------------------------------------------
