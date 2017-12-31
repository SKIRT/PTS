#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.animation.distribution Contains the DistributionAnimation class.

# -----------------------------------------------------------------

# Import standard modules
import io
import numpy as np
import copy
import imageio
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.animation import Animation
from ...core.plot.distribution import DistributionPlotter
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

class DistributionAnimation(Animation):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, variable_name, label):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(DistributionAnimation, self).__init__()

        # Set the number of frames per second
        self.fps = 5

        # Properties
        self.min_value = min_value
        self.max_value = max_value

        self.variable_name = variable_name
        self.label = label

        # The plotter
        self._plotter = DistributionPlotter()

        self.values = []

        # Reference distributions
        self.reference_distributions = OrderedDict()

    # -----------------------------------------------------------------

    def add_reference_distribution(self, label, distribution):

        """
        This function ...
        :param distribution:
        :return:
        """

        # Add the distribution to the dictionary
        self.reference_distributions[label] = distribution

    # -----------------------------------------------------------------

    def add_value(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Add the value to the list
        self.values.append(value)

        # Create the new (normalized) distribution
        new_distribution = Distribution.from_values(self.variable_name, self.values)
        new_distribution.normalize(1.0, method="max")
        buf = io.BytesIO()

        # Add the reference distributions
        for label in self.reference_distributions: self._plotter.add_distribution(self.reference_distributions[label], label)

        # Add the new distribution
        self._plotter.add_distribution(new_distribution, self.label)
        self._plotter.set_variable_name(self.variable_name)
        self._plotter.run(buf, format="png", min_value=self.min_value, max_value=self.max_value, max_count=1., logscale=True)
        buf.seek(0)
        im = imageio.imread(buf)
        buf.close()
        self.add_frame(im)

        # Clear the plotter
        self._plotter.clear()

    # -----------------------------------------------------------------

    def add_point(self, x, y, z):

        """
        This function ...
        :return:
        """

        # Add a point to the plotter
        self._plotter.add_point(x, y, z)

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
