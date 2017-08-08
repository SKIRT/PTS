#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.geometryplotter Contains the GeometryPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from textwrap import wrap
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as plt_Ellipse
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..basics.models import SersicModel3D, ExponentialDiskModel3D, DeprojectionModel3D

# -----------------------------------------------------------------

pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class GeometryPlotter(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(GeometryPlotter, self).__init__()

        # -- Attributes --

        # The geometries
        self.geometries = OrderedDict()

        # The patches
        self.patches = OrderedDict()

        # The figure
        self._figure = None

        self._min_x = None
        self._max_x = None
        self._min_y = None
        self._max_y = None

        # Properties
        self.title = None
        self.format = None
        self.transparent = False

    # -----------------------------------------------------------------

    def add_geometry(self, geometry, label):

        """
        This function ...
        :param geometry:
        :param label:
        :return:
        """

        self.geometries[label] = geometry

    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create matplotlib patches from the geometries
        self.create_patches()

        # Plot
        self.plot(path)

    # -----------------------------------------------------------------

    def create_patches(self):

        """
        This function ...
        :return:
        """

        colors = iter(pretty_colors)

        # Loop over the geometries
        for label in self.geometries:

            geometry = self.geometries[label]

            x_center = 0.0
            y_center = 0.0
            major = None # 2 * major axis radius
            minor = None # 2 * minor axis radius
            angle = None # in degrees

            if isinstance(geometry, SersicModel3D):

                major = 2.0 * geometry.effective_radius.to("pc").value
                minor = geometry.z_flattening * major
                angle = geometry.tilt.to("deg").value

            elif isinstance(geometry, ExponentialDiskModel3D):

                major = 2.0 * geometry.radial_scale.to("pc").value
                minor = 2.0 * geometry.axial_scale.to("pc").value
                angle = geometry.tilt.to("deg").value

            elif isinstance(geometry, DeprojectionModel3D):

                minor = 2.0 * geometry.scale_height.to("pc").value
                major = 0.3 * (geometry.pixelscale * geometry.x_size).to("pc").value
                angle = 0.0

            if self._min_x is None or 0.5*major > abs(self._min_x): self._min_x = - 0.5*major
            if self._max_x is None or 0.5*major > self._max_x: self._max_x = 0.5*major
            if self._min_y is None or 0.5*minor > abs(self._min_y): self._min_y = - 0.5*minor
            if self._max_y is None or 0.5*minor > self._max_y: self._max_y = 0.5*minor

            # Create the patch
            color = next(colors)
            ell = plt_Ellipse((x_center, y_center), major, minor, angle, edgecolor='none', facecolor=color, lw=3, alpha=0.7)

            # Add the patch
            self.patches[label] = ell

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Create the figure
        self._figure = plt.figure()
        ax = self._figure.add_subplot(111, aspect='equal')

        for label in self.patches:

            ax.add_patch(self.patches[label])

            # TODO: add text for label

        #plt.grid('on')

        ax.set_xlim(self._min_x, self._max_x)
        ax.set_ylim(self._min_y, self._max_y)

        # Set the title
        if self.title is not None: self._figure.suptitle("\n".join(wrap(self.title, 60)))

        # Finish
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def finish_plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        if type(path).__name__ == "BytesIO":
            log.debug("Saving the SED plot to a buffer ...")
        elif path is None: log.debug("Showing the SED plot ...")
        else: log.debug("Saving the SED plot to " + str(path) + " ...")

        # Save the figure
        if path is not None: plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------
