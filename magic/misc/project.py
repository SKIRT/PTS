#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.project Contains the Projector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...modeling.basics.models import load_2d_model
from ...core.basics.configurable import Configurable
from ...modeling.basics.properties import GalaxyProperties
from ..basics.coordinatesystem import CoordinateSystem
from ...modeling.basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from ...modeling.basics.instruments import SimpleInstrument, FrameInstrument, FullInstrument, SEDInstrument

# -----------------------------------------------------------------

class Projector(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(Projector, self).__init__(config)

        # The galaxy properties object
        self.properties = None

        # The projections
        self.projections = dict()

        # The disk model
        self.disk = None

        # The coordinate system
        self.wcs = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Show
        #self.show()

        # Create the projections
        self.create_projections()

        # 3. Writing
        if self.config.output is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Projector, self).setup(**kwargs)

        # Determine path to properties file
        properties_path = fs.join(self.config.parameters, "properties.dat")

        # Create the galaxy properties object
        self.properties = GalaxyProperties.from_file(properties_path)

        # Load the disk model
        disk_path = fs.join(self.config.parameters, "disk.mod")
        self.disk = load_2d_model(disk_path)

        # Load the WCS
        self.wcs = CoordinateSystem.from_file(self.config.wcs)

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        disk_pa = self.disk.position_angle

        # Create the 'earth' projection system
        azimuth = 0.0
        self.projections["earth"] = GalaxyProjection.from_wcs(self.wcs, self.properties.center, self.properties.distance, self.properties.inclination, azimuth, disk_pa)

        # Create the face-on projection system
        self.projections["faceon"] = FaceOnProjection.from_wcs(self.wcs, self.properties.center, self.properties.distance)

        # Create the edge-on projection system
        self.projections["edgeon"] = EdgeOnProjection.from_wcs(self.wcs, self.properties.center, self.properties.distance)

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projection systems
        for name in self.projections:

            # Create the instrument from the projection system
            self.instruments[name] = SimpleInstrument.from_projection(self.projections[name])

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print(self.properties)

        print(self.components)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the components
        self.write_projections()

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Loop over the components
        for name in self.projections:

            # Determine the path
            path = self.output_path_file(name + ".proj")

            # Save the projection
            self.projections[name].saveto(path)

# -----------------------------------------------------------------
