#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.general Contains the GeneralBuilder class, a base class for StarsBuilder and DustBuilder.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools.logging import log
from ...core.prep.smile import SKIRTSmileSchema
from ..basics.models import DeprojectionModel3D

# -----------------------------------------------------------------

class GeneralBuilder(BuildComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(GeneralBuilder, self).__init__(config, interactive)

        # The parameters
        self.parameters = dict()

        # The maps
        self.maps = dict()

        # The stellar components
        self.components = dict()

        # The deprojections
        self.deprojections = dict()

        # The SKIRT smile schema
        self.smile = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GeneralBuilder, self).setup()

        # Create the SKIRT smile schema
        self.smile = SKIRTSmileSchema()

    # -----------------------------------------------------------------

    def create_deprojection_for_map(self, map):

        """
        This function ...
        :param map:
        :return:
        """

        # Get the WCS
        reference_wcs = map.wcs

        filename = None
        hz = None

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        pa = self.earth_projection.position_angle

        # Get center coordinate of galaxy
        galaxy_center = self.galaxy_properties.center

        # Create deprojection
        # wcs, galaxy_center, distance, pa, inclination, filepath, scale_height
        deprojection = DeprojectionModel3D.from_wcs(reference_wcs, galaxy_center, distance, pa, inclination, filename, hz)

        # Return the deprojection
        return deprojection

    # -----------------------------------------------------------------

    @abstractmethod
    def build_additional(self):

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
        log.info("Writing ...")

# -----------------------------------------------------------------
