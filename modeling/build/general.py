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
from ...core.tools import filesystem as fs
from ...core.tools.serialization import write_dict

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
        #self.components = dict()

        # The deprojections
        self.deprojections = dict()

        # The models
        self.models = dict()

        # The properties
        self.properties = dict()

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

        # Write components
        self.write_components()

    # -----------------------------------------------------------------

    def write_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the components ...")

        # Loop over the components
        for name in self.parameters:

            # Create a directory
            component_path = self.output_path_file(name)
            fs.create_directory(component_path)

            # Save parameters
            path = fs.join(component_path, "parameters.cfg")
            self.parameters[name].saveto(path)

            # Save deprojection
            if name in self.deprojections:

                path = fs.join(component_path, "deprojection.mod")
                self.deprojections[name].saveto(path)

            # Save map
            if name in self.maps:

                path = fs.join(component_path, "map.fits")
                self.maps[name].saveto(path)

            # Save model
            if name in self.models:

                path = fs.join(component_path, "model.mod")
                self.models[name].saveto(path)

            # Save properties
            if name in self.properties:

                path = fs.join(component_path, "properties.dat")
                write_dict(self.properties[name], path)

# -----------------------------------------------------------------
