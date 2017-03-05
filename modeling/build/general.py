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

        # The deprojections
        self.deprojections = dict()

        # The models
        self.models = dict()

        # The properties
        self.properties = dict()

        # The paths to the component directories
        self.paths = dict()

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

    def create_deprojection_for_map(self, map, filename, scaleheight):

        """
        This function ...
        :param map:
        :param filename:
        :param scaleheight
        :return:
        """

        # Get the WCS
        reference_wcs = map.wcs

        #filename = None
        #hz = None

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        pa = self.earth_projection.position_angle

        # Get center coordinate of galaxy
        galaxy_center = self.galaxy_properties.center

        # Create deprojection
        # wcs, galaxy_center, distance, pa, inclination, filepath, scale_height
        deprojection = DeprojectionModel3D.from_wcs(reference_wcs, galaxy_center, distance, pa, inclination, filename, scaleheight)

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
        self.write_component_directories()

        # Write parameters
        self.write_parameters()

        # Write deprojections
        self.write_deprojections()

        # Write maps
        self.write_maps()

        # Write models
        self.write_models()

        # Write properties
        self.write_properties()

    # -----------------------------------------------------------------

    def write_component_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the component directories ...")

        # Loop over the components
        for name in self.parameters:

            # Create a directory
            component_path = self.output_path_file(name)
            fs.create_directory(component_path)

            # Set the path
            self.paths[name] = component_path

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the component parameters ...")

        # Loop over the components
        for name in self.parameters:

            # Save parameters
            path = fs.join(self.paths[name], "parameters.cfg")
            self.parameters[name].saveto(path)

    # -----------------------------------------------------------------

    def write_deprojections(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Writing the deprojections ...")

        # Loop over the components
        for name in self.parameters:

            # Save deprojection
            if name not in self.deprojections: continue

            # Save
            path = fs.join(self.paths[name], "deprojection.mod")
            self.deprojections[name].saveto(path)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the components
        for name in self.parameters:

            # Save map
            if name not in self.maps: continue

            # Save the map
            path = fs.join(self.paths[name], "map.fits")
            self.maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_models(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Writing the models ...")

        for name in self.parameters:

            # Save model
            if name not in self.models: continue

            # Save the model
            path = fs.join(self.paths[name], "model.mod")
            self.models[name].saveto(path)

    # -----------------------------------------------------------------

    def write_properties(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Writing the component properties ...")

        # Loop over the parameters
        for name in self.parameters:

            # Save properties
            if name not in self.properties: continue

            # Write the properties
            path = fs.join(self.paths[name], "properties.dat")
            write_dict(self.properties[name], path)

# -----------------------------------------------------------------
