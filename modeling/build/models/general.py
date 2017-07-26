#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.models.general Contains the GeneralBuilder class, a base class for StarsBuilder and DustBuilder.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..component import BuildComponent
from ....core.tools.logging import log
from ....core.prep.smile import SKIRTSmileSchema
from ....core.tools import filesystem as fs
from ....core.tools.serialization import write_dict
from ..suite import parameters_filename, deprojection_filename, model_map_filename, model_filename, properties_filename

# -----------------------------------------------------------------

class GeneralBuilder(BuildComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(GeneralBuilder, self).__init__(*args, **kwargs)

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

        # File paths
        self.parameter_paths = dict()
        self.deprojection_paths = dict()
        self.map_paths = dict()
        self.model_paths = dict()
        self.properties_paths = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GeneralBuilder, self).setup(**kwargs)

        # Create the SKIRT smile schema
        self.smile = SKIRTSmileSchema()

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
            path = fs.join(self.paths[name], parameters_filename)
            self.parameters[name].saveto(path)

            # Set path
            self.parameter_paths[name] = path

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
            path = fs.join(self.paths[name], deprojection_filename)
            self.deprojections[name].saveto(path)

            # Set path
            self.deprojection_paths[name] = path

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
            path = fs.join(self.paths[name], model_map_filename)
            self.maps[name].saveto(path)

            # Set path
            self.map_paths[name] = path

    # -----------------------------------------------------------------

    def write_models(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Writing the models ...")

        # Loop over the names
        for name in self.parameters:

            # Save model
            if name not in self.models: continue

            # Save the model
            path = fs.join(self.paths[name], model_filename)
            self.models[name].saveto(path)

            # Set path
            self.model_paths[name] = path

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
            path = fs.join(self.paths[name], properties_filename)
            write_dict(self.properties[name], path)

            # Set path
            self.properties_paths[name] = path

# -----------------------------------------------------------------
