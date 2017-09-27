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
from abc import ABCMeta, abstractmethod, abstractproperty

# Import the relevant PTS classes and modules
from ..component import BuildComponent
from ....core.basics.log import log
from ....core.prep.smile import SKIRTSmileSchema
from ....core.tools import filesystem as fs
from ....core.tools.serialization import write_dict
from ..suite import parameters_filename, deprojection_filename, model_map_filename, model_filename, properties_filename
from ....magic.core.frame import Frame
from ....core.basics.configuration import save_mapping

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

    @property
    def component_names(self):

        """
        This function ...
        :return:
        """

        return self.parameters.keys()

    # -----------------------------------------------------------------

    @property
    def with_model(self):

        """
        This function ...
        :return:
        """

        # Loop over all components
        for name in self.component_names:

            # Check
            if name not in self.models: continue
            else: yield name

    # -----------------------------------------------------------------

    @property
    def with_deprojection(self):

        """
        This function ...
        :return:
        """

        # Loop over all components
        for name in self.component_names:

            # Check
            if name not in self.deprojections: continue
            else: yield name

    # -----------------------------------------------------------------

    @property
    def with_map(self):

        """
        This function ...
        :return:
        """

        # Loop over all components
        for name in self.component_names:

            # Check
            if name not in self.maps: continue
            else: yield name

    # -----------------------------------------------------------------

    @property
    def with_properties(self):

        """
        This function ...
        :return:
        """

        # Loop over all components
        for name in self.component_names:

            # Check
            if name not in self.properties: continue
            else: yield name

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
        for name in self.parameters: self.parameter_paths[name] = write_parameters(self.parameters[name], self.paths[name])

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
        for name in self.with_deprojection: self.deprojection_paths[name] = write_deprojection(self.deprojections[name], self.paths[name])

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the components
        for name in self.with_map: self.map_paths[name] = write_map(self.maps[name], self.paths[name])

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
        for name in self.with_model: self.model_paths[name] = write_model(self.models[name], self.paths[name])

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
        for name in self.with_properties: self.properties_paths[name] = write_properties(self.properties[name], self.paths[name])

    # -----------------------------------------------------------------

    @abstractproperty
    def additional_names(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def write_component_alt(directory, name, parameters, model, deprojection, map, properties):

    """
    This function ...
    :param directory:
    :param name:
    :param parameters:
    :param model:
    :param deprojection:
    :param map:
    :param properties:
    :return:
    """

    # Debugging
    log.debug("Writing the '" + name + "' component to the '" + directory + "' directory ...")

    # Determine the path for the component directory
    component_path = fs.join(directory, name)

    # Check
    if fs.is_directory(component_path): raise IOError("Already a directory with the name '" + name + "'")
    else: fs.create_directory(component_path)

    # Save parameters
    if parameters is not None: write_parameters(parameters, component_path)

    # Save model
    if model is not None: write_model(model, component_path)

    # Save deprojection
    if deprojection is not None: write_deprojection(deprojection, component_path)

    # Save map
    if map is not None: write_map(map, component_path)

    # Save properties
    if properties is not None: write_properties(properties, component_path)

    # Return the component path
    return component_path

# -----------------------------------------------------------------

def write_parameters(parameters, path):

    """
    Thisf unction ...
    :param parameters:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Writing the component parameters ...")

    # Determine the path
    parameters_path = fs.join(path, parameters_filename)

    # Save the parameters
    #parameters.saveto(parameters_path)
    save_mapping(parameters_path, parameters)

    # Return the path
    return parameters_path

# -----------------------------------------------------------------

def write_model(model, path):

    """
    This function ...
    :param model:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Writing the component model ...")

    # Determine the path
    model_path = fs.join(path, model_filename)

    # Save the model
    model.saveto(model_path)

    # Return the path
    return model_path

# -----------------------------------------------------------------

def write_deprojection(deprojection, path):

    """
    This function ...
    :param deprojection:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Writing the component deprojection ...")

    # Determine the path
    deprojection_path = fs.join(path, deprojection_filename)

    # Save the deprojection
    deprojection.saveto(deprojection_path)

    # Return the path
    return deprojection_path

# -----------------------------------------------------------------

def write_map(map, path):

    """
    This function ...
    :param map:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Writing the component map ...")

    # Determine the path
    map_path = fs.join(path, model_map_filename)

    # Save the map
    map.saveto(map_path)

    # Return the path
    return map_path

# -----------------------------------------------------------------

def write_properties(properties, path):

    """
    Thsi function ...
    :param properties:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Writing the component properties ...")

    # Determine the path
    properties_path = fs.join(path, properties_filename)

    # Write
    write_dict(properties, properties_path)

    # Return the path
    return properties_path

# -----------------------------------------------------------------

def write_component(directory, name, component):

    """
    This function ...
    :param directory: directory for the component to be placed (either dust or stellar directory of a certain model)
    :param name: name of the component
    :param component: the component mapping
    :return:
    """

    # Get
    parameters = component.parameters if "parameters" in component else None
    model = component.model if "model" in component else None
    deprojection = component.deprojection if "deprojection" in component else None
    map_path = component.map_path if "map_path" in component else None
    map = component.map if "map" in component else None
    properties = component.properties if "properties" in component else None

    # Load the map
    if map_path is not None and map is None: map = Frame.from_file(map_path)

    # Write, return the component path
    return write_component_alt(directory, name, parameters, model, deprojection, map, properties)

# -----------------------------------------------------------------
