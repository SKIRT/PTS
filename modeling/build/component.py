#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.component Contains the BuildComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from .table import ModelsTable
from ...core.basics.map import Map
from ...core.basics.configuration import load_mapping
from ..basics.models import DeprojectionModel3D, load_3d_model
from ...core.tools.serialization import load_dict
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class BuildComponent(GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(BuildComponent, self).__init__(config, interactive)

        # Path to the models table
        self.models_table_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuildComponent, self).setup(**kwargs)

        # Determine the path to the models table
        self.models_table_path = fs.join(self.models_path, "models.dat")

        # Initialize the models table if necessary
        if not fs.is_file(self.models_table_path):
            table = ModelsTable()
            table.saveto(self.models_table_path)

    # -----------------------------------------------------------------

    @property
    def models_table(self):

        """
        This function ...
        :return:
        """

        # Open the table
        return ModelsTable.from_file(self.models_table_path)

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return:
        """

        return self.models_table.names

    # -----------------------------------------------------------------

    def get_model_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return fs.join(self.models_path, model_name)

    # -----------------------------------------------------------------

    def get_model_stellar_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return fs.join(self.get_model_path(model_name), "stellar")

    # -----------------------------------------------------------------

    def get_model_dust_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return fs.join(self.get_model_path(model_name), "dust")

    # -----------------------------------------------------------------

    def get_stellar_component_names(self, model_name):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.get_model_stellar_path(model_name), returns="name")

    # -----------------------------------------------------------------

    def get_dust_component_names(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return fs.directories_in_path(self.get_model_dust_path(model_name), returns="name")

# -----------------------------------------------------------------

def get_models_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "models", "models.dat")

# -----------------------------------------------------------------

def get_models_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return ModelsTable.from_file(get_models_table_path(modeling_path))

# -----------------------------------------------------------------

def get_model_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return get_models_table(modeling_path).names

# -----------------------------------------------------------------

def get_model_path(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.join(modeling_path, "models", model_name)

# -----------------------------------------------------------------

def get_model_stellar_path(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.join(get_model_path(modeling_path, model_name), "stellar")

# -----------------------------------------------------------------

def get_model_dust_path(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.join(get_model_path(modeling_path, model_name), "dust")

# -----------------------------------------------------------------

def get_stellar_component_names(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.directories_in_path(get_model_stellar_path(modeling_path, model_name), returns="name")

# -----------------------------------------------------------------

def get_dust_component_names(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.directories_in_path(get_model_dust_path(modeling_path, model_name), returns="name")

# -----------------------------------------------------------------

def get_stellar_component_path(modeling_path, model_name, component_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :return:
    """

    return fs.join(get_model_stellar_path(modeling_path, model_name), component_name)

# -----------------------------------------------------------------

def get_dust_component_path(modeling_path, model_name, component_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :return:
    """

    return fs.join(get_model_dust_path(modeling_path, model_name), component_name)

# -----------------------------------------------------------------

def load_component(path):

    """
    This function
    :param path:
    :return:
    """

    # Create a map
    component = Map()

    # Load the parameters
    parameters_path = fs.join(path, "parameters.cfg")
    parameters = Map()
    with open(parameters_path, "r") as fh:
        load_mapping(fh, parameters)
    component.parameters = parameters

    # Load the deprojection
    deprojection_path = fs.join(path, "deprojection.mod")
    if fs.is_file(deprojection_path):
        deprojection = DeprojectionModel3D.from_file(deprojection_path)
        component.deprojection = deprojection

    # Load the map
    map_path = fs.join(path, "map.fits")
    if fs.is_file(map_path):
        map = Frame.from_file(map_path)
        component.map = map

    # Load the model
    model_path = fs.join(path, "model.mod")
    if fs.is_file(model_path):
        model = load_3d_model(model_path)
        component.model = model

    # Load the properties
    properties_path = fs.join(path, "properties.dat")
    if fs.is_file(model_path):
        properties = load_dict(properties_path)
        component.properties = properties

    # Return the component
    return component

# -----------------------------------------------------------------

def load_stellar_component(modeling_path, model_name, component_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :return:
    """

    # Determine the path
    path = get_stellar_component_path(modeling_path, model_name, component_name)

    # Load the component
    return load_component(path)

# -----------------------------------------------------------------

def load_dust_component(modeling_path, model_name, component_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :return:
    """

    # Determine the path
    path = get_dust_component_path(modeling_path, model_name, component_name)

    # Load the component
    return load_component(path)

# -----------------------------------------------------------------
