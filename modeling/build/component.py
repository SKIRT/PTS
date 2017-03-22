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
from .tables import ModelsTable, RepresentationsTable
from ...core.basics.map import Map
from ...core.basics.configuration import load_mapping
from ..basics.models import DeprojectionModel3D, load_3d_model
from ...core.tools.serialization import load_dict
from ...magic.core.frame import Frame
from .definition import ModelDefinition

# -----------------------------------------------------------------

model_map_filename = "map.fits"

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

        # Path to the models directory
        self.models_path = None

        # Path to the representations directory
        self.representations_path = None

        # Path to the models table
        self.models_table_path = None

        # Path to the representations table
        self.representations_table_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuildComponent, self).setup(**kwargs)

        # Determine the path to the models directory
        self.models_path = fs.create_directory_in(self.build_path, "models")

        # Determine the path to the models table
        self.models_table_path = fs.join(self.models_path, "models.dat")

        # Initialize the models table if necessary
        if not fs.is_file(self.models_table_path):
            table = ModelsTable()
            table.saveto(self.models_table_path)

        # Determine the path to the representations directory
        self.representations_path = fs.create_directory_in(self.build_path, "representations")

        # Determine the path to the representations table
        self.representations_table_path = fs.join(self.representations_path, "representations.dat")

        # Initialize the representations table if necessary
        if not fs.is_file(self.representations_table_path):
            table = RepresentationsTable()
            table.saveto(self.representations_table_path)

    # -----------------------------------------------------------------

    def create_deprojection_for_wcs(self, wcs, filename, scaleheight):

        """
        This function ...
        :param wcs:
        :param filename:
        :param scaleheight:
        :return:
        """

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        position_angle = self.disk_position_angle

        # Get center coordinate of galaxy
        galaxy_center = self.galaxy_properties.center

        # Create deprojection
        # wcs, galaxy_center, distance, pa, inclination, filepath, scale_height
        deprojection = DeprojectionModel3D.from_wcs(wcs, galaxy_center, distance, position_angle, inclination, filename, scaleheight)

        # Return the deprojection
        return deprojection

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

        # Create the deprojection
        return self.create_deprojection_for_wcs(reference_wcs, filename, scaleheight)

    # -----------------------------------------------------------------

    def get_model_definition(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        path = self.get_model_path(model_name)
        if not fs.is_directory(path): raise ValueError("Model does not exist")
        return ModelDefinition(model_name, path)

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

    def get_representation(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        path = self.get_representation_path(representation_name)
        if not fs.is_directory(path): raise ValueError("Representation does not exist")

    # -----------------------------------------------------------------

    def get_representation_path(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        return fs.join(self.representations_path, representation_name)

    # -----------------------------------------------------------------

    @property
    def representations_table(self):

        """
        This function ...
        :return:
        """

        return RepresentationsTable.from_file(self.representations_table_path)

    # -----------------------------------------------------------------

    @property
    def representation_names(self):

        """
        This function ...
        :return:
        """

        return self.representations_table.names

    # -----------------------------------------------------------------

    def representations_for_model(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.representations_table.representations_for_model(model_name)

# -----------------------------------------------------------------

def get_models_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "build", "models", "models.dat")

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

    if not fs.is_file(get_models_table_path(modeling_path)): return []
    return get_models_table(modeling_path).names

# -----------------------------------------------------------------

def get_model_path(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return fs.join(modeling_path, "build", "models", model_name)

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

def load_component(path, add_map=False):

    """
    This function
    :param path:
    :param add_map:
    :return:
    """

    # Create a map
    component = Map()

    # Set the name
    component.name = fs.name(path)

    # Load the parameters
    parameters_path = fs.join(path, "parameters.cfg")
    if fs.is_file(parameters_path):
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
        component.map_path = map_path
        if add_map:
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

def load_stellar_component(modeling_path, model_name, component_name, add_map=False):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :param add_map:
    :return:
    """

    # Determine the path
    path = get_stellar_component_path(modeling_path, model_name, component_name)

    # Load the component
    return load_component(path, add_map=add_map)

# -----------------------------------------------------------------

def load_dust_component(modeling_path, model_name, component_name, add_map=False):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    :param add_map:
    :return:
    """

    # Determine the path
    path = get_dust_component_path(modeling_path, model_name, component_name)

    # Load the component
    return load_component(path, add_map=add_map)

# -----------------------------------------------------------------

def get_representations_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "build", "representations", "representations.dat")

# -----------------------------------------------------------------

def get_representations_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_representations_table_path(modeling_path)
    return RepresentationsTable.from_file(path)

# -----------------------------------------------------------------

def get_representation_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return get_representations_table(modeling_path).names

# -----------------------------------------------------------------

def get_model_name_for_representation(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    return get_representations_table(modeling_path).model_for_representation(name)

# -----------------------------------------------------------------

def get_representation_path(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    return fs.join(modeling_path, "build", "representations", name)

# -----------------------------------------------------------------

def get_representation(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    model_name = get_model_name_for_representation(modeling_path, name)
    path = get_representation_path(modeling_path, name)

    # Create and return the representation
    from .representation import Representation
    return Representation(name, model_name, path)

# -----------------------------------------------------------------

def get_representations_for_model(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    return get_representations_table(modeling_path).representations_for_model(model_name)

# -----------------------------------------------------------------

def get_earth_projection_for_representation(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    return get_representation(modeling_path, name).earth_projection

# -----------------------------------------------------------------

def get_pixelscale_for_representation(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    return get_representation(modeling_path, name).earth_projection.pixelscale

# -----------------------------------------------------------------
