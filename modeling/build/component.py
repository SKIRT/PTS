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
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from .tables import ModelsTable, RepresentationsTable
from ...core.basics.map import Map
from ...core.basics.configuration import open_mapping
from ..basics.models import DeprojectionModel3D, load_3d_model
from ...core.tools.serialization import load_dict
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from .suite import ModelSuite

# -----------------------------------------------------------------

class BuildComponent(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(BuildComponent, self).__init__(*args, **kwargs)

        # The model suite
        self.suite = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuildComponent, self).setup(**kwargs)

        # Create the model suite
        self.suite = ModelSuite(self.build_path)

    # -----------------------------------------------------------------

    @property
    def models_path(self):

        """
        This function ...
        :return:
        """

        return self.suite.models_path

    # -----------------------------------------------------------------

    @property
    def representations_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.suite.representations_path

    # -----------------------------------------------------------------

    @property
    def models_table_path(self):

        """
        This function ...
        :return:
        """

        return self.suite.models_table_path

    # -----------------------------------------------------------------

    @property
    def representations_table_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.suite.representations_table_path

    # -----------------------------------------------------------------

    def get_model_definition(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.get_model_definition(model_name)

    # -----------------------------------------------------------------

    @property
    def models_table(self):

        """
        This function ...
        :return:
        """

        return self.suite.models_table

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return:
        """

        return self.suite.model_names

    # -----------------------------------------------------------------

    def get_model_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.get_model_path(model_name)

    # -----------------------------------------------------------------

    def get_model_stellar_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.get_model_stellar_path(model_name)

    # -----------------------------------------------------------------

    def get_model_dust_path(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.get_model_dust_path(model_name)

    # -----------------------------------------------------------------

    def get_stellar_component_names(self, model_name):

        """
        This function ...
        :return:
        """

        return self.suite.get_stellar_component_names(model_name)

    # -----------------------------------------------------------------

    def get_dust_component_names(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.get_dust_component_names(model_name)

    # -----------------------------------------------------------------

    def get_representation(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        return self.suite.get_representation(representation_name)

    # -----------------------------------------------------------------

    def get_representation_path(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        return self.suite.get_representation_path(representation_name)

    # -----------------------------------------------------------------

    @property
    def representations_table(self):

        """
        This function ...
        :return:
        """

        return self.suite.representations_table

    # -----------------------------------------------------------------

    @property
    def representation_names(self):

        """
        This function ...
        :return:
        """

        return self.suite.representation_names

    # -----------------------------------------------------------------

    def representations_for_model(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        return self.suite.representations_for_model(model_name)

    # -----------------------------------------------------------------

    def create_deprojection_for_wcs(self, galaxy_properties, disk_position_angle, wcs, filename, scaleheight):

        """
        This function ...
        :param galaxy_properties: has to be passed since this class is GENERAL (BUT THIS FUNCTION CAN ONLY BE CALLED FOR A GALAXY MODELING ENVIRONMENT)
        :param disk_position_angle:
        :param wcs:
        :param filename:
        :param scaleheight:
        :return:
        """

        return self.suite.create_deprojection_for_wcs(galaxy_properties, disk_position_angle, wcs, filename, scaleheight)

    # -----------------------------------------------------------------

    def create_deprojection_for_map(self, galaxy_properties, disk_position_angle, map, filename, scaleheight):

        """
        This function ...
        :param galaxy_properties:
        :param disk_position_angle:
        :param map:
        :param filename:
        :param scaleheight
        :return:
        """

        return self.suite.create_deprojection_for_map(galaxy_properties, disk_position_angle, map, filename, scaleheight)

# -----------------------------------------------------------------

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

def get_model_definition(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    from .definition import ModelDefinition

    path = get_model_path(modeling_path, model_name)
    if not fs.is_directory(path): raise ValueError("Model '" + model_name + "' does not exist")
    return ModelDefinition(model_name, path)

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
    parameters_path = fs.join(path, parameters_filename)
    if fs.is_file(parameters_path):
        parameters = open_mapping(parameters_path)
        component.parameters = parameters

    # Load the deprojection
    deprojection_path = fs.join(path, deprojection_filename)
    if fs.is_file(deprojection_path):
        deprojection = DeprojectionModel3D.from_file(deprojection_path)
        component.deprojection = deprojection

    # Load the map
    map_path = fs.join(path, model_map_filename)
    if fs.is_file(map_path):
        component.map_path = map_path
        if add_map:
            map = Frame.from_file(map_path)
            component.map = map

    # Load the model
    model_path = fs.join(path, model_filename)
    if fs.is_file(model_path):
        model = load_3d_model(model_path)
        component.model = model

    # Load the properties
    properties_path = fs.join(path, properties_filename)
    if fs.is_file(properties_path):
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

def load_stellar_component_deprojection(modeling_path, model_name, component_name):

    """
    This function ...
    :return: 
    """

    from ..component.galaxy import get_disk_position_angle

    # Load galaxy properties
    from ..component.galaxy import get_galaxy_properties
    properties = get_galaxy_properties(modeling_path)

    # Load component
    component = load_stellar_component(modeling_path, model_name, component_name, add_map=False)

    ## Set deprojection
    if "deprojection" in component:

        # Get title
        title = component.parameters.title

        # Return
        return title, component.deprojection

    # Check if this is a new component, add geometry, SED and normalization all at once
    if "geometry" in component.parameters:

        # Get title
        title = component.parameters.title

        # Check whether this is a read FITS geometry
        geometry_type = component.parameters.geometry
        if geometry_type != "ReadFitsGeometry": return component.parameters.title, None

        # Get properties for each of the three classes
        geometry_properties = component.properties["geometry"]

        # Get the path of the input map
        filepath = geometry_properties["filename"]

        # Get the scale height
        scale_height = geometry_properties["axialScale"]

        # Get properties
        wcs = CoordinateSystem.from_file(filepath)

        # Get the galaxy distance, the inclination and position angle
        distance = properties.distance
        inclination = properties.inclination
        position_angle = get_disk_position_angle(modeling_path)
        # Get center coordinate of galaxy
        galaxy_center = properties.center

        # Create
        deprojection = DeprojectionModel3D.from_wcs(wcs, galaxy_center, distance, position_angle, inclination, filepath, scale_height)

        # Return
        return title, deprojection

    # No deprojection
    return component.parameters.title, None

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

def load_dust_component_deprojection(modeling_path, model_name, component_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :param component_name:
    """

    from ..component.galaxy import get_disk_position_angle

    # Load galaxy properties
    from ..component.galaxy import get_galaxy_properties
    properties = get_galaxy_properties(modeling_path)

    # Load the component
    component = load_dust_component(modeling_path, model_name, component_name, add_map=False)

    # Set deprojection
    if "deprojection" in component:

        # Get title
        title = component.parameters.title

        # Return
        return title, component.deprojection

    # Check if this is a new dust component, add geometry, mix and normalization all at once
    if "geometry" in component.parameters:

        # Get title
        title = component.parameters.title

        # Check whether this is a read FITS geometry
        geometry_type = component.parameters.geometry
        if geometry_type != "ReadFitsGeometry": return title, None

        # Get properties for each of the three classes
        geometry_properties = component.properties["geometry"]

        # Get the path of the input map
        filepath = geometry_properties["filename"]

        # Get the scale height
        scale_height = geometry_properties["axialScale"]

        # Get properties
        wcs = CoordinateSystem.from_file(filepath)

        # Get the galaxy distance, the inclination and position angle
        distance = properties.distance
        inclination = properties.inclination
        position_angle = get_disk_position_angle(modeling_path)
        # Get center coordinate of galaxy
        galaxy_center = properties.center

        # Create
        deprojection = DeprojectionModel3D.from_wcs(wcs, galaxy_center, distance, position_angle, inclination, filepath, scale_height)

        # Return
        return title, deprojection

    # No deprojection for this component
    return component.parameters.title, None

# -----------------------------------------------------------------

def get_models_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "build", models_name)

# -----------------------------------------------------------------

def get_representations_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "build", representations_name)

# -----------------------------------------------------------------

def get_definitions_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(get_models_path(modeling_path), models_table_filename)

# -----------------------------------------------------------------

def get_definitions_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return

def get_representations_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "build", representations_name, representations_table_filename)

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

    return fs.join(modeling_path, "build", representations_name, name)

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
