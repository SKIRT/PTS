#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.builder Contains the ModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ..component.galaxy import GalaxyModelingComponent
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from .tables import ModelsTable, RepresentationsTable
from ...core.basics.map import Map
from ...core.basics.configuration import open_mapping
from ..basics.models import DeprojectionModel3D, load_3d_model
from ...core.tools.serialization import load_dict
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

parameters_filename = "parameters.cfg"
deprojection_filename = "deprojection.mod"
model_map_filename = "map.fits"
model_filename = "model.mod"
properties_filename = "properties.dat"

# -----------------------------------------------------------------

models_name = "models"
representations_name = "representations"

# -----------------------------------------------------------------

models_table_filename = "models.dat"
representations_table_filename = "representations.dat"

# -----------------------------------------------------------------

class ModelSuite(object):

    """
    This function ...
    """

    def __init__(self, path):

        """
        This function ...
        :param path:
        """

        # Set the build path
        self.path = path

        # Determine the path to the models directory
        self.models_path = fs.create_directory_in(self.build_path, models_name)

        # Determine the path to the models table
        self.models_table_path = fs.join(self.models_path, models_table_filename)

        # Initialize the models table if necessary
        if not fs.is_file(self.models_table_path):
            table = ModelsTable()
            table.saveto(self.models_table_path)

        # Determine the path to the representations directory
        self.representations_path = fs.create_directory_in(self.build_path, representations_name)

        # Determine the path to the representations table
        self.representations_table_path = fs.join(self.representations_path, representations_table_filename)

        # Initialize the representations table if necessary
        if not fs.is_file(self.representations_table_path):
            table = RepresentationsTable()
            table.saveto(self.representations_table_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_modeling_path(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        return cls(fs.join(path, "build"))

    # -----------------------------------------------------------------

    @property
    def build_path(self):

        """
        This function ...
        :return:
        """

        return self.path

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    def get_model_definition(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        from .definition import ModelDefinition

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
        else: return get_representation(self.modeling_path, representation_name)

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

        # Get the galaxy distance, the inclination and position angle
        distance = galaxy_properties.distance
        inclination = galaxy_properties.inclination
        position_angle = disk_position_angle

        # Get center coordinate of galaxy
        galaxy_center = galaxy_properties.center

        # Create deprojection
        # wcs, galaxy_center, distance, pa, inclination, filepath, scale_height
        deprojection = DeprojectionModel3D.from_wcs(wcs, galaxy_center, distance, position_angle, inclination, filename, scaleheight)

        # Return the deprojection
        return deprojection

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

        # Get the WCS
        reference_wcs = map.wcs

        # Create the deprojection
        return self.create_deprojection_for_wcs(galaxy_properties, disk_position_angle, reference_wcs, filename, scaleheight)

# -----------------------------------------------------------------
