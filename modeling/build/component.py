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
from .tables import RepresentationsTable
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
    def maps_table_path(self):

        """
        This function ...
        :return:
        """

        return self.suite.maps_table_path

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
    def maps_table(self):

        """
        This function ...
        :return:
        """

        return self.suite.maps_table

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

    def create_deprojection_for_map(self, galaxy_properties, disk_position_angle, map, filename, scaleheight, inclination=None):

        """
        This function ...
        :param galaxy_properties:
        :param disk_position_angle:
        :param map:
        :param filename:
        :param scaleheight:
        :param inclination:
        :return:
        """

        return self.suite.create_deprojection_for_map(galaxy_properties, disk_position_angle, map, filename, scaleheight, inclination=inclination)

# -----------------------------------------------------------------

def get_model_definition(modeling_path, model_name):

    """
    This function ...
    :param modeling_path:
    :param model_name:
    :return:
    """

    suite = ModelSuite.from_modeling_path(modeling_path)
    return suite.get_model_definition(model_name)

# -----------------------------------------------------------------

def get_representation(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    suite = ModelSuite.from_modeling_path(modeling_path)
    return suite.get_representation(name)

# -----------------------------------------------------------------

def get_models_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    from .suite import models_name
    return fs.join(modeling_path, "build", models_name)

# -----------------------------------------------------------------

def get_representations_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    from .suite import representations_name
    return fs.join(modeling_path, "build", representations_name)

# -----------------------------------------------------------------

def get_definitions_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    from .suite import models_table_filename
    return fs.join(get_models_path(modeling_path), models_table_filename)

# -----------------------------------------------------------------

def get_representations_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    from .suite import representations_name, representations_table_filename
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

    from .suite import representations_name
    return fs.join(modeling_path, "build", representations_name, name)

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
