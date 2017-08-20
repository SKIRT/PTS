#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.selection Contains the ComponentMapsSelection class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ..core.environment import GalaxyModelingEnvironment
from ..core.history import ModelingHistory
from pts.core.tools.utils import lazyproperty
from ...core.tools.utils import create_lazified_class

# -----------------------------------------------------------------

class ComponentMapsSelection(object):

    """
    This class...
    """

    def __init__(self, components_path):

        """
        The constructor ...
        :param components_path:
        :return:
        """

        # Check whether path exists
        if not fs.is_directory(components_path): raise IOError("The component maps directory '" + components_path + "' does not exist")

        # Set the components path
        self.components_path = components_path

    # -----------------------------------------------------------------

    @classmethod
    def from_modeling_path(cls, modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        components_path = fs.join(modeling_path, "maps", "components")
        return cls(components_path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(fs.directory_of(self.components_path))

    # -----------------------------------------------------------------

    @lazyproperty
    def environment(self):

        """
        This function ...
        :return:
        """

        return GalaxyModelingEnvironment(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def history(self):

        """
        This fucntion ...
        :return:
        """

        return ModelingHistory.from_file(self.environment.history_file_path)

    # -----------------------------------------------------------------

    @property
    def old_component_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_old_component_path

    # -----------------------------------------------------------------

    @property
    def young_component_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_young_component_path

    # -----------------------------------------------------------------

    @property
    def ionizing_component_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_ionizing_component_path

    # -----------------------------------------------------------------

    @property
    def dust_component_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_dust_component_path

    # -----------------------------------------------------------------

    @property
    def old_map_names(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.old_component_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @property
    def young_map_names(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.young_component_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.ionizing_component_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @property
    def dust_map_names(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_component_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    def get_old_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.old_component_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_young_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.young_component_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_ionizing_map_paths(self):

        """
        This funtion ...
        :return:
        """

        return fs.files_in_path(self.ionizing_component_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_dust_map_paths(self):

        """
        Thisj function ...
        :return:
        """

        return fs.files_in_path(self.dust_component_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    @property
    def old_map_paths(self):

        """
        This function ...
        :return:
        """

        # CAN BE LAZIFIED
        return self.get_old_map_paths()

    # -----------------------------------------------------------------

    @property
    def young_map_paths(self):

        """
        This function ...
        :return:
        """

        # CAN BE LAZIFIED
        return self.get_young_map_paths()

    # -----------------------------------------------------------------

    @property
    def ionizing_map_paths(self):

        """
        This function ...
        :return:
        """

        # CAN BE LAZIFIED
        return self.get_ionizing_map_paths()

    # -----------------------------------------------------------------

    @property
    def dust_map_paths(self):

        """
        This function ...
        :return:
        """

        # CAN BE LAZIFIED
        return self.get_dust_map_paths()

# -----------------------------------------------------------------

StaticComponentMapsSelection = create_lazified_class(ComponentMapsSelection, "StaticComponentMapsSelection")

# -----------------------------------------------------------------
