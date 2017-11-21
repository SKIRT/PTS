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
from ...core.tools import filesystem as fs
from ..core.environment import GalaxyModelingEnvironment
from ..core.history import ModelingHistory
from ...core.tools.utils import lazyproperty
from ...core.tools.utils import create_lazified_class
from ...magic.core.frame import Frame

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
    def old_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_name
        return fs.join(self.old_component_path, deprojected_name)

    # -----------------------------------------------------------------

    @property
    def young_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_name
        return fs.join(self.young_component_path, deprojected_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_deprojected_path(self):

        """
        Thisn function ...
        :return:
        """

        from .components import deprojected_name
        return fs.join(self.ionizing_component_path, deprojected_name)

    # -----------------------------------------------------------------

    @property
    def dust_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_name
        return fs.join(self.dust_component_path, deprojected_name)

    # -----------------------------------------------------------------

    @property
    def old_skirt_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_skirt_name
        return fs.join(self.old_component_path, deprojected_skirt_name)

    # -----------------------------------------------------------------

    @property
    def young_skirt_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_skirt_name
        return fs.join(self.young_component_path, deprojected_skirt_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_skirt_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_skirt_name
        return fs.join(self.ionizing_component_path, deprojected_skirt_name)

    # -----------------------------------------------------------------

    @property
    def dust_skirt_deprojected_path(self):

        """
        This function ...
        :return:
        """

        from .components import deprojected_skirt_name
        return fs.join(self.dust_component_path, deprojected_skirt_name)

    # -----------------------------------------------------------------

    @property
    def old_edgeon_path(self):

        """
        This function ...
        :return:
        """

        from .components import edgeon_name
        return fs.join(self.old_component_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def young_edgeon_path(self):

        """
        This function ...
        :return:
        """

        from .components import edgeon_name
        return fs.join(self.young_component_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_edgeon_path(self):

        """
        This function ...
        :return:
        """

        from .components import edgeon_name
        return fs.join(self.ionizing_component_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def dust_edgeon_path(self):

        """
        This function ...
        :return:
        """

        from .components import edgeon_name
        return fs.join(self.dust_component_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def old_steps_path(self):

        """
        This fucntion ...
        :return:
        """

        from .components import steps_name
        return fs.join(self.old_component_path, steps_name)

    # -----------------------------------------------------------------

    @property
    def young_steps_path(self):

        """
        Thisf unction ...
        :return:
        """

        from .components import steps_name
        return fs.join(self.young_component_path, steps_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_steps_path(self):

        """
        This function ...
        :return:
        """

        from .components import steps_name
        return fs.join(self.ionizing_component_path, steps_name)

    # -----------------------------------------------------------------

    @property
    def dust_steps_path(self):

        """
        Thisf unction ...
        :return:
        """

        from .components import steps_name
        return fs.join(self.dust_component_path, steps_name)

    # -----------------------------------------------------------------

    def get_old_steps_path_for_map(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return fs.join(self.old_steps_path, name)

    # -----------------------------------------------------------------

    def get_young_steps_path_for_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.young_steps_path, name)

    # -----------------------------------------------------------------

    def get_ionizing_steps_path_for_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.ionizing_steps_path, name)

    # -----------------------------------------------------------------

    def get_dust_steps_path_for_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.dust_steps_path, name)

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

    def get_old_maps(self):

        """
        This function ...
        :return:
        """

        paths = self.get_old_map_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_young_maps(self):

        """
        This function ...
        :return:
        """

        paths = self.get_young_map_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        paths = self.get_ionizing_map_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_dust_maps(self):

        """
        This function ...
        :return:
        """

        paths = self.get_dust_map_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_old_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.old_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_young_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.young_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_ionizing_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.ionizing_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_dust_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_old_deprojected(self):

        """
        This function ...
        :return:
        """

        paths = self.get_old_deprojected_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_young_deprojected(self):

        """
        This function ...
        :return:
        """

        paths = self.get_young_deprojected_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_ionizing_deprojected(self):

        """
        Thisn function ...
        :return:
        """

        paths = self.get_ionizing_deprojected_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_dust_deprojected(self):

        """
        This function ...
        :return:
        """

        paths = self.get_dust_deprojected_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_old_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.old_skirt_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_young_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.young_skirt_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_ionizing_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.ionizing_skirt_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_dust_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_skirt_deprojected_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_old_deprojected_skirt(self):

        """
        Thisn function ...
        :return:
        """

        paths = self.get_old_deprojected_skirt_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_young_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        paths = self.get_young_deprojected_skirt_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_ionizing_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        paths = self.get_ionizing_deprojected_skirt_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_dust_deprojected_skirt(self):

        """
        Thisf unction ...
        :return:
        """

        paths = self.get_dust_deprojected_skirt_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_old_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.old_edgeon_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_young_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.young_edgeon_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_ionizing_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.ionizing_edgeon_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_dust_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_edgeon_path, returns="dict", extension="fits")

    # -----------------------------------------------------------------

    def get_old_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        paths = self.get_old_edgeon_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_young_edgeon(self):

        """
        This function ...
        :return:
        """

        paths = self.get_young_edgeon_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_ionizing_edgeon(self):

        """
        This function ...
        :return:
        """

        paths = self.get_ionizing_edgeon_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_dust_edgeon(self):

        """
        This function ...
        :return:
        """

        paths = self.get_dust_edgeon_paths()
        return paths_to_frames_dict(paths)

    # -----------------------------------------------------------------

    def get_old_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.old_component_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Old stellar map '" + map_name + "' does not exist in selection")
        return path

    # -----------------------------------------------------------------

    def get_old_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.old_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected old stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_old_skirt_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.old_skirt_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected with SKIRT old stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_old_edgeon_map_path(self, map_name):

        """
        This function ...
        :return:
        """

        path = fs.join(self.old_edgeon_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Edgeon old stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_young_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.young_component_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Young stellar map '" + map_name + "' does not exist in selection")
        return path

    # -----------------------------------------------------------------

    def get_young_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.young_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected young stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_young_skirt_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.young_skirt_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected with SKIRT young stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_young_edgeon_map_path(self, map_name):

        """
        THis function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.young_edgeon_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Edgeon young stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_ionizing_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.ionizing_component_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Ionizing stellar map '" + map_name + "' does not exist in selection")
        return path

    # -----------------------------------------------------------------

    def get_ionizing_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.ionizing_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected ionizing stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_ionizing_skirt_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.ionizing_skirt_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected with SKIRT ionizing stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_ionizing_edgeon_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.ionizing_edgeon_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Edgeon ionizing stellar map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_dust_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.dust_component_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Dust map '" + map_name + "' does not exist in selection")
        return path

    # -----------------------------------------------------------------

    def get_dust_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.dust_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected dust map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_dust_skirt_deprojected_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.dust_skirt_deprojected_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Deprojected with SKIRT dust map '" + map_name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def get_dust_edgeon_map_path(self, map_name):

        """
        This function ...
        :param map_name:
        :return:
        """

        path = fs.join(self.dust_edgeon_path, map_name + ".fits")
        if not fs.is_file(path): raise ValueError("Edgeon dust map '" + map_name + "' does not exist")
        return path

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

def paths_to_frames_dict(paths):

    """
    This function ...
    :param paths:
    :return:
    """

    maps = dict()
    for name in paths:
        path = paths[name]
        the_map = Frame.from_file(path)
        maps[name] = the_map
    return maps

# -----------------------------------------------------------------
