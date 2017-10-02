#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.definition Contains the ModelDefinition class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .tables import ModelsTable
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.list import NamedCoordinateSystemList
from ...core.basics.containers import NamedFileList
#from .stars import basic_stellar_maps_names
#from .dust import basic_dust_maps_names
from .models.stars import basic_old_map_name, basic_young_map_name, basic_ionizing_map_name
from .models.dust import basic_dust_map_name
from .models.general import parameters_filename, deprojection_filename
from ...core.basics.configuration import open_mapping
from ..basics.models import DeprojectionModel3D
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class ModelDefinition(object):
    
    """
    This class...
    """

    def __init__(self, name, path, stellar_paths, dust_paths):

        """
        This function ..
        :param name:
        :param path:
        :param stellar_paths: stellar component paths (dictionary)
        :param dust_paths: dust component paths (dictionary)
        """

        # Set model name and directory path
        self.name = name
        self.path = path

        # Set the component paths
        self.stellar_paths = stellar_paths
        self.dust_paths = dust_paths

        # Subdirectories
        self.stellar_path = fs.create_directory_in(self.path, "stellar")
        self.dust_path = fs.create_directory_in(self.path, "dust")

        # Other input
        self.input_path = fs.create_directory_in(self.path, "input")

    # -----------------------------------------------------------------

    @property
    def own_stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.stellar_path, returns="name")

    # -----------------------------------------------------------------

    @property
    def own_dust_component_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.dust_path, returns="name")

    # -----------------------------------------------------------------

    @property
    def stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return self.stellar_paths.keys()

    # -----------------------------------------------------------------

    @property
    def dust_component_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_paths.keys()

    # -----------------------------------------------------------------

    def get_stellar_component_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.stellar_paths[name]

    # -----------------------------------------------------------------

    def get_dust_component_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_paths[name]

    # -----------------------------------------------------------------

    def get_stellar_component(self, name, add_map=False):

        """
        This function ...
        :param name:
        :param add_map:
        :return:
        """

        # Determine component_path
        #component_path = fs.join(self.stellar_path, name)
        component_path = self.get_stellar_component_path(name)

        from .suite import load_component
        return load_component(component_path, add_map=add_map)

    # -----------------------------------------------------------------

    def get_dust_component(self, name, add_map=False):

        """
        This function ...
        :param name:
        :param add_map:
        :return:
        """

        # Determine component path
        #component_path = fs.join(self.dust_path, name)
        component_path = self.get_dust_component_path(name)

        from .suite import load_component
        return load_component(component_path, add_map=add_map)

    # -----------------------------------------------------------------

    @property
    def own_stellar_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.stellar_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @property
    def own_dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @property
    def stellar_map_paths(self):

        """
        This function ...
        :return:
        """

        paths = []
        for name in self.stellar_paths:
            component_path = self.stellar_paths[name]
            paths.extend(fs.files_in_path(component_path, recursive=True, exact_name="map", extension="fits"))
        return paths

    # -----------------------------------------------------------------

    @property
    def dust_map_paths(self):

        """
        Thisf unction ...
        :return:
        """

        paths = []
        for name in self.dust_paths:
            component_path = self.dust_paths[name]
            paths.extend(fs.files_in_path(component_path, recursive=True, exact_name="map", extension="fits"))
        return paths

    # -----------------------------------------------------------------

    @property
    def paths_in_input(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.input_path)

    # -----------------------------------------------------------------

    @property
    def input_paths(self):

        """
        This function ...
        :return:
        """

        return self.paths_in_input + self.stellar_map_paths + self.dust_map_paths

    # -----------------------------------------------------------------

    @property
    def models_path(self):

        """
        This function ...
        :return: 
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def models_table_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.models_path, "models.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def models_table(self):

        """
        This function ...
        :return: 
        """

        return ModelsTable.from_file(self.models_table_path)

    # -----------------------------------------------------------------

    @property
    def description(self):

        """
        This function ...
        :return: 
        """

        return self.models_table.description_for_model(self.name)

    # -----------------------------------------------------------------

    @property
    def old_stars_map_name(self):

        """
        This function ...
        :return: 
        """

        return basic_old_map_name

    # -----------------------------------------------------------------

    @property
    def old_stars_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.old_stars_map_path))

    # -----------------------------------------------------------------

    def get_stellar_component_map_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        from .suite import model_map_basename

        # Get component path
        path = self.get_stellar_component_path(component_name)

        # Find the map path
        return fs.find_file_in_path(path, exact_name=model_map_basename, extension="fits", return_none=True)

    # -----------------------------------------------------------------

    @property
    def old_stars_map_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import old_component_name
        return self.get_stellar_component_map_path(old_component_name)

    # -----------------------------------------------------------------

    @property
    def young_stars_map_name(self):

        """
        This function ...
        :return: 
        """

        return basic_young_map_name

    # -----------------------------------------------------------------

    @property
    def young_stars_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.young_stars_map_path))

    # -----------------------------------------------------------------

    @property
    def young_stars_map_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import young_component_name
        return self.get_stellar_component_map_path(young_component_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_map_name(self):

        """
        This function ...
        :return: 
        """

        return basic_ionizing_map_name

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.ionizing_stars_map_path))

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_map_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import ionizing_component_name
        return self.get_stellar_component_map_path(ionizing_component_name)

    # -----------------------------------------------------------------

    @property
    def dust_map_name(self):

        """
        This function ...
        :return: 
        """

        return basic_dust_map_name

    # -----------------------------------------------------------------

    @property
    def dust_map_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.dust_map_path))

    # -----------------------------------------------------------------

    def get_dust_component_map_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        from .suite import model_map_basename

        # Get component path
        path = self.get_dust_component_path(component_name)

        # Return the map path
        return fs.find_file_in_path(path, exact_name=model_map_basename, extension="fits", return_none=True)

    # -----------------------------------------------------------------

    @property
    def dust_map_path(self):

        """
        This function ...
        :return: 
        """

        from .models.dust import disk_component_name
        return self.get_dust_component_map_path(disk_component_name)

    # -----------------------------------------------------------------

    @property
    def old_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.old_stars_map_path)

    # -----------------------------------------------------------------

    @property
    def young_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.young_stars_map_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.ionizing_stars_map_path)

    # -----------------------------------------------------------------

    @property
    def dust_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @property
    def basic_map_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_name, self.young_stars_map_name, self.ionizing_stars_map_name, self.dust_map_name]

    # -----------------------------------------------------------------

    @property
    def basic_component_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_component_name, self.young_stars_component_name, self.ionizing_stars_component_name, self.dust_map_component_name]

    # -----------------------------------------------------------------

    @property
    def basic_map_paths(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path, self.dust_map_path]

    # -----------------------------------------------------------------

    @property
    def basic_stellar_map_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_name, self.young_stars_map_name, self.ionizing_stars_map_name]

    # -----------------------------------------------------------------

    @property
    def basic_stellar_component_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_component_name, self.young_stars_component_name, self.ionizing_stars_component_name]

    # -----------------------------------------------------------------

    @property
    def basic_stellar_map_paths(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path]

    # -----------------------------------------------------------------

    @property
    def basic_component_names_dict(self):

        """
        This function ...
        :return: 
        """

        names = dict()
        for name, component_name in zip(self.basic_map_names, self.basic_component_names): names[name] = component_name
        return names

    # -----------------------------------------------------------------

    @property
    def basic_map_file_list(self):

        """
        This function ...
        :return: 
        """

        # Initialize the list
        lst = NamedFileList()

        # Old
        lst.append(self.old_stars_map_name, self.old_stars_map_path)

        # Young
        lst.append(self.young_stars_map_name, self.young_stars_map_path)

        # Ionizing
        lst.append(self.ionizing_stars_map_name, self.ionizing_stars_map_path)

        # Dust
        lst.append(self.dust_map_name, self.dust_map_path)

        # Return the list
        return lst

    # -----------------------------------------------------------------

    @property
    def basic_map_wcs_list(self):

        """
        This function ...
        :return: 
        """

        # Initialize the list
        lst = NamedCoordinateSystemList()

        # Old
        lst.append(self.old_stars_map_wcs, self.old_stars_map_name)

        # Young
        lst.append(self.young_stars_map_wcs, self.young_stars_map_name)

        # Ionizing
        lst.append(self.ionizing_stars_map_wcs, self.ionizing_stars_map_name)

        # Dust
        lst.append(self.dust_map_wcs, self.dust_map_name)

        # Return the list
        return lst

    # -----------------------------------------------------------------

    @property
    def basic_stellar_map_file_list(self):

        """
        This function ...
        :return: 
        """

        # Initialize the list
        lst = NamedFileList()

        # Old
        lst.append(self.old_stars_map_name, self.old_stars_map_path)

        # Young
        lst.append(self.young_stars_map_name, self.young_stars_map_path)

        # Ionizing
        lst.append(self.ionizing_stars_map_name, self.ionizing_stars_map_path)

        # Return the list
        return lst

    # -----------------------------------------------------------------

    @property
    def basic_dust_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return self.dust_map_wcs

    # -----------------------------------------------------------------

    @property
    def basic_stellar_map_wcs_list(self):

        """
        This function ...
        :return: 
        """

        # Initialize the list
        lst = NamedCoordinateSystemList()

        # Old
        lst.append(self.old_stars_map_wcs, self.old_stars_map_name)

        # Young
        lst.append(self.young_stars_map_wcs, self.young_stars_map_name)

        # Ionizing
        lst.append(self.ionizing_stars_map_wcs, self.ionizing_stars_map_name)

        # Return the list
        return lst

    # -----------------------------------------------------------------

    @property
    def basic_maps_minimum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_map_wcs_list
        return coordinate_systems.min_pixelscale

    # -----------------------------------------------------------------

    @property
    def basic_maps_maximum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_map_wcs_list
        return coordinate_systems.max_pixelscale

    # -----------------------------------------------------------------

    @property
    def basic_stellar_maps_minimum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_stellar_map_wcs_list
        return coordinate_systems.min_pixelscale

    # -----------------------------------------------------------------

    @property
    def basic_stellar_maps_maximum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_stellar_map_wcs_list
        return coordinate_systems.max_pixelscale

    # -----------------------------------------------------------------

    @property
    def basic_dust_map_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        return self.basic_dust_map_wcs.average_pixelscale

    # -----------------------------------------------------------------

    @property
    def old_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import old_component_name
        #return fs.directory_of(self.old_stars_map_path)
        return self.stellar_paths[old_component_name]

    # -----------------------------------------------------------------

    @property
    def young_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import young_component_name
        #return fs.directory_of(self.young_stars_map_path)
        return self.stellar_paths[young_component_name]

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import ionizing_component_name
        #return fs.directory_of(self.ionizing_stars_map_path)
        return self.stellar_paths[ionizing_component_name]

    # -----------------------------------------------------------------

    @property
    def dust_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.dust import disk_component_name
        #return fs.directory_of(self.dust_map_path)
        return self.dust_paths[disk_component_name]

    # -----------------------------------------------------------------

    @property
    def old_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.old_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @property
    def young_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.young_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.ionizing_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @property
    def dust_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.dust_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @property
    def old_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.old_stars_parameters_path)

    # -----------------------------------------------------------------

    @property
    def young_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.young_stars_parameters_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.ionizing_stars_parameters_path)

    # -----------------------------------------------------------------

    @property
    def dust_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.dust_parameters_path)

    # -----------------------------------------------------------------

    @property
    def old_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @property
    def young_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @property
    def dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.dust_parameters.mass

    # -----------------------------------------------------------------

    @property
    def old_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.old_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @property
    def young_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.young_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.ionizing_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @property
    def dust_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.dust_parameters.scale_height

    # -----------------------------------------------------------------

    @property
    def old_stars_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.old_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @property
    def young_stars_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.young_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_deprojection_path(self):

        """
        Tihs function ...
        :return: 
        """

        return fs.join(self.ionizing_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @property
    def dust_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.dust_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @property
    def old_stars_deprojection(self):

        """
        This function ...
        :return: 
        """

        return DeprojectionModel3D.from_file(self.old_stars_deprojection_path)

    # -----------------------------------------------------------------

    @property
    def young_stars_deprojection(self):

        """
        This function ...
        :return: 
        """

        return DeprojectionModel3D.from_file(self.young_stars_deprojection_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_deprojection(self):

        """
        This function ...
        :return: 
        """

        return DeprojectionModel3D.from_file(self.ionizing_stars_deprojection_path)

    # -----------------------------------------------------------------

    @property
    def dust_deprojection(self):

        """
        THis function ...
        :return: 
        """

        return DeprojectionModel3D.from_file(self.dust_deprojection_path)

    # -----------------------------------------------------------------

    def load_stellar_component(self, name, add_map=False):

        """
        Thisf unction ..
        :param name:
        :param add_map: 
        :return: 
        """

        from .suite import load_component
        path = self.get_stellar_component_path(name)
        return load_component(path, add_map=add_map)

    # -----------------------------------------------------------------

    def load_dust_component(self, name, add_map=False):

        """
        This function ...
        :param name:
        :param add_map:
        :return:
        """

        from .suite import load_component
        path = self.get_dust_component_path(name)
        return load_component(path, add_map=add_map)

    # -----------------------------------------------------------------

    def load_old_stars_component(self, add_map=False):

        """
        This function ...
        :param add_map: 
        :return: 
        """

        from .suite import load_component
        return load_component(self.old_stars_component_path, add_map=add_map)

    # -----------------------------------------------------------------

    def load_young_stars_component(self, add_map=False):

        """
        This function ...
        :param add_map: 
        :return: 
        """

        from .suite import load_component
        return load_component(self.young_stars_component_path, add_map=add_map)

    # -----------------------------------------------------------------

    def load_ionizing_stars_component(self, add_map=False):

        """
        This function ...
        :param add_map: 
        :return: 
        """

        from .suite import load_component
        return load_component(self.ionizing_stars_component_path, add_map=add_map)

    # -----------------------------------------------------------------

    def load_dust_disk_component(self, add_map=False):

        """
        This function ...
        :param add_map: 
        :return: 
        """

        from .suite import load_component
        return load_component(self.dust_component_path, add_map=add_map)

    # -----------------------------------------------------------------

    @property
    def additional_stellar_names(self):

        """
        This function ...
        :return:
        """

        from .models.stars import basic_stellar_component_names
        for name in self.stellar_component_names:
            if name in basic_stellar_component_names: continue
            else: yield name

    # -----------------------------------------------------------------

    @property
    def additional_dust_names(self):

        """
        This function ...
        :return:
        """

        from .models.dust import basic_dust_component_names
        for name in self.dust_component_names:
            if name in basic_dust_component_names: continue
            else: yield name

# -----------------------------------------------------------------
