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
from .tables import ModelsTable, ModelMapsTable
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.list import NamedCoordinateSystemList
from ...core.basics.containers import NamedFileList
from .models.stars import basic_old_map_name, basic_young_map_name, basic_ionizing_map_name
from .models.dust import basic_dust_map_name
from .models.general import parameters_filename, deprojection_filename
from ...core.basics.configuration import open_mapping, save_mapping
from ..basics.models import DeprojectionModel3D
from ...core.tools.utils import lazyproperty
from ..core.mappings import Mappings
from ...core.filter.filter import parse_filter
from ...core.tools import numbers
from ...core.units.unit import parse_unit as u
from ..core.bruzualcharlot import create_bruzual_charlot_sed

# -----------------------------------------------------------------

stellar_name = "stellar"
dust_name = "dust"
input_name = "input"

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
        self.stellar_path = fs.create_directory_in(self.path, stellar_name)
        self.dust_path = fs.create_directory_in(self.path, dust_name)

        # Other input
        self.input_path = fs.create_directory_in(self.path, input_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get name
        name = fs.name(path)

        # Set paths
        stellar_path = fs.join(path, stellar_name)
        dust_path = fs.join(path, dust_name)

        # Get stellar and dust component paths
        stellar_paths = fs.directories_in_path(stellar_path, returns="dict")
        dust_paths = fs.directories_in_path(dust_path, returns="dict")

        # Return
        return cls(name, path, stellar_paths=stellar_paths, dust_paths=dust_paths)

    # -----------------------------------------------------------------

    @lazyproperty
    def own_stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.stellar_path, returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def own_dust_component_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.dust_path, returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return self.stellar_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_paths.keys()

    # -----------------------------------------------------------------

    def is_stellar_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.stellar_component_names

    # -----------------------------------------------------------------

    def is_dust_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.dust_component_names

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

    @lazyproperty
    def own_stellar_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.stellar_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def own_dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
    def paths_in_input(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.input_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def input_paths(self):

        """
        This function ...
        :return:
        """

        return self.paths_in_input + self.stellar_map_paths + self.dust_map_paths

    # -----------------------------------------------------------------

    @lazyproperty
    def models_path(self):

        """
        This function ...
        :return: 
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
    def maps_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.models_path, "maps.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_table(self):

        """
        This function ...
        :return:
        """

        return ModelMapsTable.from_file(self.maps_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def original_old_map_name(self):

        """
        This property ...
        """

        return self.maps_table.old_stars_map_name_for_model(self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def original_young_map_name(self):

        """
        This property ...
        """

        return self.maps_table.young_stars_map_name_for_model(self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def original_ionizing_map_name(self):

        """
        This property ...
        :return:
        """

        return self.maps_table.ionizing_stars_map_name_for_model(self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def original_dust_map_name(self):

        """
        This function ...
        :return:
        """

        return self.maps_table.dust_map_name_for_model(self.name)

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
    def young_stars_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.young_stars_map_path))

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
    def ionizing_stars_component_name(self):

        """
        This function ...
        :return: 
        """

        return fs.name(fs.directory_of(self.ionizing_stars_map_path))

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
    def dust_map_path(self):

        """
        This function ...
        :return: 
        """

        from .models.dust import disk_component_name
        return self.get_dust_component_map_path(disk_component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.old_stars_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.young_stars_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.ionizing_stars_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return CoordinateSystem.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_map_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_name, self.young_stars_map_name, self.ionizing_stars_map_name, self.dust_map_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_component_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_component_name, self.young_stars_component_name, self.ionizing_stars_component_name, self.dust_map_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_map_paths(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path, self.dust_map_path]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_stellar_map_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_name, self.young_stars_map_name, self.ionizing_stars_map_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_stellar_component_names(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_component_name, self.young_stars_component_name, self.ionizing_stars_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_stellar_map_paths(self):

        """
        This function ...
        :return: 
        """

        return [self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path]

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_component_names_dict(self):

        """
        This function ...
        :return: 
        """

        names = dict()
        for name, component_name in zip(self.basic_map_names, self.basic_component_names): names[name] = component_name
        return names

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
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

    @lazyproperty
    def basic_dust_map_wcs(self):

        """
        This function ...
        :return: 
        """

        return self.dust_map_wcs

    # -----------------------------------------------------------------

    @lazyproperty
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

    @lazyproperty
    def basic_maps_minimum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_map_wcs_list
        return coordinate_systems.min_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_maps_maximum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_map_wcs_list
        return coordinate_systems.max_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_stellar_maps_minimum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_stellar_map_wcs_list
        return coordinate_systems.min_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_stellar_maps_maximum_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        coordinate_systems = self.basic_stellar_map_wcs_list
        return coordinate_systems.max_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_dust_map_average_pixelscale(self):

        """
        This function ...
        :return: 
        """

        return self.basic_dust_map_wcs.average_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_component_path(self):

        """
        This function ...
        :return:
        """

        from .models.stars import bulge_component_name
        return self.stellar_paths[bulge_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import old_component_name
        return self.stellar_paths[old_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import young_component_name
        return self.stellar_paths[young_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.stars import ionizing_component_name
        return self.stellar_paths[ionizing_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component_path(self):

        """
        This function ...
        :return: 
        """

        from .models.dust import disk_component_name
        return self.dust_paths[disk_component_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.bulge_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.old_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.young_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.ionizing_stars_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_parameters_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.dust_component_path, parameters_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_parameters(self):

        """
        This function ...
        :return:
        """

        return open_mapping(self.bulge_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.old_stars_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.young_stars_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.ionizing_stars_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_sed(self):

        """
        This function ...
        :return:
        """

        return create_bruzual_charlot_sed(metallicity=self.metallicity, age=self.old_stars_age)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed(self):

        """
        This function ...
        :return:
        """

        return create_bruzual_charlot_sed(metallicity=self.metallicity, age=self.young_stars_age)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_sfr(self):

        """
        This function ...
        :return:
        """

        # Get scalar SFR from parameters
        sfr_scalar = self.ionizing_stars_parameters.sfr

        # Get the SFR
        sfr_from_lum = Mappings.sfr_for_luminosity(self.metallicity, self.ionizing_stars_compactness, self.ionizing_stars_pressure, self.ionizing_stars_covering_factor, self.ionizing_stars_luminosity, self.ionizing_stars_normalization_wavelength)

        # Check
        if not numbers.is_close(sfr_scalar, sfr_from_lum.to("Msun/yr").value): raise ValueError("Inconsistent SFR and FUV luminosity: " + str(sfr_scalar) + " =/= " + str(sfr_from_lum.to("Msun/yr").value))

        # Return
        return sfr_from_lum

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_pressure(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.pressure

    # -----------------------------------------------------------------

    @ionizing_stars_pressure.setter
    def ionizing_stars_pressure(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.ionizing_stars_parameters.pressure = value

        # Save
        save_mapping(self.ionizing_stars_parameters_path, self.ionizing_stars_parameters)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_covering_factor(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.covering_factor

    # -----------------------------------------------------------------

    @ionizing_stars_covering_factor.setter
    def ionizing_stars_covering_factor(self, value):

        """
        Thisf unction ...
        :param value:
        :return:
        """

        # Set the new value
        self.ionizing_stars_parameters.covering_factor = value

        # Save
        save_mapping(self.ionizing_stars_parameters_path, self.ionizing_stars_parameters)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_compactness(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.compactness

    # -----------------------------------------------------------------

    @ionizing_stars_compactness.setter
    def ionizing_stars_compactness(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.ionizing_stars_parameters.compactness = value

        # Save
        save_mapping(self.ionizing_stars_parameters_path, self.ionizing_stars_parameters)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_mappings(self):

        """
        This function ...
        :return:
        """

        return Mappings(self.metallicity, self.ionizing_stars_compactness, self.ionizing_stars_pressure, self.ionizing_stars_covering_factor, sfr=self.ionizing_stars_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_mappings_normalized(self):

        """
        This function ...
        :return:
        """

        return Mappings(self.metallicity, self.ionizing_stars_compactness, self.ionizing_stars_pressure, self.ionizing_stars_covering_factor)

    # -----------------------------------------------------------------

    @property
    def ionizing_sed(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_mappings.sed

    # -----------------------------------------------------------------

    @property
    def ionizing_sed_normalized(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_mappings_normalized.sed

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_parameters(self):

        """
        This function ...
        :return: 
        """

        return open_mapping(self.dust_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.bulge_parameters.luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_normalization_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter(self.bulge_parameters.filter)

    # -----------------------------------------------------------------

    @property
    def bulge_normalization_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.bulge_parameters.wavelength

    # -----------------------------------------------------------------

    @property
    def bulge_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.bulge_parameters.neutral_luminosity

    # -----------------------------------------------------------------

    @property
    def bulge_fluxdensity(self):

        """
        This function ...
        :return:
        """

        return self.bulge_parameters.fluxdensity

    # -----------------------------------------------------------------

    @property
    def old_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @old_stars_luminosity.setter
    def old_stars_luminosity(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.old_stars_parameters.luminosity = value

        # Save
        save_mapping(self.old_stars_parameters_path, self.old_stars_parameters)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_normalization_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter(self.old_stars_parameters.filter)

    # -----------------------------------------------------------------

    @property
    def old_stars_normalization_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.wavelength

    # -----------------------------------------------------------------

    @property
    def old_stars_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.neutral_luminosity

    # -----------------------------------------------------------------

    @property
    def old_stars_fluxdensity(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.old_stars_fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_age(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.age * u("Gyr")

    # -----------------------------------------------------------------

    @property
    def young_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @young_stars_luminosity.setter
    def young_stars_luminosity(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.young_stars_parameters.luminosity = value

        # Save
        save_mapping(self.young_stars_parameters_path, self.young_stars_parameters)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_normalization_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter(self.young_stars_parameters.filter)

    # -----------------------------------------------------------------

    @property
    def young_stars_normalization_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.wavelength

    # -----------------------------------------------------------------

    @property
    def young_stars_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.neutral_luminosity

    # -----------------------------------------------------------------

    @property
    def young_stars_fluxdensity(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_age(self):

        """
        This function ...
        :return:
        """

        return self.young_stars_parameters.age * u("Gyr")

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.luminosity

    # -----------------------------------------------------------------

    @ionizing_stars_luminosity.setter
    def ionizing_stars_luminosity(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.ionizing_stars_parameters.luminosity = value

        # Save
        save_mapping(self.ionizing_stars_parameters_path, self.ionizing_stars_parameters)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_normalization_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter(self.ionizing_stars_parameters.filter)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_normalization_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.wavelength

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.neutral_luminosity

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_fluxdensity(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_stars_parameters.fluxdensity

    # -----------------------------------------------------------------

    @property
    def dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.dust_parameters.mass

    # -----------------------------------------------------------------

    @dust_mass.setter
    def dust_mass(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.dust_parameters.mass = value

        # Save
        save_mapping(self.dust_parameters_path, self.dust_parameters)

    # -----------------------------------------------------------------

    @property
    def old_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.old_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @old_stars_scaleheight.setter
    def old_stars_scaleheight(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.old_stars_parameters.scale_height = value

        # Save
        save_mapping(self.old_stars_parameters_path, self.old_stars_parameters)

    # -----------------------------------------------------------------

    @property
    def young_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.young_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @young_stars_scaleheight.setter
    def young_stars_scaleheight(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new value
        self.young_stars_parameters.scale_height = value

        # Save
        save_mapping(self.young_stars_parameters_path, self.young_stars_parameters)

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.ionizing_stars_parameters.scale_height

    # -----------------------------------------------------------------

    @ionizing_stars_scaleheight.setter
    def ionizing_stars_scaleheight(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the new scale height
        self.ionizing_stars_parameters.scale_height = value

        # Save
        save_mapping(self.ionizing_stars_parameters_path, self.ionizing_stars_parameters)

    # -----------------------------------------------------------------

    @property
    def dust_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.dust_parameters.scale_height

    # -----------------------------------------------------------------

    @dust_scaleheight.setter
    def dust_scaleheight(self, value):

        """
        Thisf unction ...
        :param value:
        :return:
        """

        # Set the new value
        self.dust_parameters.scale_height = value

        # Save
        save_mapping(self.dust_parameters_path, self.dust_parameters)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.old_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.young_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_deprojection_path(self):

        """
        Tihs function ...
        :return: 
        """

        return fs.join(self.ionizing_stars_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_deprojection_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.dust_component_path, deprojection_filename)

    # -----------------------------------------------------------------

    def get_stellar_deprojection_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return fs.join(self.stellar_paths[component_name], deprojection_filename)

    # -----------------------------------------------------------------

    def get_dust_deprojection_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return fs.join(self.dust_paths[component_name], deprojection_filename)

    # -----------------------------------------------------------------

    def has_stellar_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return fs.is_file(self.get_stellar_deprojection_path(component_name))

    # -----------------------------------------------------------------

    def has_dust_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return fs.is_file(self.get_dust_deprojection_path(component_name))

    # -----------------------------------------------------------------

    def get_stellar_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return DeprojectionModel3D.from_file(self.get_stellar_deprojection_path(component_name))

    # -----------------------------------------------------------------

    def get_dust_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        return DeprojectionModel3D.from_file(self.get_dust_deprojection_path(component_name))

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_deprojection.galaxy_distance

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # TODO: REDO THE CONVERSION FROM FLUX TO LUMINOSITIES BASED ON THE NEW DISTANCE?

        # Loop over the stellar components
        for component_name in self.stellar_component_names:
            if not self.has_stellar_deprojection(component_name): continue

            # Load the deprojection
            deprojection = self.get_stellar_deprojection(component_name)

            # Set the new distance and save
            deprojection.galaxy_distance = value
            deprojection.save()

        # Loop over the dust components
        for component_name in self.dust_component_names:
            if not self.has_dust_deprojection(component_name): continue

            # Load the deprojection
            deprojection = self.get_dust_deprojection(component_name)

            # Set the new distance and save
            deprojection.galaxy_distance = value
            deprojection.save()

    # -----------------------------------------------------------------

    @property
    def inclination(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_deprojection.inclination

    # -----------------------------------------------------------------

    @inclination.setter
    def inclination(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Loop over the stellar components
        for component_name in self.stellar_component_names:
            if not self.has_stellar_deprojection(component_name): continue

            # Get the deprojection
            deprojection = self.get_stellar_deprojection(component_name)

            # Set the new value and save
            deprojection.inclination = value
            deprojection.save()

        # Loop over the dust components
        for component_name in self.dust_component_names:
            if not self.has_dust_deprojection(component_name): continue

            # Get the deprojection
            deprojection = self.get_dust_deprojection(component_name)

            # Set the new value and save
            deprojection.inclination = value
            deprojection.save()

    # -----------------------------------------------------------------

    @property
    def position_angle(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_deprojection.position_angle

    # -----------------------------------------------------------------

    @position_angle.setter
    def position_angle(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Loop over the stellar components
        for component_name in self.stellar_component_names:
            if not self.has_stellar_deprojection(component_name): continue

            # Get the deprojection
            deprojection = self.get_stellar_deprojection(component_name)

            # Set the new value and save
            deprojection.position_angle = value
            deprojection.save()

        # Loop over the dust components
        for component_name in self.dust_component_names:
            if not self.has_dust_deprojection(component_name): continue

            # Get the deprojection
            deprojection = self.get_dust_deprojection(component_name)

            # Set the new value and save
            deprojection.position_angle = value
            deprojection.save()

    # -----------------------------------------------------------------

    @property
    def metallicity(self):

        """
        This function ...
        :return:
        """

        return self.old_stars_parameters.metallicity

    # -----------------------------------------------------------------

    @metallicity.setter
    def metallicity(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Loop over the stellar components
        for component_name in self.stellar_component_names:

            # Get the parameters path
            path = self.get_stellar_component_parameters_path(component_name)

            # Load the parameters
            parameters = self.get_stellar_component_parameters(component_name)

            # Set the new metallicity and save
            parameters.metallicity = value
            save_mapping(path, parameters)

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

    def load_bulge_component(self):

        """
        This function ...
        :return:
        """

        from .suite import load_component
        return load_component(self.bulge_component_path)

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

    def get_stellar_component_parameters_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get component path
        path = self.get_stellar_component_path(component_name)

        # Determine the parameters file path
        parameters_path = fs.join(path, parameters_filename)

        # Check
        if fs.is_file(parameters_path): return parameters_path
        else: return None

    # -----------------------------------------------------------------

    def get_dust_component_parameters_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get component path
        path = self.get_dust_component_path(component_name)

        # Determine the parameters file path
        parameters_path = fs.join(path, parameters_filename)

        # Check
        if fs.is_file(parameters_path): return parameters_path
        else: return None

    # -----------------------------------------------------------------

    def get_stellar_component_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the parameters filepath
        path = self.get_stellar_component_parameters_path(name)

        # Return
        if path is None: return None
        else: return open_mapping(path)

    # -----------------------------------------------------------------

    def get_stellar_component_normalization_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        parameters = self.get_stellar_component_parameters(name)
        if parameters is None: return None
        else: return parameters.filter

    # -----------------------------------------------------------------

    def get_stellar_component_normalization_wavelength(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        parameters = self.get_stellar_component_parameters(name)
        if parameters is None: return None
        else: return parameters.wavelength

    # -----------------------------------------------------------------

    def get_dust_component_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the parameters filepath
        path = self.get_dust_component_parameters_path(name)

        # Return
        if path is None: return None
        else: return open_mapping(path)

    # -----------------------------------------------------------------

    def get_stellar_component_deprojection_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get component path
        path = self.get_stellar_component_path(component_name)

        # Determine the deprojection file path
        deprojection_path = fs.join(path, deprojection_filename)

        # Return
        if fs.is_file(deprojection_path): return deprojection_path
        else: return None

    # -----------------------------------------------------------------

    def get_dust_component_deprojection_path(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get component path
        path = self.get_dust_component_path(component_name)

        # Determine the deprojection file path
        deprojection_path = fs.join(path, deprojection_filename)

        # Return
        if fs.is_file(deprojection_path): return deprojection_path
        else: return None

    # -----------------------------------------------------------------

    def get_stellar_component_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get deprojection filepath
        path = self.get_stellar_component_deprojection_path(component_name)

        # Return
        if path is None: return None
        else: return DeprojectionModel3D.from_file(path)

    # -----------------------------------------------------------------

    def get_dust_component_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get deprojection filepath
        path = self.get_dust_component_deprojection_path(component_name)

        # Return
        if path is None: return None
        else: return DeprojectionModel3D.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_filters(self):

        """
        This function ...
        :return:
        """

        # Create a set of the filters
        filters = set()

        # Add the filters
        filters.add(self.bulge_normalization_filter)
        filters.add(self.old_stars_normalization_filter)
        filters.add(self.young_stars_normalization_filter)
        filters.add(self.ionizing_stars_normalization_filter)
        for name in self.additional_stellar_names:
            fltr = self.get_stellar_component_normalization_filter(name)
            if fltr is not None: filters.add(fltr)

        # Return the list of filters
        return list(filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Create a set of the wavelengths
        wavelengths = set()

        # Add the wavelengths
        wavelengths.add(self.bulge_normalization_wavelength)
        wavelengths.add(self.old_stars_normalization_wavelength)
        wavelengths.add(self.young_stars_normalization_wavelength)
        wavelengths.add(self.ionizing_stars_normalization_wavelength)
        for name in self.additional_stellar_names:
            wavelength = self.get_stellar_component_normalization_wavelength(name)
            if wavelength is not None: wavelengths.add(wavelength)

        # Return the list of wavelengths
        return list(wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_center_wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.center for fltr in self.normalization_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_effective_wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.effective for fltr in self.normalization_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_pivot_wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.pivot for fltr in self.normalization_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_mean_wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.mean for fltr in self.normalization_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_peak_wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.peak for fltr in self.normalization_filters]

# -----------------------------------------------------------------
