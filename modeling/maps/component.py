#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.component Contains the MapsComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.tools.colours import get_filters_for_colour
from ...magic.core.frame import Frame
from ...core.tools import types
from ...core.filter.filter import parse_filter
from ...magic.core.list import NamedFrameList
from ...core.tools.serialization import load_dict, write_dict

# -----------------------------------------------------------------

origins_filename = "origins.txt"

# -----------------------------------------------------------------

class MapsComponent(GalaxyModelingComponent):
    
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
        super(MapsComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The maps
        self.maps = dict()
        self.error_maps = dict()

        # The origins
        self.origins = dict()

        # The paths to the maps
        self.paths = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def sub_path_for_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        name = self.sub_name_for_command(command)
        return fs.join(self.maps_path, name)

    # -----------------------------------------------------------------

    def sub_name_for_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        if command == "make_colour_maps": return self.maps_colours_name
        elif command == "make_ssfr_maps": return self.maps_ssfr_name
        elif command == "make_tir_maps": return self.maps_tir_name
        elif command == "make_attenuation_maps": return self.maps_attenuation_path
        elif command == "make_old_stars_map": return self.maps_old_name
        elif command == "make_young_stars_map": return self.maps_young_name
        elif command == "make_ionizing_stars_map": return self.maps_ionizing_name
        elif command == "make_dust_map": return self.maps_dust_name
        else: raise ValueError("Invalid commands: " + command)

    # -----------------------------------------------------------------

    def command_for_sub_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name == self.maps_colours_name: return "make_colour_maps"
        elif name == self.maps_ssfr_name: return "make_ssfr_maps"
        elif name == self.maps_tir_name: return "make_tir_maps"
        elif name == self.maps_attenuation_name: return "make_attenuation_maps"
        elif name == self.maps_old_name: return "make_old_stars_map"
        elif name == self.maps_young_name: return "make_young_stars_map"
        elif name == self.maps_ionizing_name: return "make_ionizing_stars_map"
        elif name == self.maps_dust_name: return "make_dust_map"
        else: raise ValueError("Invalid sub name: " + name)

    # -----------------------------------------------------------------

    def command_for_sub_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        name = fs.relative_to(path, self.maps_path).split("/")[0]
        return self.command_for_sub_name(name)

    # -----------------------------------------------------------------

    @property
    def maps_colours_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_colours_path

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_colours_name

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_ssfr_path

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_ssfr_name

    # -----------------------------------------------------------------

    @property
    def maps_tir_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_tir_path

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):

        """
        THis function ...
        :return:
        """

        return self.environment.maps_tir_name

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_attenuation_path

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_attenuation_name

    # -----------------------------------------------------------------

    @property
    def maps_old_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_old_path

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_old_name

    # -----------------------------------------------------------------

    @property
    def maps_young_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_young_path

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_young_name

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_path(self):

        """
        This fucntion ...
        :return:
        """

        return self.environment.map_ionizing_path

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):

        """
        THis function ...
        :return:
        """

        return self.environment.maps_ionizing_name

    # -----------------------------------------------------------------

    @property
    def maps_dust_path(self):

        """
        This fucntion ...
        :return:
        """

        return self.environment.maps_dust_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):

        """
        THis function ...
        :return:
        """

        return self.environment.maps_dust_name

    # -----------------------------------------------------------------

    @property
    def maps_sub_paths(self):

        """
        This function ...
        :return: 
        """

        return [self.maps_colours_path, self.maps_ssfr_path, self.maps_tir_path, self.maps_attenuation_path, self.maps_old_path, self.maps_young_path, self.maps_ionizing_path, self.maps_dust_path]

    # -----------------------------------------------------------------

    @property
    def maps_sub_names(self):

        """
        This function ...
        :return: 
        """

        return [fs.name(path) for path in self.maps_sub_paths]

    # -----------------------------------------------------------------

    def get_origins_sub_name(self, name, flatten=False):

        """
        This function ...
        :param name:
        :param flatten:
        :return: 
        """

        # Determine path
        sub_path = fs.join(self.maps_path, name)
        if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")
        direct_origins_path = fs.join(sub_path, origins_filename)

        # No subdirectories
        if fs.is_file(direct_origins_path): origins = load_dict(direct_origins_path)

        # Subdirectories
        else:

            # Initialize
            origins = dict()

            # Loop over subdirectories
            for method_path in fs.directories_in_path(sub_path):

                origins_path = fs.join(method_path, origins_filename)
                if not fs.is_file(origins_path): raise ValueError("File '" + origins_path + "' is missing")

                #print(origins_path)

                # Determine method
                method = fs.name(method_path)

                # Load the origins for this method
                origins_method = load_dict(origins_path)

                # Flatten into a one-level dict
                if flatten:
                    for map_name in origins_method: origins[method + "_" + map_name] = origins_method[map_name]

                # Don't flatten: get nested dict
                else: origins[method] = origins_method

        # Return the origins
        return origins

    # -----------------------------------------------------------------

    def get_origins_sub_names(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        origins = dict()
        for name in self.maps_sub_names: origins[name] = self.get_origins_sub_name(name, flatten=flatten)
        return origins

    # -----------------------------------------------------------------

    def get_colours_origins(self, flatten=False):
        
        """
        This function ...
        :param flatten:
        :return: 
        """

        return self.get_origins_sub_name(self.maps_colours_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_ssfr_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_ssfr_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_tir_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_tir_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_attenuation_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_attenuation_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_old_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return: 
        """

        return self.get_origins_sub_name(self.maps_old_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_young_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return: 
        """

        return self.get_origins_sub_name(self.maps_young_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_ionizing_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_ionizing_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_dust_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return: 
        """

        return self.get_origins_sub_name(self.maps_dust_name, flatten=flatten)

    # -----------------------------------------------------------------

    @abstractproperty
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @property
    def maps_sub_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_sub_path)

    # -----------------------------------------------------------------

    @property
    def colour_map_filters(self):

        """
        This function ...
        :return: 
        """

        filters = []

        # Loop over the images in the colour maps directory
        for path, name in fs.files_in_path(self.maps_colours_path, extension="fits", returns=["path", "name"]):

            # Get the filters
            fltr_a, fltr_b = get_filters_for_colour(name)

            # Add a tuple
            filters.append((fltr_a, fltr_b))

        # Return the list
        return filters

    # -----------------------------------------------------------------

    @property
    def colour_map_filters_and_paths(self):

        """
        This function ...
        :return: 
        """

        filters = dict()

        # Loop over the images in the colour maps directory
        for path, name in fs.files_in_path(self.maps_colours_path, extension="fits", returns=["path", "name"]):

            # Get the filters
            fltr_a, fltr_b = get_filters_for_colour(name)

            # Add
            filters[(fltr_a, fltr_b)] = path

        # Return the dictionary
        return filters

    # -----------------------------------------------------------------

    def has_colour_map_for_filters(self, fltr_a, fltr_b):

        """
        This function ...
        :param fltr_a: 
        :param fltr_b: 
        :return: 
        """

        # Loop over the filters of the existing colour maps
        for fltr_ai, fltr_bi in self.colour_map_filters:

            if (fltr_ai, fltr_bi) == (fltr_a, fltr_b): return True
            if (fltr_bi, fltr_ai) == (fltr_a, fltr_b): return True

        # No colour map encountered
        return False

    # -----------------------------------------------------------------

    def get_colour_map_and_name_for_filters(self, fltr_a, fltr_b):

        """
        This function ...
        :param fltr_a:
        :param fltr_b:
        :return:
        """

        # Loop over the existing colour maps
        filters = self.colour_map_filters_and_paths
        for fltr_ai, fltr_bi in filters:

            # Get the path
            path = filters[(fltr_ai, fltr_bi)]

            # Determine the name
            name = fs.strip_extension(fs.name(path))

            # Check colours
            if (fltr_ai, fltr_bi) == (fltr_a, fltr_b): return Frame.from_file(path), name
            if (fltr_bi, fltr_ai) == (fltr_a, fltr_b): return -1. * Frame.from_file(path), name

        # No colour map encountered
        return None, None

    # -----------------------------------------------------------------

    def get_colour_map_for_filters(self, fltr_a, fltr_b):

        """
        This function ...
        :param fltr_a: 
        :param fltr_b: 
        :return: 
        """

        # Loop over the existing colour maps
        filters = self.colour_map_filters_and_paths
        for fltr_ai, fltr_bi in filters:

            # Get the path
            path = filters[(fltr_ai, fltr_bi)]

            # Check colours
            if (fltr_ai, fltr_bi) == (fltr_a, fltr_b): return Frame.from_file(path)
            if (fltr_bi, fltr_ai) == (fltr_a, fltr_b): return -1. * Frame.from_file(path)

        # No colour map encountered
        return None

    # -----------------------------------------------------------------

    def get_colour_map(self, colour):

        """
        This function ...
        :param colour:
        :return:
        """

        fltr_a, fltr_b = get_filters_for_colour(colour)
        return self.get_colour_map_for_filters(fltr_a, fltr_b)

    # -----------------------------------------------------------------

    def get_colour_map_and_name(self, colour):

        """
        This function ...
        :param colour:
        :return:
        """

        fltr_a, fltr_b = get_filters_for_colour(colour)
        return self.get_colour_map_and_name_for_filters(fltr_a, fltr_b)

    # -----------------------------------------------------------------

    def has_frame_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # OLD
        #return self.dataset.has_frame_for_filter(fltr)

        # NEW
        return fltr in self.frame_list.filters

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr):

        """
        THis function ...
        :param fltr: 
        :return: 
        """

        # OLD
        #return self.dataset.get_frame_for_filter(fltr)

        # NEW
        return self.frame_list[fltr]

    # -----------------------------------------------------------------

    def has_errormap_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # OLD
        #return self.dataset.has_errormap_for_filter(fltr)

        # NEW
        return fltr in self.errormap_list.filters

    # -----------------------------------------------------------------

    def get_errormap_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # OLD
        #return self.dataset.get_errormap_for_filter(fltr)

        # NEW
        return self.errormap_list[fltr]

    # -----------------------------------------------------------------

    def get_colour_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_colours_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_ssfr_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return: 
        """

        #return NamedFrameList.from_directory(self.maps_ssfr_path).to_dictionary()

        return self.get_maps_sub_name(self.maps_ssfr_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_tir_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return: 
        """

        return self.get_maps_sub_name(self.maps_tir_name, flatten=flatten, framelist=framelist)

        #single = self.get_tir_single_maps()
        #multi = self.get_tir_multi_maps()

        #maps = dict()
        #maps["single"] = single
        #maps["multi"] = multi

    # -----------------------------------------------------------------

    def get_tir_single_maps(self):

        """
        This function ...
        :return: 
        """

        path = fs.join(self.maps_tir_path, "single")
        return NamedFrameList.from_directory(path).to_dictionary()

    # -----------------------------------------------------------------

    def get_tir_multi_maps(self):

        """
        This function ...
        :return: 
        """

        path = fs.join(self.maps_tir_path, "multi")
        return NamedFrameList.from_directory(path).to_dictionary()

    # -----------------------------------------------------------------

    def get_tir_single_origins(self):

        """
        This function ...
        :return: 
        """

        path = fs.join(self.maps_tir_path, "single", origins_filename)
        return load_dict(path)

    # -----------------------------------------------------------------

    def get_tir_multi_origins(self):

        """
        This function ...
        :return: 
        """

        path = fs.join(self.maps_tir_path, "multi", origins_filename)
        return load_dict(path)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return: 
        """

        cortese = self.get_cortese_fuv_attenuation_maps()
        buat = self.get_buat_fuv_attenuation_maps()

        if flatten:

            maps = dict()
            for name in cortese: maps["cortese_" + name] = cortese[name]
            for name in buat: maps["buat_" + name] = buat[name]
            return maps

        else:

            maps = dict()
            maps["cortese"] = cortese
            maps["buat"] = buat
            return maps

    # -----------------------------------------------------------------

    def get_fuv_attenuation_origins(self):

        """
        This function ...
        :return: 
        """

        cortese = self.get_cortese_fuv_attenuation_origins()
        buat = self.get_buat_fuv_attenuation_origins()

        origins = dict()
        for name in cortese: origins["cortese_" + name] = cortese[name]
        for name in buat: origins["cortese_" + name] = buat[name]
        return origins

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_maps(self):

        """
        This function ...
        :return: 
        """

        buat_path = fs.join(self.maps_attenuation_path, "buat")
        return NamedFrameList.from_directory(buat_path, contains="FUV").to_dictionary()

    # -----------------------------------------------------------------

    def get_cortese_fuv_attenuation_maps(self):

        """
        This function ...
        :return: 
        """

        cortese_path = fs.join(self.maps_attenuation_path, "cortese")
        return NamedFrameList.from_directory(cortese_path).to_dictionary()

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_origins(self):

        """
        This function ...
        :return: 
        """

        buat_path = fs.join(self.maps_attenuation_path, "buat")
        table_path = fs.join(buat_path, origins_filename)
        return load_dict(table_path)

    # -----------------------------------------------------------------

    def get_cortese_fuv_attenuation_origins(self):

        """
        This function ...
        :return: 
        """

        cortese_path = fs.join(self.maps_attenuation_path, "cortese")
        table_path = fs.join(cortese_path, origins_filename)
        return load_dict(table_path)

    # -----------------------------------------------------------------

    def get_old_stellar_disk_map(self, fltr):

        """
        This fucntion ...
        :param fltr:
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        path = fs.join(self.maps_old_path, "disk", str(fltr) + ".fits")
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    def get_hot_dust_maps(self):

        """
        This function ...
        :return: 
        """

        hot_dust_path = fs.join(self.maps_dust_path, "hot")
        return NamedFrameList.from_directory(hot_dust_path).to_dictionary()

    # -----------------------------------------------------------------

    def get_hot_dust_origins(self):

        """
        This function ...
        :return: 
        """

        hot_dust_path = fs.join(self.maps_dust_path, "hot")
        origins_path = fs.join(hot_dust_path, origins_filename)
        return load_dict(origins_path)

    # -----------------------------------------------------------------

    def get_current_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_sub_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_map_paths_sub_name(self, name, flatten=False):

        """
        This function ...
        :param name:
        :param flatten:
        :return:
        """

        # Determine path
        sub_path = fs.join(self.maps_path, name)
        if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")
        #direct_origins_path = fs.join(sub_path, origins_filename)
        # No subdirectories
        #if fs.is_file(direct_origins_path): origins = load_dict(direct_origins_path)

        # Subdirectories
        if fs.contains_directories(sub_path):

            paths = dict()

            # Loop over the subdirectories
            for method_path, method in fs.directories_in_path(sub_path, returns=["path", "name"]):

                # Get dictionary of file paths, but only FITS files
                files = fs.files_in_path(method_path, returns="dict", extension="fits")

                # Set the map paths, as a dictionary with the filename as keys
                if flatten: paths[method] = files
                else:
                    for map_name in files: paths[method + "_" + map_name] = files[map_name]

            # Return the paths
            return paths

        # Files present
        elif fs.contains_files(sub_path): return fs.files_in_path(sub_path, returns="dict", extension="fits")

        # Nothing present
        else: return dict()

    # -----------------------------------------------------------------

    def get_current_map_paths(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_sub_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_maps_sub_name(self, name, flatten=False, framelist=False):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :return:
        """

        # Initialize the maps dictionary
        maps = dict()

        # Get map paths
        paths = self.get_map_paths_sub_name(name, flatten=flatten)

        # Loop over the entries
        for method_or_name in paths:

            # Methods
            if types.is_dictionary(paths[method_or_name]):

                method = method_or_name
                maps[method] = dict()

                # Loop over the paths, load the maps and add to dictionary
                for name in paths[method_or_name]:

                    map_path = paths[method_or_name][name]
                    try: maps[method][name] = Frame.from_file(map_path)
                    except IOError:
                        command = self.command_for_sub_name(name)
                        log.warning("The " + method + "/" + name + " map is probably damaged. Run the '" + command + "' command again.")
                        log.warning("Removing the " + map_path + " map ...")
                        fs.remove_file(map_path)
                        self.history.remove_entries_and_save(command)

            # Just maps
            elif types.is_string_type(paths[method_or_name]):

                name = method_or_name
                map_path = paths[method_or_name]
                try: maps[name] = Frame.from_file(map_path)
                except IOError:
                    command = self.command_for_sub_name(name)
                    log.warning("The " + name + " map is probably damaged. Run the '" + command + "' command again.")
                    log.warning("Removing the " + map_path + " map ...")
                    fs.remove_file(map_path)
                    self.history.remove_entries_and_save(command)

            # Something wrong
            else: raise RuntimeError("Something went wrong")

        # Return the maps
        if framelist: return NamedFrameList(**maps)
        else: return maps

    # -----------------------------------------------------------------

    def get_current_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_sub_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    @lazyproperty
    def current_maps(self):

        """
        This function is a memoized property implementation of the get_current_maps method
        :return:
        """

        return self.get_current_maps()

    # -----------------------------------------------------------------

    def get_path_for_map(self, name, method=None):

        """
        This function ...
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Determine path
            map_path = fs.join(path, name + ".fits")

        # Determine path
        else: map_path = fs.join(self.maps_sub_path, name + ".fits")

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Create directory
                #path = fs.create_directory_in(self.maps_sub_path, method)

                # Loop over the maps
                for name in self.maps[method]:

                    # Determine path
                    #map_path = fs.join(path, name + ".fits")

                    map_path = self.get_path_for_map(name, method)

                    if fs.is_file(map_path): continue

                    # Save
                    self.maps[method][name].saveto(map_path)

            # No different methods
            else:

                # Determine path
                #map_path = fs.join(self.maps_sub_path, method + ".fits")

                map_path = self.get_path_for_map(method)

                if fs.is_file(map_path): continue

                # Save
                self.maps[method].saveto(map_path)

    # -----------------------------------------------------------------

    #def write_error_maps(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Inform the user
        #log.info("Writing the error maps (with different methods) ...")

        # Loop over the methods
        #for method in self.maps:

            # Create a directory
            #path = fs.create_directory_in(self.maps_dust_path, method)

            # Loop over the maps
            #for name in self.error_maps[method]:

                # Determine path
                #map_path = fs.join(path, name + "_error.fits")

                # Save the map
                #self.maps[method][name].saveto(map_path)

    # -----------------------------------------------------------------

    def write_origins(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the methods
        for method in self.origins:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Directory path
                path = fs.join(self.maps_sub_path, method)

                # Origins path
                origins_path = fs.join(path, origins_filename)

                # Write
                write_dict(self.origins[method], origins_path)

            # No different methods
            else:

                # Determine the origins file path
                origins_path = fs.join(self.maps_sub_path, origins_filename)

                # Write
                write_dict(self.origins, origins_path)

# -----------------------------------------------------------------

def get_dust_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "dust")

# -----------------------------------------------------------------

def get_dust_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_dust_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_old_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "old")

# -----------------------------------------------------------------

def get_old_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_old_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_young_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "young")

# -----------------------------------------------------------------

def get_young_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_young_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_ionizing_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "ionizing")

# -----------------------------------------------------------------

def get_ionizing_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_ionizing_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------
