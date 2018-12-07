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
from collections import defaultdict
from abc import ABCMeta, abstractproperty, abstractmethod

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.tools.colours import get_filters_for_colour
from ...magic.core.frame import Frame
from ...core.tools import types
from ...core.filter.filter import parse_filter
from ...core.tools.serialization import write_dict
from ...core.tools import sequences
from ...core.basics.configuration import prompt_proceed
from ...core.tools.stringify import tostr
from ...core.launch.pts import find_match
from ...core.tools import introspection
from .collection import MapsCollection, StaticMapsCollection
from .selection import ComponentMapsSelection, StaticComponentMapsSelection
from ...core.tools.utils import lazyproperty
from ..core.environment import colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ...core.basics.configuration import prompt_string_list
from ...core.basics.containers import create_subdict
from ...magic.tools import plotting
from ...magic.core.image import Image
from ...core.tools import strings

# -----------------------------------------------------------------

# For counting the number of negative pixels
default_central_ellipse_factor = 0.4

# -----------------------------------------------------------------

origins_filename = "origins.txt"
methods_filename = "methods.txt"

# -----------------------------------------------------------------

maps_commands = ["make_colours_maps", "make_ssfr_maps", "make_tir_maps", "make_attenuation_maps", "make_old_stellar_maps", "make_dust_map", "make_young_stellar_maps", "make_ionizing_stellar_maps"]

# -----------------------------------------------------------------

titles = dict()
titles[colours_name] = "Colour"
titles[ssfr_name] = "sSFR"
titles[tir_name] = "TIR"
titles[attenuation_name] = "UV dust attenuation"
titles[old_name] = "Old stars"
titles[young_name] = "Young stars"
titles[ionizing_name] = "Ionizing stars"
titles[dust_name] = "Dust"

# -----------------------------------------------------------------

def maps_commands_before(command):

    """
    This function ...
    :param command:
    :return:
    """

    # FInd index of command
    index = sequences.find_unique(maps_commands, command)

    # Return
    if index == 0: return []
    else: return maps_commands[:index]

# -----------------------------------------------------------------

def maps_commands_after(command):

    """
    This function ...
    :param command:
    :return:
    """

    # Find index of command
    index = sequences.find_unique(maps_commands, command)

    # Return
    if index == len(maps_commands) - 1: return []
    else: return maps_commands[index+1:]

# -----------------------------------------------------------------

def maps_commands_before_and_including(command):

    """
    This function ...
    :param command:
    :return:
    """

    return maps_commands_before(command) + [command]

# -----------------------------------------------------------------

def maps_commands_after_and_including(command):

    """
    This function ...
    :param command:
    :return:
    """

    return [command] + maps_commands_after(command)

# -----------------------------------------------------------------

class MapMakerBase(GalaxyModelingComponent):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MapMakerBase, self).__init__(*args, **kwargs)

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # Extra maps
        self.extra_maps = dict()

        # Name for extra maps directory
        self.extra_maps_name = "extra"

    # -----------------------------------------------------------------

    @abstractproperty
    def maps_sub_path(self):
        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def current_maps(self):

        """
        This function is a memoized property implementation of the get_current_maps method
        :return:
        """

        return self.get_current_maps()

    # -----------------------------------------------------------------

    def get_current_maps_method(self, method):

        """
        This
        :param method:
        :return:
        """

        maps = self.get_current_maps()
        if maps is None: return None

        if method not in maps: return dict()
        else: return maps[method]

    # -----------------------------------------------------------------

    def get_current_extra_maps_method(self, method):

        """
        Thisf unction ...
        :param method:
        :return:
        """

        maps = self.get_current_extra_maps()
        if maps is None: return None

        if method not in maps: return dict()
        else: return maps[method]

    # -----------------------------------------------------------------

    @property
    def maps_sub_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_sub_path)

    # -----------------------------------------------------------------

    def get_maps_sub_name(self, name, flatten=False, framelist=False, method=None, factors=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :param factors:
        :return:
        """

        return self.collection.get_maps_sub_name(name, flatten=flatten, framelist=framelist, method=method, factors=factors)

    # -----------------------------------------------------------------

    def get_extra_maps_sub_name(self, name, flatten=False, framelist=False, method=None):

        """
        Thisf unction ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.collection.get_extra_maps_sub_name(name, flatten=flatten, framelist=framelist, method=method)

    # -----------------------------------------------------------------

    def get_current_maps(self, flatten=False, framelist=False, method=None, factors=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param method:
        :param factors:
        :return:
        """

        return self.get_maps_sub_name(self.maps_sub_name, flatten=flatten, framelist=framelist, method=method, factors=factors)

    # -----------------------------------------------------------------

    def get_current_extra_maps(self, flatten=False, framelist=False, method=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.get_extra_maps_sub_name(self.maps_sub_name, flatten=flatten, framelist=framelist, method=method)

    # -----------------------------------------------------------------

    def get_path_for_map(self, name, method=None, add_extension=True, extension="fits"):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Determine path
            if add_extension: map_path = fs.join(path, name + "." + extension)
            else: map_path = fs.join(path, name)

        # Determine path
        else:

            if add_extension: map_path = fs.join(self.maps_sub_path, name + "." + extension)
            else: map_path = fs.join(self.maps_sub_path, name)

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_colour_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        Thisf unction ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_colour_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_colour_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_colour_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_colour_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_colour_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def select_ssfr_maps(self, names=None, method=None, methods=None, prompt=False, title="sSFR maps"):

        """
        This function ...
        :param names:
        :param method:
        :param methods:
        :param prompt:
        :param title:
        :return:
        """

        # Get maps
        ssfrs = self.get_ssfr_maps(flatten=True, method=method, methods=methods)

        # Origins
        ssfrs_origins = self.get_ssfr_origins(flatten=True, method=method, methods=methods)

        # Methods
        ssfrs_methods = self.get_ssfr_methods(flatten=True, method=method, methods=methods)

        # Nans
        ssfrs_nans = self.get_ssfr_nans(flatten=True, method=method, methods=methods)

        # Get only certain sSFR maps
        if names is not None:
            ssfrs = create_subdict(ssfrs, names)
            ssfrs_origins = create_subdict(ssfrs_origins, names)
            ssfrs_methods = create_subdict(ssfrs_methods, names)
            ssfrs_nans = create_subdict(ssfrs_nans, names)

        # Select interactively
        if prompt:
            ssfrs, ssfr_names = select_maps(ssfrs, title, return_names=True)
            ssfrs_origins = create_subdict(ssfrs_origins, ssfr_names)
            ssfrs_methods = create_subdict(ssfrs_methods, ssfr_names)
            ssfrs_nans = create_subdict(ssfrs_nans, ssfr_names)

        # Return
        return ssfrs, ssfrs_origins, ssfrs_methods, ssfrs_nans

    # -----------------------------------------------------------------

    def get_ssfr_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_ssfr_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_ssfr_maps(self, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        return self.collection.get_ssfr_maps(flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_ssfr_nans(self, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        Thisfunction ...
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        return self.collection.get_ssfr_nans(flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def select_tir_maps(self, names=None, methods=None, prompt=False, title="TIR maps"):

        """
        This function ...
        :param names:
        :param methods:
        :param prompt:
        :param title:
        :return:
        """

        # Get maps
        tirs = self.get_tir_maps(flatten=True, methods=methods)

        # Origins
        tirs_origins = self.get_tir_origins(flatten=True, methods=methods)

        # Methods
        tirs_methods = self.get_tir_methods(flatten=True, methods=methods)

        # Nans
        tirs_nans = self.get_tir_nans(flatten=True, methods=methods)

        # Get only certain TIR maps
        if names is not None:
            tirs = create_subdict(tirs, names)
            tirs_origins = create_subdict(tirs_origins, names)
            tirs_methods = create_subdict(tirs_methods, names)
            tirs_nans = create_subdict(tirs_nans, names)

        # Select interactively
        if prompt:
            tirs, tir_names = select_maps(tirs, title, return_names=True)
            tirs_origins = create_subdict(tirs_origins, tir_names)
            tirs_methods = create_subdict(tirs_methods, tir_names)
            tirs_nans = create_subdict(tirs_nans, tir_names)

        # Return
        return tirs, tirs_origins, tirs_methods, tirs_nans

    # -----------------------------------------------------------------

    def get_tir_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_tir_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_tir_maps(self, flatten=False, framelist=False, methods=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param methods:
        :return:
        """

        return self.collection.get_tir_maps(flatten=flatten, framelist=framelist, methods=methods)

    # -----------------------------------------------------------------

    def get_tir_nans(self, flatten=False, framelist=False, methods=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param methods:
        :return:
        """

        return self.collection.get_tir_nans(flatten=flatten, framelist=framelist, methods=methods)

    # -----------------------------------------------------------------

    def get_attenuation_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_attenuation_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_attenuation_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist::
        :return:
        """

        return self.collection.get_attenuation_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_attenuation_extra_maps(self, flatten=False, framelist=False):

        """
        Thisf ucntion ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_attenuation_extra_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_attenuation_nans(self, flatten=False, framelist=False):

        """
        Thisfunction ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_attenuation_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        Thisf unction ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_old_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_old_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_old_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_nans(self, flatten=False, framelist=False):

        """
        Thisnfunction ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_old_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_young_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_young_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_young_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_young_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_young_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_young_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_ionizing_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_ionizing_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_ionizing_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_ionizing_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_ionizing_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_ionizing_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_dust_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        Thisf ucntion ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.collection.get_dust_map_paths(flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_dust_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_dust_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_dust_negatives(self, flatten=False, framelist=False, method=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.collection.get_dust_negatives(flatten=flatten, framelist=framelist, method=method)

    # -----------------------------------------------------------------

    def get_dust_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_dust_nans(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_stellar_disk_maps(self, framelist=False):

        """
        This function ...
        :param framelist:
        :return:
        """

        return self.collection.get_old_stellar_disk_maps(framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_stellar_total_origins(self):

        """
        Thisfunction ...
        :return:
        """

        return self.collection.get_old_stellar_total_origins()

    # -----------------------------------------------------------------

    def get_old_stellar_disk_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_disk_origins()

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_bulge_origins()

    # -----------------------------------------------------------------

    def get_old_stellar_total_methods(self):

        """
        Thisn function ...
        :return:
        """

        return self.collection.get_old_stellar_total_methods()

    # -----------------------------------------------------------------

    def get_old_stellar_disk_methods(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_disk_methods()

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_methods(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_bulge_methods()

    # -----------------------------------------------------------------

    def get_old_stellar_total_map(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.collection.get_old_stellar_total_map(fltr)

    # -----------------------------------------------------------------

    def get_old_stellar_disk_map(self, fltr):

        """
        This fucntion ...
        :param fltr:
        :return:
        """

        return self.collection.get_old_stellar_disk_map(fltr)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        COPY FROM FUNCTION IN MAPSCOMPONENT
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Loop over the maps
                for name in self.maps[method]:

                    # Determine path
                    map_path = self.get_path_for_map(name, method)

                    # If map already exists and we don't have to remake
                    if fs.is_file(map_path) and not self.config.remake: continue

                    # Save
                    self.maps[method][name].saveto(map_path)

            # No different methods
            else:

                # Determine path
                map_path = self.get_path_for_map(method)

                # If map already exists and we don't have to remake
                if fs.is_file(map_path) and not self.config.remake: continue

                # Save
                self.maps[method].saveto(map_path)

    # -----------------------------------------------------------------

    def get_path_for_extra_map(self, name, method=None, add_extension=True, extension="fits"):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create extra path
            extra_path = fs.create_directory_in(path, self.extra_maps_name)

            # Determine path
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Determine path
        else:

            extra_path = fs.create_directory_in(self.maps_sub_path, self.extra_maps_name)
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def write_extra_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the extra maps ...")

        # Loop over the methods
        for method in self.extra_maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.extra_maps[method]):

                # Loop over the maps
                for name in self.extra_maps[method]:

                    # Determine path
                    map_path = self.get_path_for_extra_map(name, method)

                    # If map already exists and we don't have to remake
                    if fs.is_file(map_path) and not self.config.remake: continue

                    # Save
                    self.extra_maps[method][name].saveto(map_path)

            # No different methods
            else:

                # Determine path
                map_path = self.get_path_for_extra_map(method)

                # If map already exists and we don't have to remake
                if fs.is_file(map_path) and not self.config.remake: continue

                # Save
                self.extra_maps[method].saveto(map_path)

    # -----------------------------------------------------------------

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]): return True

            # No subdivision: method is actually the map name
            else: return False

        # We shouldn't get here
        raise ValueError("No maps present yet")

    # -----------------------------------------------------------------

    def write_origins(self):

        """
        COPY FROM FUNCTION IN MAPSCOMPONENT
        :return:
        """

        # Inform the user
        log.info("Writing the map origins ...")

        # CHECK WHETHER ORIGIN IS DEFINED FOR EACH MAP
        self.check_origins()

        # If has different methods
        if self.has_methods:

            # Loop over the methods
            for method in self.origins:

                # Depending on whether subdictionaries
                # if types.is_dictionary(self.origins[method]):

                # Directory path
                path = fs.create_directory_in(self.maps_sub_path, method)

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

    def write_methods(self):

        """
        COPY FROM THE FUNCTION IN MAPSCOMPONENT
        :return:
        """

        # Inform the user
        log.info("Writing the map methods ...")

        # CHECK WHETHER METHOD IS DEFINED FOR EACH MAP
        self.check_methods()

        # If has different methods
        if self.has_methods:

            # Loop over the methods OF THIS STEP
            for method in self.methods:

                # Depending on whether subdirectories
                # if types.is_dictionary(self.methods[method]):

                # Directory path
                path = fs.create_directory_in(self.maps_sub_path, method)

                # Methods path
                methods_path = fs.join(path, methods_filename)

                # Write
                write_dict(self.methods[method], methods_path)

        # No different methods
        else:

            # Determine the file path
            methods_path = fs.join(self.maps_sub_path, methods_filename)

            # Write
            write_dict(self.methods, methods_path)

    # -----------------------------------------------------------------
    # PLOTTING
    # -----------------------------------------------------------------

    @property
    def map_plots_name(self):

        """
        This function ...
        :return:
        """

        return "plots"

    # -----------------------------------------------------------------

    @property
    def negatives_plots_name(self):

        """
        This function ...
        :return:
        """

        return "negatives"

    # -----------------------------------------------------------------

    @property
    def nans_plots_name(self):

        """
        This function ...
        :return:
        """

        return "nans"

    # -----------------------------------------------------------------

    @property
    def contour_plots_name(self):

        """
        This function ...
        :return:
        """

        return "contours"

    # -----------------------------------------------------------------

    @property
    def profile_plots_name(self):

        """
        Thisfunction ...
        :return:
        """

        return "profiles"

    # -----------------------------------------------------------------

    def get_path_for_map_plot(self, name, method=None, add_extension=True, extension="pdf", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create path
            plots_path = fs.create_directory_in(path, self.map_plots_name)

            # Set base path
            map_base_path = fs.join(plots_path, name)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plots_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Determine path
        else:

            plots_path = fs.create_directory_in(self.maps_sub_path, self.map_plots_name)

            # Set base path
            map_base_path = fs.join(plots_path, name)

            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plots_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_negatives_plot(self, name, method=None, add_extension=True, extension="pdf", suffix="", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param suffix:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create plot path
            plot_path = fs.create_directory_in(path, self.negatives_plots_name)

            # Set base path
            map_base_path = fs.join(plot_path, name + suffix)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plot_path, exact_name=name+suffix, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Determine path
        else:

            # Create plot path
            plot_path = fs.create_directory_in(self.maps_sub_path, self.negatives_plots_name)

            # Set base path
            map_base_path = fs.join(plot_path, name + suffix)

            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plot_path, exact_name=name+suffix, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_nnegatives_curves_plot(self, name, method=None, extension="pdf"):

        """
        This function ...
        :param name:
        :param method:
        :param extension:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create plot path
            plot_path = fs.create_directory_in(path, self.negatives_plots_name)

            # Set path
            map_base_path = fs.join(plot_path, name + "__counts")
            filepath = map_base_path + "." + extension

        # Determine path
        else:

            # Create plot path
            plot_path = fs.create_directory_in(self.maps_sub_path, self.negatives_plots_name)

            # Set path
            map_base_path = fs.join(plot_path, name + "__counts")
            filepath = map_base_path + "." + extension

        # Return the path
        return filepath

    # -----------------------------------------------------------------

    def get_path_for_nans_plot(self, name, method=None, add_extension=True, extension="pdf", suffix="", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param suffix:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create plot path
            plot_path = fs.create_directory_in(path, self.nans_plots_name)

            # Set base path
            map_base_path = fs.join(plot_path, name + suffix)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plot_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Determine path
        else:

            plot_path = fs.create_directory_in(self.maps_sub_path, self.nans_plots_name)

            # Set base path
            map_base_path = fs.join(plot_path, name + suffix)

            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(plot_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_contour_plot(self, name, method=None, add_extension=True, extension="pdf", suffix="", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param suffix:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create path
            contours_path = fs.create_directory_in(path, self.contour_plots_name)

            # Set base path
            map_base_path = fs.join(contours_path, name + suffix)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(contours_path, exact_name=name+suffix, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Determine path
        else:

            contours_path = fs.create_directory_in(self.maps_sub_path, self.contour_plots_name)

            # Set base path
            map_base_path = fs.join(contours_path, name + suffix)

            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(contours_path, exact_name=name+suffix, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_profile_plot(self, name, method=None, add_extension=True, extension="pdf", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create path
            profiles_path = fs.create_directory_in(path, self.profile_plots_name)

            # Set base path
            map_base_path = fs.join(profiles_path, name)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(profiles_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Determine path
        else:

            profiles_path = fs.create_directory_in(self.maps_sub_path, self.profile_plots_name)

            # Set base path
            map_base_path = fs.join(profiles_path, name)

            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(profiles_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                map_path = map_base_path + "." + extension
            else: map_path = map_base_path

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def plot_maps(self, cmap="viridis", scale="log", format="pdf", cropping_factor=1.3, scales=None, share_limits=True,
                  mask_negatives=False, clear_other_formats=True, show_axes=False, transparent=True, interval=None,
                  strict_vmin=None, strict_vmax=None, soft_vmin=None, soft_vmax=None, cmaps=None, methods=None,
                  not_methods=None, soft_max_scaling=1.5):

        """
        Thisfunction ...
        :param cmap:
        :param scale:
        :param format:
        :param cropping_factor:
        :param scales:
        :param share_limits:
        :param mask_negatives:
        :param clear_other_formats:
        :param show_axes:
        :param transparent:
        :param interval:
        :param strict_vmin:
        :param strict_vmax:
        :param soft_vmin:
        :param soft_vmax:
        :param cmaps:
        :param methods:
        :param not_methods:
        :param soft_max_scaling:
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # Determine min, max etc.
        vmin, vmax, passed_min, passed_max, soft_min, soft_max = set_vminvmax(interval=interval, strict_vmin=strict_vmin, strict_vmax=strict_vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax, share_limits=share_limits)

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Plot method
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting maps for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting the '" + name + "' map ...")

                    # Determine path
                    plot_path = self.get_path_for_map_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If the plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    # Determine scale for this image
                    if scales is not None and method in scales:
                        if isinstance(scales[method], dict):
                            if name in scales[method]: frame_scale = scales[method][name]
                            else: frame_scale = scale
                        else: frame_scale = scales[method]
                    else: frame_scale = scale

                    # Determine color map for this image
                    if cmaps is not None and method in cmaps:
                        if isinstance(cmaps[method], dict):
                            if name in cmaps[method]: frame_cmap = cmaps[method][name]
                            else: frame_cmap = cmap
                        else: frame_cmap = cmaps[method]
                    else: frame_cmap = cmap

                    # Plot
                    minmax = get_minmax(interval=interval, vmin=vmin, vmax=vmax, passed_min=passed_min, passed_max=passed_max, share_limits=share_limits, default_interval="pts")
                    frame = self.maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary
                    vmini, vmaxi = plotting.plot_frame(frame, crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                     truncate_outside=self.truncation_ellipse, path=plot_path, format=format, interval=minmax,
                                                     scale=frame_scale, cmap=frame_cmap, normalize_in=self.truncation_ellipse, colorbar=True,
                                                     mask_negatives=mask_negatives, show_axes=show_axes, transparent=transparent,
                                                     soft_min=soft_min, soft_max=soft_max, soft_max_scaling=soft_max_scaling)
                    if share_limits: vmin, vmax = vmini, vmaxi

                # End of method: reset vmin and vmax
                vmin, vmax, soft_min, soft_max = reset_vminvmax(soft_min, soft_max, interval=interval, strict_vmin=strict_vmin, strict_vmax=strict_vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax, share_limits=share_limits)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting the '" + method + "' map ...")

                # Determine path
                plot_path = self.get_path_for_map_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Determine scale for this image
                if scales is not None and method in scales: frame_scale = scales[method]
                else: frame_scale = scale

                # Determine cmap for this image
                if cmaps is not None and method in cmaps: frame_cmap = cmaps[method]
                else: frame_cmap = cmap

                # Plot
                minmax = get_minmax(interval=interval, vmin=vmin, vmax=vmax, passed_min=passed_min, passed_max=passed_max, share_limits=share_limits, default_interval="pts")
                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary
                vmini, vmaxi = plotting.plot_frame(frame, crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                 truncate_outside=self.truncation_ellipse, path=plot_path, format=format,
                                                 interval=minmax, scale=frame_scale, cmap=frame_cmap, normalize_in=self.truncation_ellipse,
                                                 colorbar=True, mask_negatives=mask_negatives, show_axes=show_axes,
                                                 transparent=transparent, soft_min=soft_min, soft_max=soft_max, soft_max_scaling=soft_max_scaling)
                if share_limits: vmin, vmax = vmini, vmaxi

    # -----------------------------------------------------------------

    @lazyproperty
    def central_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * default_central_ellipse_factor

    # -----------------------------------------------------------------

    def plot_negatives(self, format="pdf", cropping_factor=1.3, clear_other_formats=True, show_axes=False,
                       transparent=True, methods=None, not_methods=None, curves=True, count_within=None):

        """
        This function ...
        :param format:
        :param cropping_factor:
        :param clear_other_formats:
        :param show_axes:
        :param transparent:
        :param methods:
        :param not_methods:
        :param curves:
        :param count_within:
        :return:
        """

        from ...magic.core.mask import Mask
        from ...magic.basics.mask import Mask as oldMask

        # Inform the user
        log.info("Plotting the negative pixel masks ...")

        # Initialize dictionary
        nnegatives = dict()

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Do method?
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting negative pixel masks for the '" + method + "' method ...")

                # Initialize dictionary
                nnegatives_method = dict()

                # Loop over the maps
                for name in self.maps[method]:

                    # Is there a mask?
                    if not (isinstance(self.maps[method][name], Image) and self.maps[method][name].has_mask("negatives")): continue
                    mask = self.maps[method][name].masks["negatives"]

                    # Make sure it's a new mask object
                    if isinstance(mask, oldMask): mask = Mask(mask, wcs=self.maps[method][name].wcs)
                    elif mask.wcs is None: mask.wcs = self.maps[method][name].wcs

                    # Count the number of negatives
                    if count_within is not None: relative_nnegatives = mask.relative_nmasked_in(count_within)
                    else: relative_nnegatives = mask.relative_nmasked
                    if relative_nnegatives > 0.2: log.warning("The number of negative values in the '" + name + "' map is higher than 20%")
                    nnegatives_method[name] = relative_nnegatives

                    # Debugging
                    log.debug("Plotting the negative pixel mask of the '" + name + "' map ...")

                    # Determine path
                    plot_path = self.get_path_for_negatives_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If the plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    # Plot the mask
                    plotting.plot_mask(mask, crop_to=self.truncation_box,
                                     cropping_factor=cropping_factor,
                                     truncate_outside=self.truncation_ellipse, path=plot_path,
                                     format=format, show_axes=show_axes, transparent=transparent)

                # Set the nnegatives dictionary for this method
                nnegatives[method] = nnegatives_method

            # No different methods
            else:

                # Is there a mask?
                if not (isinstance(self.maps[method], Image) and self.maps[method].has_mask("negatives")): continue
                mask = self.maps[method].masks["negatives"]

                # Make sure it's a new mask object
                if isinstance(mask, oldMask): mask = Mask(mask, wcs=self.maps[method].wcs)
                elif mask.wcs is None: mask.wcs = self.maps[method].wcs

                # Count the number of negatives
                if count_within is not None: relative_nnegatives = mask.relative_nmasked_in(count_within)
                else: relative_nnegatives = mask.relative_nmasked
                if relative_nnegatives > 0.2: log.warning("The number of negative values in the '" + method + "' map is higher than 20%")
                nnegatives[method] = relative_nnegatives

                # Debugging
                log.debug("Plotting the negative pixel mask of the '" + method + "' map ...")

                # Determine path
                plot_path = self.get_path_for_negatives_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Plot the mask
                plotting.plot_mask(mask, crop_to=self.truncation_box,
                                   cropping_factor=cropping_factor,
                                   truncate_outside=self.truncation_ellipse, path=plot_path,
                                   format=format, show_axes=show_axes, transparent=transparent)

        from ...core.tools import numbers
        from ...core.basics import containers
        from ...core.tools.parsing import real

        # Subdictionaries for different methods
        if types.is_dictionary_of_dictionaries(nnegatives):

            # Loop over the methods
            for method in nnegatives:

                # Initialize
                grouped = defaultdict(dict)

                # Loop over the maps
                for name in nnegatives[method]:

                    # Split
                    if "__" not in name: continue
                    base_name, last_part = strings.split_at_last(name, "__")

                    # Get factor
                    if not numbers.is_number(last_part): continue
                    factor = real(last_part)

                    # Add
                    grouped[base_name][factor] = nnegatives[method][name]

                # Sort for each base name based on the factors values
                for base_name in grouped: grouped[base_name] = containers.ordered_by_key(grouped[base_name])

                # Loop over the base names
                for base_name in grouped:

                    # Skip if only one datapoint
                    if len(grouped[base_name]) == 1: continue

                    # Determine plot path
                    plot_path = self.get_path_for_nnegatives_curves_plot(base_name, method, extension=format)

                    # Plot
                    plotting.plot_xy(grouped[base_name].keys(), grouped[base_name].values(), path=plot_path)

        # No subdictionaries
        else:

            # Initialize
            grouped = defaultdict(dict)

            # Loop over the maps
            for name in nnegatives:

                # Split
                if "__" not in name: continue
                base_name, last_part = strings.split_at_last(name, "__")

                # Get factor
                if not numbers.is_number(last_part): continue
                factor = real(last_part)

                # Add
                grouped[base_name][factor] = nnegatives[name]

            # Sort for each base name based on the factors values
            for base_name in grouped: grouped[base_name] = containers.ordered_by_key(grouped[base_name])

            # Loop over the base names
            for base_name in grouped:

                # Skip if only one datapoint
                if len(grouped[base_name]) == 1: continue

                # Determine plot path
                plot_path = self.get_path_for_nnegatives_curves_plot(base_name, extension=format)

                # Plot
                plotting.plot_xy(grouped[base_name].keys(), grouped[base_name].values(), path=plot_path)

    # -----------------------------------------------------------------

    def plot_nans(self, format="pdf", cropping_factor=1.3, clear_other_formats=True, show_axes=False, transparent=True,
                  methods=None, not_methods=None):

        """
        This finction ...
        :param format:
        :param cropping_factor:
        :param clear_other_formats:
        :param show_axes:
        :param transparent:
        :param methods:
        :param not_methods:
        :return:
        """

        from ...magic.core.mask import Mask
        from ...magic.basics.mask import Mask as oldMask

        # Inform the user
        log.info("Plotting the NaN pixel masks ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Plot method?
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting NaN pixel masks for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Is there a mask?
                    if not (isinstance(self.maps[method][name], Image) and self.maps[method][name].has_mask("nans")): continue
                    mask = self.maps[method][name].masks["nans"]

                    # Make sure it's a new mask object
                    if isinstance(mask, oldMask): mask = Mask(mask, wcs=self.maps[method][name].wcs)
                    elif mask.wcs is None: mask.wcs = self.maps[method][name].wcs

                    # Debugging
                    log.debug("Plotting the NaN pixel mask of the '" + name + "' map ...")

                    # Determine path
                    plot_path = self.get_path_for_nans_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If the plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    # Plot the mask
                    plotting.plot_mask(mask, crop_to=self.truncation_box,
                                       cropping_factor=cropping_factor,
                                       truncate_outside=self.truncation_ellipse, path=plot_path,
                                       format=format, show_axes=show_axes, transparent=transparent)

            # No different methods
            else:

                # Is there a mask?
                if not (isinstance(self.maps[method], Image) and self.maps[method].has_mask("nans")): continue
                mask = self.maps[method].masks["nans"]

                # Make sure it's a new mask object
                if isinstance(mask, oldMask): mask = Mask(mask, wcs=self.maps[method].wcs)
                elif mask.wcs is None: mask.wcs = self.maps[method].wcs

                # Debugging
                log.debug("Plotting the NaN pixel mask of the '" + method + "' map ...")

                # Determine path
                plot_path = self.get_path_for_nans_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Plot the mask
                plotting.plot_mask(mask, crop_to=self.truncation_box,
                                   cropping_factor=cropping_factor,
                                   truncate_outside=self.truncation_ellipse, path=plot_path,
                                   format=format, show_axes=show_axes, transparent=transparent)

    # -----------------------------------------------------------------

    def get_path_for_extra_map_plot(self, name, method=None, add_extension=True, extension="pdf", clear_other_formats=False):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param clear_other_formats:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create extra path
            extra_path = fs.create_directory_in(path, self.extra_maps_name)

            # Create plot path
            extra_plot_path = fs.create_directory_in(extra_path, "plots")

            # Set base path
            base_path = fs.join(extra_plot_path, name)

            # Determine path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(extra_plot_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                plot_path = base_path + "." + extension
            else: plot_path = base_path

        # Determine path
        else:

            # Get paths
            extra_path = fs.create_directory_in(self.maps_sub_path, self.extra_maps_name)
            extra_plot_path = fs.create_directory_in(extra_path, "plots")

            # Set base path
            base_path = fs.join(extra_plot_path, name)

            # Determine plot file path
            if add_extension:
                if clear_other_formats: fs.remove_files_in_path(extra_plot_path, exact_name=name, extension=sequences.all_except_one(["pdf", "png"], extension))
                plot_path = base_path + "." + extension
            else: plot_path = base_path

        # Return the plot file path
        return plot_path

    # -----------------------------------------------------------------

    def plot_extra_maps(self, cmap="viridis", scale="log", format="pdf", cropping_factor=1.3, scales=None,
                        share_limits=True, mask_negatives=False, clear_other_formats=True, show_axes=False,
                        transparent=True, interval=None, strict_vmin=None, strict_vmax=None, soft_vmin=None,
                        soft_vmax=None, methods=None, not_methods=None, soft_max_scaling=1.5):

        """
        This function ...
        :param cmap:
        :param scale:
        :param format:
        :param cropping_factor:
        :param scales:
        :param share_limits:
        :param mask_negatives:
        :param clear_other_formats:
        :param show_axes:
        :param transparent:
        :param interval:
        :param strict_vmin:
        :param strict_vmax:
        :param soft_vmin:
        :param soft_vmax:
        :param methods:
        :param not_methods:
        :param soft_max_scaling:
        :return:
        """

        # Inform the user
        log.info("Plotting the extra maps ...")

        # Determine min, max etc.
        vmin, vmax, passed_min, passed_max, soft_min, soft_max = set_vminvmax(interval=interval, strict_vmin=strict_vmin, strict_vmax=strict_vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax, share_limits=share_limits)

        # Loop over the methods
        for method in self.extra_maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.extra_maps[method]):

                # Plot method?
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting extra maps for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.extra_maps[method]:

                    # Debugging
                    log.debug("Plotting the '" + name + "' extra map ...")

                    # Determine path
                    plot_path = self.get_path_for_extra_map_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If the plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    # Determine scale for this image
                    if scales is not None and method in scales:
                        if isinstance(scales[method], dict):
                            if name in scales[method]: frame_scale = scales[method][name]
                            else: frame_scale = scale
                        else: frame_scale = scales[method]
                    else: frame_scale = scale

                    # Plot
                    minmax = get_minmax(interval=interval, vmin=vmin, vmax=vmax, passed_min=passed_min, passed_max=passed_max, share_limits=share_limits, default_interval="pts")
                    frame = self.extra_maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary
                    vmini, vmaxi = plotting.plot_frame(frame, crop_to=self.truncation_box,
                                                     cropping_factor=cropping_factor,
                                                     truncate_outside=self.truncation_ellipse, path=plot_path,
                                                     format=format, interval=minmax,
                                                     scale=frame_scale, cmap=cmap, normalize_in=self.truncation_ellipse,
                                                     colorbar=True, mask_negatives=mask_negatives, show_axes=show_axes,
                                                     transparent=transparent, soft_min=soft_min, soft_max=soft_max,
                                                     soft_max_scaling=soft_max_scaling)
                    if share_limits: vmin, vmax = vmini, vmaxi

                # End of method: reset vmin and vmax
                vmin, vmax, soft_min, soft_max = reset_vminvmax(soft_min, soft_max, interval=interval, strict_vmin=strict_vmin, strict_vmax=strict_vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax, share_limits=share_limits)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting the '" + method + "' extra map ...")

                # Determine path
                plot_path = self.get_path_for_extra_map_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Determine scale for this image
                if scales is not None and method in scales: frame_scale = scales[method]
                else: frame_scale = scale

                # Plot
                minmax = get_minmax(interval=interval, vmin=vmin, vmax=vmax, passed_min=passed_min, passed_max=passed_max, share_limits=share_limits, default_interval="pts")
                frame = self.extra_maps[method]
                if isinstance(frame, Image): frame = frame.primary
                vmini, vmaxi = plotting.plot_frame(frame, crop_to=self.truncation_box,
                                                 cropping_factor=cropping_factor,
                                                 truncate_outside=self.truncation_ellipse, path=plot_path,
                                                 format=format,
                                                 interval=minmax, scale=frame_scale, cmap=cmap,
                                                 normalize_in=self.truncation_ellipse, colorbar=True,
                                                 mask_negatives=mask_negatives, show_axes=show_axes,
                                                 transparent=transparent, soft_min=soft_min, soft_max=soft_max,
                                                 soft_max_scaling=soft_max_scaling)
                if share_limits: vmin, vmax = vmini, vmaxi

    # -----------------------------------------------------------------

    def plot_contours(self, filled=False, nlevels=10, format="pdf", cropping_factor=1.3, clear_other_formats=True,
                      show_axes=False, transparent=True, methods=None, not_methods=None):

        """
        This function ...
        :param filled:
        :param nlevels:
        :param format:
        :param cropping_factor:
        :param clear_other_formats:
        :param show_axes:
        :param transparent:
        :param methods:
        :param not_methods:
        :return:
        """

        # Inform the user
        log.info("Plotting the contours of the maps ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Plot method?
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting contours for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting contours for the '" + name + "' map ...")

                    # Determine the path
                    plot_path = self.get_path_for_contour_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If plot already exists and we don't have to replot
                    #if fs.is_file(plot_path) and not self.config.replot: continue

                    frame = self.maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary

                    # Plot
                    # if filled: plotting.plot_filled_frame_contours(self.maps[method][name], path=plot_path, nlevels=nlevels,
                    #                                                crop_to=self.truncation_box, cropping_factor=cropping_factor,
                    #                                                truncate_outside=self.truncation_ellipse)
                    # else: plotting.plot_frame_contours(self.maps[method][name], path=plot_path, nlevels=nlevels,
                    #                                    crop_to=self.truncation_box, cropping_factor=cropping_factor,
                    #                                    truncate_outside=self.truncation_ellipse)

                    if not fs.is_file(plot_path) or self.config.replot:

                        plotting.plot_frame_contours(frame, path=plot_path, nlevels=nlevels,
                                                    crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                    truncate_outside=self.truncation_ellipse, show_axes=show_axes,
                                                    transparent=transparent)

                    if not filled: continue
                    filled_plot_path = self.get_path_for_contour_plot(name, method, extension=format, suffix="_filled", clear_other_formats=clear_other_formats)
                    if not fs.is_file(filled_plot_path) or self.config.replot:

                        plotting.plot_filled_frame_contours(frame, path=filled_plot_path, nlevels=nlevels,
                                                            crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                            truncate_outside=self.truncation_ellipse, show_axes=show_axes,
                                                            transparent=transparent)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting contours for the '" + method + "' map ...")

                # Determine the path
                plot_path = self.get_path_for_contour_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If plot already exists and we don't have to remake
                #if fs.is_file(plot_path) and not self.config.replot: continue

                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary

                if not fs.is_file(plot_path) or self.config.replot:

                    plotting.plot_frame_contours(frame, path=plot_path, nlevels=nlevels,
                                                 crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                 truncate_outside=self.truncation_ellipse, show_axes=show_axes,
                                                 transparent=transparent)

                if not filled: continue
                filled_plot_path = self.get_path_for_contour_plot(method, extension=format, suffix="_filled", clear_other_formats=clear_other_formats)
                if not fs.is_file(filled_plot_path) or self.config.replot:

                    plotting.plot_filled_frame_contours(frame, path=filled_plot_path, nlevels=nlevels,
                                                        crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                        truncate_outside=self.truncation_ellipse, show_axes=show_axes,
                                                        transparent=transparent)

    # -----------------------------------------------------------------

    def plot_profiles(self, nbins=20, format="pdf", clear_other_formats=True, transparent=True, methods=None, not_methods=None):

        """
        Thisn function ...
        :param nbins:
        :param format:
        :param clear_other_formats:
        :param transparent:
        :param methods:
        :param not_methods:
        :return:
        """

        # Inform the user
        log.info("Plotting the radial profiles of the maps ...")

        # Get properties
        angle = self.disk_ellipse.angle # or truncation ellipse
        ratio = self.disk_ellipse.semiminor / self.disk_ellipse.semimajor

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # plot method?
                if methods is not None and method not in methods: continue
                if not_methods is not None and method in not_methods: continue

                # Debugging
                log.debug("Plotting radial profiles for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting radial profiles for the '" + name + "' map ...")

                    # Determine the path
                    plot_path = self.get_path_for_profile_plot(name, method, extension=format, clear_other_formats=clear_other_formats)

                    # If plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    frame = self.maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary

                    # Plot
                    plotting.plot_radial_profile(frame, self.galaxy_center, angle, ratio, nbins=nbins, path=plot_path,
                                                 max_radius=self.truncation_ellipse.semimajor, transparent=transparent)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting radial profiles for the '" + method + "' map ...")

                # Determine the path
                plot_path = self.get_path_for_profile_plot(method, extension=format, clear_other_formats=clear_other_formats)

                # If plot already exist and we don't have to remake
                if fs.is_file(plot_path) and not self.config.replot: continue

                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary

                # Plot
                plotting.plot_radial_profile(frame, self.galaxy_center, angle, ratio, nbins=nbins, path=plot_path,
                                             max_radius=self.truncation_ellipse.semimajor, transparent=transparent)

    # -----------------------------------------------------------------

    def check_origins(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the origins dictionary ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                if method not in self.origins: raise ValueError("'" + method + "' section of the origins is missing")

                # Loop over the maps
                for name in self.maps[method]:
                    if name not in self.origins[method]: raise ValueError("Origin for '" + method + "/" + name + "' map is not defined")

            # No different methods
            else:
                if method not in self.origins: raise ValueError("Origin for '" + method + "' map is not defined")

    # -----------------------------------------------------------------

    def check_methods(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Checking the methods dictionary ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                if method not in self.methods: raise ValueError("'" + method + "' section of the methods is missing")

                # Loop over the maps
                for name in self.maps[method]:
                    if name not in self.methods[method]: raise ValueError("Method for '" + method + "/" + name + "' map is not defined")

            # No different methods
            else:
                if method not in self.methods: raise ValueError("Method for '" + method + "' map is not defined")

    # -----------------------------------------------------------------

    @abstractmethod
    def load_collection(self):
        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def collection(self):
        return self.load_collection()

    # -----------------------------------------------------------------

    @abstractmethod
    def load_static_collection(self):
        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def static_collection(self):
        return self.load_static_collection()

    # -----------------------------------------------------------------

    # ORIGINS

    def get_colours_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_colours_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_ssfr_origins(self, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        return self.collection.get_ssfr_origins(flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_tir_origins(self, flatten=False, methods=None):

        """
        This function ...
        :param flatten:
        :param methods:
        :return:
        """

        return self.collection.get_tir_origins(flatten=flatten, methods=methods)

    # -----------------------------------------------------------------

    def get_attenuation_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_attenuation_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_old_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_old_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_young_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_young_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_ionizing_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_ionizing_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_dust_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_dust_origins(flatten=flatten)

    # -----------------------------------------------------------------

    # FILTERS

    def get_all_colours_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_colours_filters()

    # -----------------------------------------------------------------

    def get_all_ssfr_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_ssfr_filters()

    # -----------------------------------------------------------------

    def get_all_tir_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_tir_filters()

    # -----------------------------------------------------------------

    def get_all_attenuation_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_attenuation_filters()

    # -----------------------------------------------------------------

    def get_all_old_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_old_filters()

    # -----------------------------------------------------------------

    def get_all_young_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_young_filters()

    # -----------------------------------------------------------------

    def get_all_ionizing_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_ionizing_filters()

    # -----------------------------------------------------------------

    def get_all_dust_filters(self):

        """
        This function ..
        :return:
        """

        return self.collection.get_all_dust_filters()

    # -----------------------------------------------------------------

    # METHODS

    def get_colours_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_colours_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_ssfr_methods(self, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        return self.collection.get_ssfr_methods(flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_tir_methods(self, flatten=False, methods=None):

        """
        This function ...
        :param flatten:
        :param methods:
        :return:
        """

        return self.collection.get_tir_methods(flatten=flatten, methods=methods)

    # -----------------------------------------------------------------

    def get_attenuation_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_attenuation_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_old_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_old_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_young_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_young_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_ionizing_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_ionizing_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_dust_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_dust_methods(flatten=flatten)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_fuv_attenuation_maps(flatten=flatten)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_fuv_attenuation_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_and_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_fuv_attenuation_maps_and_origins(flatten=flatten)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_origins_and_methods(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        return self.collection.get_fuv_attenuation_maps_origins_and_methods(flatten=flatten, cortese=cortese, buat=buat)

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_origins_methods_and_nans(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        return self.collection.get_fuv_attenuation_maps_origins_methods_and_nans(flatten=flatten, cortese=cortese, buat=buat)

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

    @property
    def maps_colours_path(self):
        return self.collection.maps_colours_path

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):
        return self.collection.maps_colours_name

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_path(self):
        return self.collection.maps_ssfr_path

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):
        return self.collection.maps_ssfr_name

    # -----------------------------------------------------------------

    @property
    def maps_tir_path(self):
        return self.collection.maps_tir_path

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):
        return self.collection.maps_tir_name

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_path(self):
        return self.collection.maps_attenuation_path

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):
        return self.collection.maps_attenuation_name

    # -----------------------------------------------------------------

    @property
    def maps_old_path(self):
        return self.collection.maps_old_path

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):
        return self.collection.maps_old_name

    # -----------------------------------------------------------------

    @property
    def maps_young_path(self):
        return self.collection.maps_young_path

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):
        return self.collection.maps_young_name

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_path(self):
        return self.collection.maps_ionizing_path

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):
        return self.collection.maps_ionizing_name

    # -----------------------------------------------------------------

    @property
    def maps_dust_path(self):
        return self.collection.maps_dust_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):
        return self.collection.maps_dust_name

    # -----------------------------------------------------------------

    @property
    def maps_sub_paths(self):
        return self.collection.maps_sub_paths

    # -----------------------------------------------------------------

    @property
    def maps_sub_names(self):
        return self.collection.maps_sub_names

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

    def get_frame(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse filter
        fltr = parse_filter(fltr)

        # Get frame
        return self.get_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_hot_dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_hot_dust_maps()

    # -----------------------------------------------------------------

    def get_hot_dust_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_hot_dust_origins()

    # -----------------------------------------------------------------

    def get_hot_dust_methods(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_hot_dust_methods()

# -----------------------------------------------------------------

class MapsComponent(GalaxyModelingComponent):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def load_collection(self):
        return MapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    def load_static_collection(self):
        return StaticMapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def selection(self):
        return ComponentMapsSelection(self.maps_components_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_selection(self):
        return StaticComponentMapsSelection(self.maps_components_path)

    # -----------------------------------------------------------------
    # MAPS
    # -----------------------------------------------------------------

    def get_component_old_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_map_paths()

    # -----------------------------------------------------------------

    def get_component_old_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_maps()

    # -----------------------------------------------------------------

    def get_component_young_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_map_paths()

    # -----------------------------------------------------------------

    def get_component_young_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_maps()

    # -----------------------------------------------------------------

    def get_component_ionizing_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_ionizing_map_paths()

    # -----------------------------------------------------------------

    def get_component_ionizing_maps(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.static_selection.get_ionizing_maps()

    # -----------------------------------------------------------------

    def get_component_dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_map_paths()

    # -----------------------------------------------------------------

    def get_component_dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_maps()

    # -----------------------------------------------------------------

    # DEPROJECTED

    # -----------------------------------------------------------------

    def get_component_old_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_deprojected_paths()

    # -----------------------------------------------------------------

    def get_component_old_deprojected(self):

        """
        Thisnf unction ...
        :return:
        """

        return self.static_selection.get_old_deprojected()

    # -----------------------------------------------------------------

    def get_component_young_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_deprojected_paths()

    # -----------------------------------------------------------------

    def get_component_young_deprojected(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_deprojected()

    # -----------------------------------------------------------------

    def get_component_ionizing_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_ionizing_deprojected_paths()

    # -----------------------------------------------------------------

    def get_component_ionizing_deprojected(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_selection.get_ionizing_deprojected()

    # -----------------------------------------------------------------

    def get_component_dust_deprojected_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_deprojected_paths()

    # -----------------------------------------------------------------

    def get_component_dust_deprojected(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_selection.get_dust_deprojected()

    # -----------------------------------------------------------------

    # DEPROJECTED WITH SKIRT

    # -----------------------------------------------------------------

    def get_component_old_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_deprojected_skirt_paths()

    # -----------------------------------------------------------------

    def get_component_old_deprojected_skirt(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_selection.get_old_deprojected_skirt()

    # -----------------------------------------------------------------

    def get_component_young_deprojected_skirt_paths(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_selection.get_young_deprojected_skirt_paths()

    # -----------------------------------------------------------------

    def get_component_young_deprojected_skirt(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_selection.get_young_deprojected_skirt()

    # -----------------------------------------------------------------

    def get_component_ionizing_deprojected_skirt_paths(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_selection.get_ionizing_deprojected_skirt_paths()

    # -----------------------------------------------------------------

    def get_component_ionizing_deprojected_skirt(self):

        """
        Thisn function ...
        :return:
        """

        return self.static_selection.get_ionizing_deprojected_skirt()

    # -----------------------------------------------------------------

    def get_component_dust_deprojected_skirt_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_deprojected_skirt_paths()

    # -----------------------------------------------------------------

    def get_component_dust_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_deprojected_skirt()

    # -----------------------------------------------------------------

    # EDGE-ON

    # -----------------------------------------------------------------

    def get_component_old_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_edgeon_paths()

    # -----------------------------------------------------------------

    def get_component_old_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_old_edgeon()

    # -----------------------------------------------------------------

    def get_component_young_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_edgeon_paths()

    # -----------------------------------------------------------------

    def get_component_young_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_young_edgeon()

    # -----------------------------------------------------------------

    def get_component_ionizing_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_ionizing_edgeon_paths()

    # -----------------------------------------------------------------

    def get_component_ionizing_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.static_selection.get_ionizing_edgeon()

    # -----------------------------------------------------------------

    def get_component_dust_edgeon_paths(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.get_dust_edgeon_paths()

    # -----------------------------------------------------------------

    def get_component_dust_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.static_selection.get_dust_edgeon()

    # -----------------------------------------------------------------

    @property
    def old_component_map_names(self):

        """
        Thisfunction ...
        :return:
        """

        return self.static_selection.old_map_names

    # -----------------------------------------------------------------

    @property
    def young_component_map_names(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.young_map_names

    # -----------------------------------------------------------------

    @property
    def ionizing_component_map_names(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.ionizing_map_names

    # -----------------------------------------------------------------

    @property
    def dust_component_map_names(self):

        """
        This function ...
        :return:
        """

        return self.static_selection.dust_map_names

    # -----------------------------------------------------------------

    def get_component_map_paths(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        paths = dict()

        # Flattened?
        if flatten:

            # Add old
            old_paths = self.get_component_old_map_paths()
            for name in old_paths:
                new_name = "old_" + name
                paths[new_name] = old_paths[name]

            # Add young
            young_paths = self.get_component_young_map_paths()
            for name in young_paths:
                new_name = "young_" + name
                paths[new_name] = young_paths[name]

            # Add ionizing
            ionizing_paths = self.get_component_ionizing_map_paths()
            for name in ionizing_paths:
                new_name = "ionizing_" + name
                paths[new_name] = ionizing_paths[name]

            # Add dust
            dust_paths = self.get_component_dust_map_paths()
            for name in dust_paths:
                new_name = "dust_" + name
                paths[new_name] = dust_paths[name]

        # Not flattened
        else:

            paths["old"] = self.get_component_old_map_paths()
            paths["young"] = self.get_component_young_map_paths()
            paths["ionizing"] = self.get_component_ionizing_map_paths()
            paths["dust"] = self.get_component_dust_map_paths()

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    @property
    def old_component_maps_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_old_component_path

    # -----------------------------------------------------------------

    @property
    def young_component_maps_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_young_component_path

    # -----------------------------------------------------------------

    @property
    def ionizing_component_maps_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_ionizing_component_path

    # -----------------------------------------------------------------

    @property
    def dust_component_maps_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_dust_component_path

    # -----------------------------------------------------------------

    @property
    def maps_html_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_html_path

    # -----------------------------------------------------------------

    @property
    def all_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.all_maps_html_page_path

    # -----------------------------------------------------------------

    @property
    def maps_summary_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_summary_html_page_path

    # -----------------------------------------------------------------

    @property
    def clip_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.clip_maps_html_page_path

    # -----------------------------------------------------------------

    @property
    def old_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.old_maps_html_page_path

    # -----------------------------------------------------------------

    @property
    def young_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.young_maps_html_page_path

    # -----------------------------------------------------------------

    @property
    def ionizing_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.ionizing_maps_html_page_path

    # -----------------------------------------------------------------

    @property
    def dust_maps_html_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.dust_maps_html_page_path

# -----------------------------------------------------------------

class MapMakingComponent(MapMakerBase, MapsComponent):
    
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
        #super(MapMakingComponent, self).__init__(*args, **kwargs)
        MapsComponent.__init__(self, no_config=True)
        MapMakerBase.__init__(self, *args, **kwargs)

        # -- Attributes --

        # The error maps
        self.error_maps = dict()

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
        #super(MapMakingComponent, self).setup(**kwargs)
        MapMakerBase.setup(self, **kwargs)
        MapsComponent.setup(self, **kwargs)

        # Perform some checks to see whether previous output is present
        if self.config.checks: self.perform_checks()

    # -----------------------------------------------------------------

    def type_for_map_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine map type (colour/sSFR/TIR/attenuation/old/young/ionizing/dust)
        relative = fs.relative_to(path, self.maps_path)
        map_type = fs.base_directory(relative)
        return map_type

    # -----------------------------------------------------------------

    def method_for_map_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get type
        map_type = self.type_for_map_path(path)

        map_directory = fs.directory_of(path)
        map_directory_name = fs.name(map_directory)

        if map_directory_name == map_type: return None
        else: return map_directory_name

    # -----------------------------------------------------------------

    def method_or_type_for_map_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        method = self.method_for_map_path(path)

        if method is None: return self.type_for_map_path(path)
        else: return method

    # -----------------------------------------------------------------

    def perform_checks(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking existing maps ...")

        # Loop over the finished maps commands
        for command_name in self.history.finished_maps_commmands:

            # Get sub path for command
            path = self.sub_path_for_command(command_name)
            name = fs.name(path)

            # Check whether the path is present and not empty
            if not fs.is_directory(path):

                # Warning
                log.error("The '" + name + "' directory does not exist from the '" + command_name + "' output")

                # Remove
                self.remove_output_and_history_after_and_including(command_name)

            # If the directory is empty
            elif fs.is_empty(path):

                # Warning
                log.error("The '" + name + "' directory from the '" + command_name + "' output is empty")

                # Remove
                self.remove_output_and_history_after_and_including(command_name)

            # Directory present and not empty, loop over the methods
            else:

                methods = self.methods_for_sub_name(name)
                if methods is None: continue
                for method in methods:

                    # Determine the path
                    method_path = fs.join(path, method)

                    #print(method_path)

                    # NO: WE DON'T DO ALL METHODS (E.G. DUST->BLACK-BODY)
                    # Check if the path if present
                    #if not fs.is_directory(method_path):
                        # Warning
                        #log.error("The '" + name + "/" + method + "' directory does not exist from the '" + command_name + "' output")
                        # Remove
                        #self.remove_output_and_history_after_and_including(command_name)

                    # If the directory is empty
                    #elif fs.is_empty(method_path):
                    if fs.is_directory(method_path) and fs.is_empty(method_path):

                        # Warning
                        log.error("The '" + name + "/" + method + "' directory from the '" + command_name + "' output is empty")

                        # Remove
                        self.remove_output_and_history_after_and_including(command_name)

    # -----------------------------------------------------------------

    def remove_output_and_history_after_and_including(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        log.warning("Removing the '" + command_name + "' and successive commands from the history ...")
        commands = maps_commands_after_and_including(command_name)
        existing_commands = [command for command in commands if command in self.history.commands]
        yn = prompt_proceed("Remove the history and output of the following commands: " + tostr(existing_commands) + "?")
        if yn:

            # Remove from history
            for command in existing_commands:

                # Debugging
                log.debug("Removing the '" + command + "' command from the modeling history ...")

                # Remove
                self.history.remove_entries(command)

            self.history.save()

            # Remove output
            for command in existing_commands:

                path = self.sub_path_for_command(command)

                # Debugging
                log.debug("Removing the output path [" + path + "] for the '" + command + "' command ...")

                # Remove
                fs.remove_directory(path)

            # Exit
            exit()

        # User answered no
        else: raise RuntimeError("Cannot proceed")

    # -----------------------------------------------------------------

    def sub_path_for_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        name = self.sub_name_for_command(command)
        return fs.join(self.maps_raw_path, name)

    # -----------------------------------------------------------------

    def sub_name_for_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        if command == "make_colours_maps": return self.maps_colours_name
        elif command == "make_ssfr_maps": return self.maps_ssfr_name
        elif command == "make_tir_maps": return self.maps_tir_name
        elif command == "make_attenuation_maps": return self.maps_attenuation_name
        elif command == "make_old_stellar_maps": return self.maps_old_name
        elif command == "make_young_stellar_maps": return self.maps_young_name
        elif command == "make_ionizing_stellar_maps": return self.maps_ionizing_name
        elif command == "make_dust_map": return self.maps_dust_name
        else: raise ValueError("Invalid commands: " + command)

    # -----------------------------------------------------------------

    def command_for_sub_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name == self.maps_colours_name: return "make_colours_maps"
        elif name == self.maps_ssfr_name: return "make_ssfr_maps"
        elif name == self.maps_tir_name: return "make_tir_maps"
        elif name == self.maps_attenuation_name: return "make_attenuation_maps"
        elif name == self.maps_old_name: return "make_old_stellar_maps"
        elif name == self.maps_young_name: return "make_young_stellar_maps"
        elif name == self.maps_ionizing_name: return "make_ionizing_stellar_maps"
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

    def class_for_sub_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the command
        command = self.command_for_sub_name(name)

        # Resolve the PTS command
        subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(command)

        # Return the class
        cls = introspection.get_class(class_module_path, class_name)

        # Return the class
        return cls

    # -----------------------------------------------------------------

    def methods_for_sub_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the command
        command = self.command_for_sub_name(name)

        # Resolve the PTS command
        subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(command)

        # Try to get the methods
        methods = introspection.try_importing_variable(class_module_path, "methods")

        # Return the methods
        return methods

    # -----------------------------------------------------------------

    def has_methods_for_sub_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.methods_for_sub_name(name) is not None

    # -----------------------------------------------------------------

    def get_origins_sub_name(self, name, flatten=False, method=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param method:
        :return:
        """

        return self.collection.get_origins_sub_name(name, flatten=flatten, method=method)

    # -----------------------------------------------------------------

    def get_origins_sub_names(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_origins_sub_names(flatten=flatten)

    # -----------------------------------------------------------------

    def get_all_filters_sub_name(self, name, method=None):

        """
        This function ...
        """

        return self.collection.get_all_filters_sub_name(name, method=method)

    # -----------------------------------------------------------------

    def get_all_filters(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_all_filters()

    # -----------------------------------------------------------------

    def get_maps_sub_names(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_maps_sub_names(flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def has_colour_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_colour_maps

    # -----------------------------------------------------------------

    @property
    def has_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_ssfr_maps

    # -----------------------------------------------------------------

    @property
    def has_tir_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_tir_maps

    # -----------------------------------------------------------------

    @property
    def has_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_attenuation_maps

    # -----------------------------------------------------------------

    @property
    def has_old_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_old_maps

    # -----------------------------------------------------------------

    @property
    def has_young_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_young_maps

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_ionizing_maps

    # -----------------------------------------------------------------

    @property
    def has_dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.has_dust_maps

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

    def get_tir_single_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_tir_single_maps()

    # -----------------------------------------------------------------

    def get_tir_multi_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_tir_multi_maps()

    # -----------------------------------------------------------------

    def get_tir_single_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_tir_single_origins()

    # -----------------------------------------------------------------

    def get_tir_multi_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_tir_multi_origins()

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_buat_fuv_attenuation_maps()

    # -----------------------------------------------------------------

    def get_cortese_fuv_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_cortese_fuv_attenuation_maps()

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_buat_fuv_attenuation_origins()

    # -----------------------------------------------------------------

    def get_buat_nuv_attenuation_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_buat_nuv_attenuation_origins()

    # -----------------------------------------------------------------

    def get_cortese_fuv_attenuation_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_cortese_fuv_attenuation_origins()

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_map(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.collection.get_old_stellar_bulge_map(fltr)

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_maps(self, framelist=False):

        """
        This function ...
        :param framelist:
        :return:
        """

        return self.collection.get_old_stellar_bulge_maps(framelist=framelist)

    # -----------------------------------------------------------------

    def get_current_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_sub_name, flatten=flatten)

    # -----------------------------------------------------------------

    def get_map_paths_sub_name(self, name, flatten=False, method=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param method:
        :return:
        """

        return self.collection.get_map_paths_sub_name(name, flatten=flatten, method=method)

    # -----------------------------------------------------------------

    def get_current_map_paths(self, flatten=False, method=None):

        """
        This function ...
        :param flatten:
        :param method:
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_sub_name, flatten=flatten, method=method)

    # -----------------------------------------------------------------

    def clear_current_all(self):

        """
        This function ...
        :return:
        """

        # Maps
        self.clear_current_maps()

        # Origins
        self.clear_current_origins()

        # Methods
        self.clear_current_methods()

        # Plots
        self.clear_current_plots()

        # Contour plots
        self.clear_current_contours()

        # Profile plots
        self.clear_current_profiles()

        # Negatives plots
        self.clear_current_negatives()

        # NaNs plots
        self.clear_current_nans()

    # -----------------------------------------------------------------

    def clear_current_all_method(self, method):

        """
        This function ....
        :param method:
        :return:
        """

        # Maps
        self.clear_current_maps_method(method)

        # Origins
        self.clear_current_origins_method(method)

        # Methods
        self.clear_current_methods_method(method)

        # Plots
        self.clear_current_plots_method(method)

        # Contour plots
        self.clear_current_contours_method(method)

        # Profile plots
        self.clear_current_profiles_method(method)

        # Negatives plots
        self.clear_current_negatives_method(method)

        # NaNs plots
        self.clear_current_nans_method(method)

    # -----------------------------------------------------------------

    def clear_current_maps(self):

        """
        This function ...
        :return:
        """

        paths = self.get_map_paths_sub_name(self.maps_sub_name)
        self.clear_maps_from_paths(paths)

    # -----------------------------------------------------------------

    def clear_current_maps_method(self, method):

        """
        Thisf ucntion ...
        :param method:
        :return:
        """

        paths = self.get_map_paths_sub_name(self.maps_sub_name, flatten=True, method=method)
        self.clear_maps_from_paths(paths)

    # -----------------------------------------------------------------

    def clear_maps_from_paths(self, paths):

        """
        Thisnf unction ...
        :param paths:
        :return:
        """

        # Loop over the files
        for name in paths:

            path = paths[name]
            log.warning("Removing the '" + name + "' map ...")
            fs.remove_file(path)

    # -----------------------------------------------------------------

    def clear_current_origins(self):

        """
        This function ...
        :return:
        """

        # Remove the origins file
        origins_path = fs.join(self.maps_sub_path, "origins.txt")
        if fs.is_file(origins_path):
            log.warning("Removing the origins file ...")
            fs.remove_file(origins_path)

    # -----------------------------------------------------------------

    def clear_current_origins_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Remove the origins file
        origins_path = fs.join(self.maps_sub_path, method, "origins.txt")
        if fs.is_file(origins_path):
            log.warning("Removing the origins file ...")
            fs.remove_file(origins_path)

    # -----------------------------------------------------------------

    def clear_current_methods(self):

        """
        This function ...
        :return:
        """

        # Remove the methods file
        methods_path = fs.join(self.maps_sub_path, "methods.txt")
        if fs.is_file(methods_path):
            log.warning("Removing the methods file ...")
            fs.remove_file(methods_path)

    # -----------------------------------------------------------------

    def clear_current_methods_method(self, method):

        """
        Thisf unction ...
        :param method:
        :return:
        """

        # Remove the methods file
        methods_path = fs.join(self.maps_sub_path, method, "methods.txt")
        if fs.is_file(methods_path):
            log.warning("Removing the methods file ...")
            fs.remove_file(methods_path)

    # -----------------------------------------------------------------

    def clear_current_plots(self):

        """
        This function ...
        :return:
        """

        # Determine plots path
        plots_path = fs.join(self.maps_sub_path, self.map_plots_name)
        if not fs.is_directory(plots_path): return

        # Clear
        log.warning("Clearing the plots directory ...")
        fs.clear_directory(plots_path)

    # -----------------------------------------------------------------

    def clear_current_plots_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Determine plots path
        plots_path = fs.join(self.maps_sub_path, method, self.map_plots_name)
        if not fs.is_directory(plots_path): return

        # Clear
        log.warning("Clearing the '" + method + "' plots directory ...")
        fs.clear_directory(plots_path)

    # -----------------------------------------------------------------

    def clear_current_contours(self):

        """
        This function ...
        :return:
        """

        # Determine contours path
        contours_path = fs.join(self.maps_sub_path, self.contour_plots_name)
        if not fs.is_directory(contours_path): return

        # Clear
        log.warning("Clearing the contours directory ...")
        fs.clear_directory(contours_path)

    # -----------------------------------------------------------------

    def clear_current_contours_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Determine contours path
        contours_path = fs.join(self.maps_sub_path, method, self.contour_plots_name)
        if not fs.is_directory(contours_path): return

        # Clear
        log.warning("Clearing the '" + method + "' contours directory ...")
        fs.clear_directory(contours_path)

    # -----------------------------------------------------------------

    def clear_current_profiles(self):

        """
        This function ...
        :return:
        """

        # Determine profiles path
        profiles_path = fs.join(self.maps_sub_path, self.profile_plots_name)
        if not fs.is_directory(profiles_path): return

        # Clear
        log.warning("Clearing the profiles directory ...")
        fs.clear_directory(profiles_path)

    # -----------------------------------------------------------------

    def clear_current_profiles_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Determine profiles path
        profiles_path = fs.join(self.maps_sub_path, method, self.profile_plots_name)
        if not fs.is_directory(profiles_path): return

        # Clear
        log.warning("Clearing the '" + method + "' profiles directory ...")
        fs.clear_directory(profiles_path)

    # -----------------------------------------------------------------

    def clear_current_negatives(self):

        """
        This function ...
        :return:
        """

        # Determine negatives path
        negatives_path = fs.join(self.maps_sub_path, self.negatives_plots_name)
        if not fs.is_directory(negatives_path): return

        # Clear
        log.warning("Clearing the negatives directory ...")
        fs.clear_directory(negatives_path)

    # -----------------------------------------------------------------

    def clear_current_negatives_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Determine negatives path
        negatives_path = fs.join(self.maps_sub_path, method, self.negatives_plots_name)
        if not fs.is_directory(negatives_path): return

        # Clear
        log.warning("Clearing the '" + method + "' negatives directory ...")
        fs.clear_directory(negatives_path)

    # -----------------------------------------------------------------

    def clear_current_nans(self):

        """
        This function ...
        :return:
        """

        # Determine nans path
        nans_path = fs.join(self.maps_sub_path, self.nans_plots_name)
        if not fs.is_directory(nans_path): return

        # Clear
        log.warning("Clearing the NaNs directory ...")
        fs.clear_directory(nans_path)

    # -----------------------------------------------------------------

    def clear_current_nans_method(self, method):

        """
        Thisf unction ...
        :param method:
        :return:
        """

        # Determine nans path
        nans_path = fs.join(self.maps_sub_path, method, self.nans_plots_name)
        if not fs.is_directory(nans_path): return

        # Clear
        log.warning("Clearing the '" + method + "' NaNs directory ...")
        fs.clear_directory(nans_path)

    # -----------------------------------------------------------------

    @property
    def colours_scale(self):

        """
        This function ...
        :return:
        """

        return "squared"

    # -----------------------------------------------------------------

    @property
    def colours_cmap(self):

        """
        This function ...
        :return:
        """

        return "viridis"

    # -----------------------------------------------------------------

    @property
    def colours_js9_cmap(self):

        """
        This function ...
        :return:
        """

        return "viridis"

    # -----------------------------------------------------------------

    # @property
    # def scales(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     scales = dict()
    #
    #     from ...magic.tools.colours import is_fir_colour, is_fir_or_submm_colour
    #
    #     # Loop over the colours
    #     for colour_name in self.maps:
    #
    #         # Check whether it is a FIR colour
    #         if is_fir_colour(colour_name) or is_fir_or_submm_colour(colour_name):
    #             #around_zero = True
    #             scale = "linear"
    #         else:
    #             #around_zero = False
    #             scale = "log"
    #
    #         # Set the scale
    #         scales[colour_name] = scale
    #
    #     # Return the scales
    #     return scales

    # -----------------------------------------------------------------

    @property
    def colours_share_limits(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    def plot_colours(self, maps=True, contours=True, profiles=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the colour maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.colours_scale, cmap=self.colours_cmap, share_limits=self.colours_share_limits, format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_scale(self):

        """
        This function ...
        :return:
        """

        return "squared"

    # -----------------------------------------------------------------

    @property
    def ssfr_cmap(self):

        """
        Thisn function ...
        :return:
        """

        return "jet"

    # -----------------------------------------------------------------

    @property
    def ssfr_js9_cmap(self):

        """
        This function ...
        :return:
        """

        return "rainbow"

    # -----------------------------------------------------------------

    @property
    def ssfr_interval(self):

        """
        This function ...
        :return:
        """

        return [0, 10]

    # -----------------------------------------------------------------

    @property
    def ssfr_min(self):

        """
        This function ...
        :return:
        """

        return 0

    # -----------------------------------------------------------------

    @property
    def ssfr_soft_max(self):

        """
        This function ...
        :return:
        """

        return 10

    # -----------------------------------------------------------------

    def plot_ssfr(self, maps=True, contours=True, profiles=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.ssfr_scale, cmap=self.ssfr_cmap, strict_vmin=self.ssfr_min,
                                              soft_vmax=self.ssfr_soft_max, share_limits=False,
                                              format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def tir_scale(self):
        return "log"

    # -----------------------------------------------------------------

    @property
    def tir_cmap(self):
        return "viridis"

    # -----------------------------------------------------------------

    @property
    def tir_js9_cmap(self):
        return "viridis"

    # -----------------------------------------------------------------

    def plot_tir(self, maps=True, contours=True, profiles=True, nans=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the TIR maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.tir_scale, cmap=self.tir_cmap, format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the NaNs masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_scale(self):
        return "linear"

    # -----------------------------------------------------------------

    @property
    def attenuation_cmap(self):
        return "cool"

    # -----------------------------------------------------------------

    @property
    def attenuation_js9_cmap(self):
        return "cool"

    # -----------------------------------------------------------------

    @property
    def tir_to_fuv_cmap(self):
        return "ds9heat"

    # -----------------------------------------------------------------

    def plot_attenuation(self, maps=True, contours=True, profiles=True, extra=True, nans=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param extra:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the attenuation maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.attenuation_scale, cmap=self.attenuation_cmap, format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the extra maps
        if extra: self.plot_extra_maps(scale="sqrt", format=format, clear_other_formats=True, cmap=self.tir_to_fuv_cmap)

        # Plot the NaNs masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def old_scale(self):
        return "log"

    # -----------------------------------------------------------------

    @property
    def old_cmap(self):
        # return "Wistia"
        return "afmhot"

    # -----------------------------------------------------------------

    @property
    def old_js9_cmap(self):
        return "heat"

    # -----------------------------------------------------------------

    def plot_old(self, maps=True, contours=True, profiles=True, negatives=True, nans=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param negatives:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the old stellar maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.old_scale, cmap=self.old_cmap, format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the negatives
        if negatives: self.plot_negatives(format=format, clear_other_formats=True, methods=["disk"], count_within=self.central_ellipse)

        # Plot the NaNs masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def dust_scales(self):

        """
        Thisf unction ...
        :return:
        """

        from .dust import blackbody, attenuation, hot
        scales = dict()
        scales[blackbody] = "linear"
        scales[attenuation] = "linear"
        scales[hot] = "log"
        return scales

    # -----------------------------------------------------------------

    @property
    def dust_cmaps(self):

        """
        Thisf unction ...
        :return:
        """

        from .dust import blackbody, attenuation, hot
        cmaps = dict()
        cmaps[blackbody] = "jet"
        cmaps[attenuation] = "gist_ncar"
        cmaps[hot] = "jet"
        return cmaps

    # -----------------------------------------------------------------

    @property
    def dust_js9_cmap(self):

        """
        Thisf unction ...
        :return:
        """

        return "a"

    # -----------------------------------------------------------------

    def plot_dust(self, maps=True, contours=True, profiles=True, negatives=True, nans=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param negatives:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the dust maps ...")

        # Plot the maps
        if maps: self.plot_maps(scales=self.dust_scales, cmaps=self.dust_cmaps, format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the negative pixel masks
        if negatives: self.plot_negatives(format=format, clear_other_formats=True, methods=["hot"], count_within=self.central_ellipse)

        # Plot the NaNs masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def young_scale(self):

        """
        This function ...
        :return:
        """

        #return "sqrt"
        return "log"

    # -----------------------------------------------------------------

    @property
    def young_cmap(self):

        """
        This function ...
        :return:
        """

        return "plasma"

    # -----------------------------------------------------------------

    @property
    def young_js9_cmap(self):

        """
        This function ...
        :return:
        """

        return "plasma"

    # -----------------------------------------------------------------

    def plot_young(self, maps=True, contours=True, profiles=True, negatives=True, nans=True, format="pdf"):

        """
        Thins function ...
        :param maps:
        :param contours:
        :param profiles:
        :param negatives:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the young stellar maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.young_scale, cmap=self.young_cmap, mask_negatives=True,
                                              format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the negative pixel masks
        if negatives: self.plot_negatives(format=format, clear_other_formats=True, count_within=self.central_ellipse)

        # Plot the NaN pixel masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

    # -----------------------------------------------------------------

    @property
    def ionizing_scale(self):
        return "log"

    # -----------------------------------------------------------------

    @property
    def ionizing_cmap(self):
        # return "ds9he"
        return "summer"

    # -----------------------------------------------------------------

    @property
    def ionizing_js9_cmap(self):
        return "cool"

    # -----------------------------------------------------------------

    @property
    def halpha_to_hot_cmap(self):
        return "jet"

    # -----------------------------------------------------------------

    @property
    def halpha_to_hot_vmax(self):
        return 1.

    # -----------------------------------------------------------------

    def plot_ionizing(self, maps=True, contours=True, profiles=True, extra=True, negatives=True, nans=True, format="pdf"):

        """
        This function ...
        :param maps:
        :param contours:
        :param profiles:
        :param extra:
        :param negatives:
        :param nans:
        :param format:
        :return:
        """

        # Inform the user
        log.info("Plotting the ionizing stellar maps ...")

        # Plot the maps
        if maps: self.plot_maps(scale=self.ionizing_scale, cmap=self.ionizing_cmap,
                                              format=format, clear_other_formats=True)

        # Plot the contours
        if contours: self.plot_contours(filled=True, format=format, clear_other_formats=True)

        # Plot the radial profiles
        if profiles: self.plot_profiles(format=format, clear_other_formats=True)

        # Plot the extra maps
        if extra: self.plot_extra_maps(scale="log", format=format, clear_other_formats=True, cmap=self.halpha_to_hot_cmap, strict_vmax=self.halpha_to_hot_vmax)

        # Plot the negative masks
        if negatives: self.plot_negatives(format=format, clear_other_formats=True, count_within=self.central_ellipse)

        # Plot the NaNs masks
        if nans: self.plot_nans(format=format, clear_other_formats=True)

# -----------------------------------------------------------------

def select_maps(maps, title, return_names=False):

    """
    This function ...
    :param maps:
    :param title:
    :param return_names:
    :return:
    """

    # Select names interactively
    names = prompt_string_list("names", title, choices=maps.keys())

    # New maps
    new_maps = dict()

    # Get selection
    for name in names: new_maps[name] = maps[name]

    # Return the selected maps
    if return_names: return new_maps, names
    else: return new_maps

# -----------------------------------------------------------------

def set_vminvmax(interval=None, strict_vmin=None, strict_vmax=None, soft_vmin=None, soft_vmax=None, share_limits=True):

    """
    This function ...
    :param interval:
    :param strict_vmin:
    :param strict_vmax:
    :param soft_vmin:
    :param soft_vmax:
    :param share_limits:
    :return:
    """

    # Initialize flags
    passed_min = False
    passed_max = False
    soft_min = False
    soft_max = False

    # The vmin and vmax
    if interval is not None and (isinstance(interval, list) or isinstance(interval, tuple)):

        if not share_limits: raise ValueError("share_limits can only be true when interval passed is vmin,vmax")
        vmin, vmax = interval
        passed_min = True
        passed_max = True

    else:

        # MIN

        if strict_vmin is not None:

            vmin = strict_vmin
            passed_min = True

        elif soft_vmin is not None:

            vmin = soft_vmin
            soft_min = True
            passed_min = True

        else: vmin = None

        # MAX

        if strict_vmax is not None:

            vmax = strict_vmax
            passed_max = True

        elif soft_vmax is not None:

            vmax = soft_vmax
            soft_max = True
            passed_max = True

        else: vmax = None

    # Return
    return vmin, vmax, passed_min, passed_max, soft_min, soft_max

# -----------------------------------------------------------------

def get_minmax(interval=None, vmin=None, vmax=None, passed_min=False, passed_max=False, share_limits=True, default_interval="pts"):

    """
    This function ...
    :param interval:
    :param vmin:
    :param vmax:
    :param passed_min:
    :param passed_max:
    :param share_limits:
    :param default_interval:
    :return:
    """

    # Share limits?
    if share_limits and vmin is not None and vmax is not None: minmax = [vmin, vmax]

    # Passed min and passed max: fixed interval
    elif passed_min and passed_max: minmax = [vmin, vmax]

    # Interval is set
    elif interval is not None: minmax = interval

    # Let PTS determine interval
    else: minmax = default_interval

    # Return
    return minmax

# -----------------------------------------------------------------

def reset_vminvmax(soft_min, soft_max, interval=None, strict_vmin=None, strict_vmax=None, soft_vmin=None, soft_vmax=None, share_limits=True):

    """
    This function ...
    :param soft_min:
    :param soft_max:
    :param interval:
    :param strict_vmin:
    :param strict_vmax:
    :param soft_vmin:
    :param soft_vmax:
    :param share_limits:
    :return:
    """

    if interval is not None and (isinstance(interval, list) or isinstance(interval, tuple)):

        if not share_limits: raise ValueError("share_limits can only be true when interval passed is vmin,vmax")
        vmin, vmax = interval

    else:

        if strict_vmin is not None: vmin = strict_vmin

        elif soft_vmin is not None:

            vmin = soft_vmin
            soft_min = True

        else: vmin = None

        if strict_vmax is not None: vmax = strict_vmax

        elif soft_vmax is not None:

            vmax = soft_vmax
            soft_max = True

        else: vmax = None

    # Return
    return vmin, vmax, soft_min, soft_max

# -----------------------------------------------------------------