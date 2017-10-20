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
from pts.core.tools.utils import lazyproperty
from ..core.environment import colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ...core.basics.configuration import prompt_string_list
from ...core.basics.containers import create_subdict
from ...magic.tools import plotting
from ...magic.core.image import Image

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

        """
        This function ...
        :return:
        """

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

    def get_maps_sub_name(self, name, flatten=False, framelist=False, method=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.collection.get_maps_sub_name(name, flatten=flatten, framelist=framelist, method=method)

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

    def get_current_maps(self, flatten=False, framelist=False, method=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.get_maps_sub_name(self.maps_sub_name, flatten=flatten, framelist=framelist, method=method)

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

    def select_ssfr_maps(self, names=None, methods=None, prompt=False, title="sSFR maps"):

        """
        This function ...
        :param names:
        :param methods:
        :param prompt:
        :param title:
        :return:
        """

        # Get maps
        ssfrs = self.get_ssfr_maps(flatten=True)

        # Origins
        ssfrs_origins = self.get_ssfr_origins(flatten=True)

        # Methods
        ssfrs_methods = self.get_ssfr_methods(flatten=True)

        # Nans
        ssfrs_nans = self.get_ssfr_nans(flatten=True)

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

        # Return
        return ssfrs, ssfrs_origins, ssfrs_methods, ssfrs_nans

    # -----------------------------------------------------------------

    def get_ssfr_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_ssfr_maps(flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_ssfr_nans(self, flatten=False, framelist=False):

        """
        Thisfunction ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_ssfr_nans(flatten=flatten, framelist=framelist)

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

        # Return
        return tirs, tirs_origins, tirs_methods, tirs_nans

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

    def get_attenuation_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist::
        :return:
        """

        return self.collection.get_attenuation_maps(flatten=flatten, framelist=framelist)

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

    def get_dust_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.collection.get_dust_maps(flatten=flatten, framelist=framelist)

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

    def get_old_stellar_disk_origins(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_disk_origins()

    # -----------------------------------------------------------------

    def get_old_stellar_disk_methods(self):

        """
        This function ...
        :return:
        """

        return self.collection.get_old_stellar_disk_methods()

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
                path = fs.join(self.maps_sub_path, method)

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

    def get_path_for_map_plot(self, name, method=None, add_extension=True, extension="pdf"):

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
            extra_path = fs.create_directory_in(path, self.map_plots_name)

            # Determine path
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Determine path
        else:

            extra_path = fs.create_directory_in(self.maps_sub_path, self.map_plots_name)
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_contour_plot(self, name, method=None, add_extension=True, extension="pdf", suffix=""):

        """
        This function ...
        :param name:
        :param method:
        :param add_extension:
        :param extension:
        :param suffix:
        :return:
        """

        # Subdivided into methods
        if method is not None:

            # Create directory, if necessary
            if not fs.contains_directory(self.maps_sub_path, method): path = fs.create_directory_in(self.maps_sub_path, method)
            else: path = fs.join(self.maps_sub_path, method)

            # Create extra path
            extra_path = fs.create_directory_in(path, self.contour_plots_name)

            # Determine path
            if add_extension: map_path = fs.join(extra_path, name + suffix + "." + extension)
            else: map_path = fs.join(extra_path, name + suffix)

        # Determine path
        else:

            extra_path = fs.create_directory_in(self.maps_sub_path, self.contour_plots_name)
            if add_extension: map_path = fs.join(extra_path, name + suffix + "." + extension)
            else: map_path = fs.join(extra_path, name + suffix)

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def get_path_for_profile_plot(self, name, method=None, add_extension=True, extension="pdf"):

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
            extra_path = fs.create_directory_in(path, self.profile_plots_name)

            # Determine path
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Determine path
        else:

            extra_path = fs.create_directory_in(self.maps_sub_path, self.profile_plots_name)
            if add_extension: map_path = fs.join(extra_path, name + "." + extension)
            else: map_path = fs.join(extra_path, name)

        # Return the map path
        return map_path

    # -----------------------------------------------------------------

    def plot_maps(self, cmap="viridis", scale="log", format="pdf", cropping_factor=1.3, scales=None, share_limits=True,
                  mask_negatives=False):

        """
        Thisfunction ...
        :param cmap:
        :param scale:
        :param format:
        :param cropping_factor:
        :param scales:
        :param share_limits:
        :param mask_negatives:
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # The vmin and vmax
        vmin = vmax = None

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Debugging
                log.debug("Plotting maps for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting the '" + name + "' map ...")

                    # Determine path
                    plot_path = self.get_path_for_map_plot(name, method, extension=format)

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
                    if share_limits and vmin is not None and vmax is not None: interval = [vmin, vmax]
                    else: interval = "pts"
                    frame = self.maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary
                    vmin, vmax = plotting.plot_frame(frame, crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                     truncate_outside=self.truncation_ellipse, path=plot_path, format=format, interval=interval,
                                                     scale=frame_scale, cmap=cmap, normalize_in=self.truncation_ellipse, colorbar=True,
                                                     mask_negatives=mask_negatives)

                # End of method: reset vmin and vmax
                vmin = vmax = None

            # No different methods
            else:

                # Debugging
                log.debug("Plotting the '" + method + "' map ...")

                # Determine path
                plot_path = self.get_path_for_map_plot(method, extension=format)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Determine scale for this image
                if scales is not None and method in scales: frame_scale = scales[method]
                else: frame_scale = scale

                # Plot
                if share_limits and vmin is not None and vmax is not None: interval = [vmin, vmax]
                else: interval = "pts"
                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary
                vmin, vmax = plotting.plot_frame(frame, crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                 truncate_outside=self.truncation_ellipse, path=plot_path, format=format,
                                                 interval=interval, scale=frame_scale, cmap=cmap, normalize_in=self.truncation_ellipse,
                                                 colorbar=True, mask_negatives=mask_negatives)

    # -----------------------------------------------------------------

    def get_path_for_extra_map_plot(self, name, method=None, add_extension=True, extension="fits"):

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

            # Create plot path
            extra_plot_path = fs.create_directory_in(extra_path, "plots")

            # Determine path
            if add_extension: plot_path = fs.join(extra_plot_path, name + "." + extension)
            else: plot_path = fs.join(extra_plot_path, name)

        # Determine path
        else:

            # Get paths
            extra_path = fs.create_directory_in(self.maps_sub_path, self.extra_maps_name)
            extra_plot_path = fs.create_directory_in(extra_path, "plots")

            # Determine plot file path
            if add_extension: plot_path = fs.join(extra_plot_path, name + "." + extension)
            else: plot_path = fs.join(extra_plot_path, name)

        # Return the plot file path
        return plot_path

    # -----------------------------------------------------------------

    def plot_extra_maps(self, cmap="viridis", scale="log", format="pdf", cropping_factor=1.3, scales=None, share_limits=True, mask_negatives=False):

        """
        This function ...
        :param cmap:
        :param scale:
        :param format:
        :param cropping_factor:
        :param scales:
        :param share_limits:
        :param mask_negatives:
        :return:
        """

        # Inform the user
        log.info("Plotting the extra maps ...")

        # The vmin and vmax
        vmin = vmax = None

        # Loop over the methods
        for method in self.extra_maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.extra_maps[method]):

                # Debugging
                log.debug("Plotting extra maps for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.extra_maps[method]:

                    # Debugging
                    log.debug("Plotting the '" + name + "' extra map ...")

                    # Determine path
                    plot_path = self.get_path_for_extra_map_plot(name, method, extension=format)

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
                    if share_limits and vmin is not None and vmax is not None: interval = [vmin, vmax]
                    else: interval = "pts"
                    frame = self.extra_maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary
                    vmin, vmax = plotting.plot_frame(frame, crop_to=self.truncation_box,
                                                     cropping_factor=cropping_factor,
                                                     truncate_outside=self.truncation_ellipse, path=plot_path,
                                                     format=format, interval=interval,
                                                     scale=frame_scale, cmap=cmap, normalize_in=self.truncation_ellipse,
                                                     colorbar=True, mask_negatives=mask_negatives)

                # End of method: reset vmin and vmax
                vmin = vmax = None

            # No different methods
            else:

                # Debugging
                log.debug("Plotting the '" + method + "' extra map ...")

                # Determine path
                plot_path = self.get_path_for_extra_map_plot(method, extension=format)

                # If the plot already exists and we don't have to replot
                if fs.is_file(plot_path) and not self.config.replot: continue

                # Determine scale for this image
                if scales is not None and method in scales: frame_scale = scales[method]
                else: frame_scale = scale

                # Plot
                if share_limits and vmin is not None and vmax is not None: interval = [vmin, vmax]
                else: interval = "pts"
                frame = self.extra_maps[method]
                if isinstance(frame, Image): frame = frame.primary
                vmin, vmax = plotting.plot_frame(frame, crop_to=self.truncation_box,
                                                 cropping_factor=cropping_factor,
                                                 truncate_outside=self.truncation_ellipse, path=plot_path,
                                                 format=format,
                                                 interval=interval, scale=frame_scale, cmap=cmap,
                                                 normalize_in=self.truncation_ellipse, colorbar=True,
                                                 mask_negatives=mask_negatives)

    # -----------------------------------------------------------------

    def plot_contours(self, filled=False, nlevels=10, format="pdf", cropping_factor=1.3):

        """
        This function ...
        :param filled:
        :param nlevels:
        :param format:
        :param cropping_factor:
        :return:
        """

        # Inform the user
        log.info("Plotting the contours of the maps ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Debugging
                log.debug("Plotting contours for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting contours for the '" + name + "' map ...")

                    # Determine the path
                    plot_path = self.get_path_for_contour_plot(name, method, extension=format)

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
                                                            truncate_outside=self.truncation_ellipse)

                    if not filled: continue
                    filled_plot_path = self.get_path_for_contour_plot(name, method, extension=format, suffix="_filled")
                    if not fs.is_file(filled_plot_path) or self.config.replot:

                        plotting.plot_filled_frame_contours(frame, path=filled_plot_path, nlevels=nlevels,
                                                                        crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                                        truncate_outside=self.truncation_ellipse)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting contours for the '" + method + "' map ...")

                # Determine the path
                plot_path = self.get_path_for_contour_plot(method, extension=format)

                # If plot already exists and we don't have to remake
                #if fs.is_file(plot_path) and not self.config.replot: continue

                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary

                # Plot
                # if filled: plotting.plot_filled_frame_contours(self.maps[method], path=plot_path, nlevels=nlevels,
                #                                                crop_to=self.truncation_box, cropping_factor=cropping_factor,
                #                                                truncate_outside=self.truncation_ellipse)
                # else: plotting.plot_frame_contours(self.maps[method], path=plot_path, nlevels=nlevels,
                #                                    crop_to=self.truncation_box, cropping_factor=cropping_factor,
                #                                    truncate_outside=self.truncation_ellipse)

                if not fs.is_file(plot_path) or self.config.replot:

                    plotting.plot_frame_contours(frame, path=plot_path, nlevels=nlevels,
                                                       crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                        truncate_outside=self.truncation_ellipse)

                if not filled: continue
                filled_plot_path = self.get_path_for_contour_plot(method, extension=format, suffix="_filled")
                if not fs.is_file(filled_plot_path) or self.config.replot:

                    plotting.plot_filled_frame_contours(frame, path=filled_plot_path, nlevels=nlevels,
                                                                    crop_to=self.truncation_box, cropping_factor=cropping_factor,
                                                                    truncate_outside=self.truncation_ellipse)

    # -----------------------------------------------------------------

    def plot_profiles(self, nbins=20, format="pdf"):

        """
        Thisn function ...
        :param nbins:
        :param format:
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

                # Debugging
                log.debug("Plotting radial profiles for the '" + method + "' method ...")

                # Loop over the maps
                for name in self.maps[method]:

                    # Debugging
                    log.debug("Plotting radial profiles for the '" + name + "' map ...")

                    # Determine the path
                    plot_path = self.get_path_for_profile_plot(name, method, extension=format)

                    # If plot already exists and we don't have to replot
                    if fs.is_file(plot_path) and not self.config.replot: continue

                    frame = self.maps[method][name]
                    if isinstance(frame, Image): frame = frame.primary

                    # Plot
                    plotting.plot_radial_profile(frame, self.galaxy_center, angle, ratio, nbins=nbins, path=plot_path, max_radius=self.truncation_ellipse.semimajor)

            # No different methods
            else:

                # Debugging
                log.debug("Plotting radial profiles for the '" + method + "' map ...")

                # Determine the path
                plot_path = self.get_path_for_profile_plot(method, extension=format)

                # If plot already exist and we don't have to remake
                if fs.is_file(plot_path) and not self.config.replot: continue

                frame = self.maps[method]
                if isinstance(frame, Image): frame = frame.primary

                # Plot
                plotting.plot_radial_profile(frame, self.galaxy_center, angle, ratio, nbins=nbins, path=plot_path, max_radius=self.truncation_ellipse.semimajor)

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

                if method not in self.methods: raise ValueError(
                    "'" + method + "' section of the methods is missing")

                # Loop over the maps
                for name in self.maps[method]:
                    if name not in self.methods[method]: raise ValueError(
                        "Method for '" + method + "/" + name + "' map is not defined")

            # No different methods
            else:
                if method not in self.methods: raise ValueError("Method for '" + method + "' map is not defined")

    # -----------------------------------------------------------------

    @abstractmethod
    def load_collection(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def collection(self):

        """
        This function ...
        :return:
        """

        return self.load_collection()

    # -----------------------------------------------------------------

    @abstractmethod
    def load_static_collection(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def static_collection(self):

        """
        This function ...
        :return:
        """

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

    def get_ssfr_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_ssfr_origins(flatten=flatten)

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

    def get_ssfr_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.collection.get_ssfr_methods(flatten=flatten)

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

        """
        This function ...
        :return:
        """

        return self.collection.maps_colours_path

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_colours_name

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_path(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_ssfr_path

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_ssfr_name

    # -----------------------------------------------------------------

    @property
    def maps_tir_path(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_tir_path

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):

        """
        THis function ...
        :return:
        """

        return self.collection.maps_tir_name

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_path(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_attenuation_path

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_attenuation_name

    # -----------------------------------------------------------------

    @property
    def maps_old_path(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_old_path

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_old_name

    # -----------------------------------------------------------------

    @property
    def maps_young_path(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_young_path

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_young_name

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_path(self):

        """
        This fucntion ...
        :return:
        """

        return self.collection.maps_ionizing_path

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):

        """
        THis function ...
        :return:
        """

        return self.collection.maps_ionizing_name

    # -----------------------------------------------------------------

    @property
    def maps_dust_path(self):

        """
        This fucntion ...
        :return:
        """

        return self.collection.maps_dust_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):

        """
        THis function ...
        :return:
        """

        return self.collection.maps_dust_name

    # -----------------------------------------------------------------

    @property
    def maps_sub_paths(self):

        """
        This function ...
        :return:
        """

        return self.collection.maps_sub_paths

    # -----------------------------------------------------------------

    @property
    def maps_sub_names(self):

        """
        This function ...
        :return:
        """

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

class MapsComponent(MapMakerBase):
    
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
        super(MapsComponent, self).setup(**kwargs)

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

    def load_collection(self):

        """
        This function ...
        :return:
        """

        return MapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    def load_static_collection(self):

        """
        Thisf unction ...
        :return:
        """

        return StaticMapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def selection(self):

        """
        This function ...
        :return:
        """

        return ComponentMapsSelection(self.maps_components_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_selection(self):

        """
        This function ...
        :return:
        """

        return StaticComponentMapsSelection(self.maps_components_path)

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
