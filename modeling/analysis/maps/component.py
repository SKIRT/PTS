#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.component Contains the MapsAnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.basics.log import log
from ...maps.component import MapMakerBase
from ....core.tools import filesystem as fs
from ...maps.collection import get_map_paths_in_sub_path
from ....core.tools import types
from ....magic.core.frame import Frame
from ....magic.core.list import NamedFrameList
from ...maps.collection import MapsCollection, StaticMapsCollection

# -----------------------------------------------------------------

class MapsAnalysisComponent(AnalysisComponent, MapMakerBase):
    
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
        #super(MapsAnalysisComponent, self).__init__(*args, **kwargs)
        AnalysisComponent.__init__(self, no_config=True)
        MapMakerBase.__init__(self, *args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MapsAnalysisComponent, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def simulated_dataset(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulated_dataset

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.analysis_run.get_simulated_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    @property
    def colours_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.colour_maps_path

    # -----------------------------------------------------------------

    @property
    def ssfr_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.ssfr_maps_path

    # -----------------------------------------------------------------

    @property
    def tir_path(self):

        """
        Thisnf unction ...
        :return:
        """

        return self.analysis_run.tir_maps_path

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):

        """
        Thisn function ...
        :return:
        """

        return self.analysis_run.attenuation_maps_path

    # -----------------------------------------------------------------

    @property
    def old_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.old_maps_path

    # -----------------------------------------------------------------

    @property
    def young_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.young_maps_path

    # -----------------------------------------------------------------

    @property
    def ionizing_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.ionizing_maps_path

    # -----------------------------------------------------------------

    @property
    def dust_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.dust_maps_path

    # -----------------------------------------------------------------

    def load_collection(self):

        """
        This function ...
        :return:
        """

        return MapsCollection.from_analysis_run(self.analysis_run)

    # -----------------------------------------------------------------

    def load_static_collection(self):

        """
        This function ...
        :return:
        """

        return StaticMapsCollection.from_analysis_run(self.analysis_run)

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

        return get_maps_sub_name(self.analysis_run.maps_path, name, flatten=flatten, framelist=framelist, method=method)

# -----------------------------------------------------------------

def get_maps_sub_name(maps_path, name, flatten=False, framelist=False, method=None, not_method=None, not_methods=None):

    """
    This function ...
    :param maps_path:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Initialize the maps dictionary
    maps = dict()

    # Get map paths
    paths = get_map_paths_sub_name(maps_path, name, flatten=flatten, method=method, not_method=not_method, not_methods=not_methods)

    # Loop over the entries
    for method_or_name in paths:

        # Methods
        if types.is_dictionary(paths[method_or_name]):

            method_name = method_or_name
            maps[method_name] = dict()

            # Loop over the paths, load the maps and add to dictionary
            for name in paths[method_or_name]:

                map_path = paths[method_or_name][name]
                maps[method_name][name] = Frame.from_file(map_path)

        # Just maps
        elif types.is_string_type(paths[method_or_name]):

            name = method_or_name
            map_path = paths[method_or_name]
            maps[name] = Frame.from_file(map_path)

        # Something wrong
        else: raise RuntimeError("Something went wrong")

    # Return the maps
    if framelist: return NamedFrameList(**maps)
    else: return maps

# -----------------------------------------------------------------

def get_map_paths_sub_name(maps_path, name, flatten=False, method=None, not_method=None, not_methods=None):

    """
    This function ...
    :param maps_path:
    :param name:
    :param flatten:
    :param method:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(maps_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Get map paths
    return get_map_paths_in_sub_path(sub_path, flatten=flatten, method=method, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------
