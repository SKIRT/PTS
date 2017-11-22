#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.collection Contains the MapsCollection class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools.serialization import load_dict
from ..core.environment import GalaxyModelingEnvironment
from ..core.history import ModelingHistory
from ...core.tools import sequences
from ...magic.core.list import NamedFrameList
from ...magic.core.frame import Frame
from ...core.filter.filter import parse_filter
from ...core.tools.stringify import tostr
from ...core.tools import types
from ..core.environment import colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ...core.tools.utils import lazyproperty
from ...core.tools.utils import create_lazified_class
from ...magic.core.fits import get_mask_names
from ...magic.core.mask import Mask
from ...core.basics.configuration import prompt_string
from ...magic.core.image import Image

# -----------------------------------------------------------------

extra_maps_directory_names = ["extra", "TIRtoFUV", "HalphaToHot"]
plot_directory_names = ["plots", "negatives", "nans", "contours", "profiles"]

# Directories that should not be recognized as different methods
directory_names_not_methods = extra_maps_directory_names + plot_directory_names

# -----------------------------------------------------------------

origins_filename = "origins.txt"
methods_filename = "methods.txt"

# -----------------------------------------------------------------

maps_commands = ["make_colours_maps", "make_ssfr_maps", "make_tir_maps", "make_attenuation_maps", "make_old_stellar_maps", "make_dust_map", "make_young_stellar_maps", "make_ionizing_stellar_maps"]

# -----------------------------------------------------------------

def sub_name_for_command(command):

    """
    This function ...
    :param command:
    :return:
    """

    if command == "make_colours_maps": return colours_name
    elif command == "make_ssfr_maps": return ssfr_name
    elif command == "make_tir_maps": return tir_name
    elif command == "make_attenuation_maps": return attenuation_name
    elif command == "make_old_stellar_maps": return old_name
    elif command == "make_young_stellar_maps": return young_name
    elif command == "make_ionizing_stellar_maps": return ionizing_name
    elif command == "make_dust_map": return dust_name
    else: raise ValueError("Invalid commands: " + command)

# -----------------------------------------------------------------

def command_for_sub_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    if name == colours_name: return "make_colours_maps"
    elif name == ssfr_name: return "make_ssfr_maps"
    elif name == tir_name: return "make_tir_maps"
    elif name == attenuation_name: return "make_attenuation_maps"
    elif name == old_name: return "make_old_stellar_maps"
    elif name == young_name: return "make_young_stellar_maps"
    elif name == ionizing_name: return "make_ionizing_stellar_maps"
    elif name == dust_name: return "make_dust_map"
    else: raise ValueError("Invalid sub name: " + name)

# -----------------------------------------------------------------

class MapsCollection(object):

    """
    This class...
    """

    def __init__(self, maps_path, analysis_run_name=None):

        """
        The constructor ...
        :param maps_path:
        :param analysis_run_name:
        :return:
        """

        # Check whether path exists
        if not fs.is_directory(maps_path): raise IOError("The maps directory '" + maps_path + "' does not exist")

        # Set the maps path
        self.maps_path = maps_path

        # The analysis run
        self.analysis_run_name = analysis_run_name

    # -----------------------------------------------------------------

    @classmethod
    def from_modeling_path(cls, modeling_path, analysis_run_name=None):

        """
        This function ...
        :param modeling_path:
        :param analysis_run_name:
        :return:
        """

        # Set the path
        if analysis_run_name is not None: maps_path = fs.join(modeling_path, "analysis", analysis_run_name, "maps")
        else: maps_path = fs.join(modeling_path, "maps", "raw")

        # Create and return the collection
        return cls(maps_path, analysis_run_name=analysis_run_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_analysis_run(cls, run):

        """
        This function ...
        :param run:
        :return:
        """

        # Create and return the collection
        return cls.from_modeling_path(run.modeling_path, analysis_run_name=run.name)

    # -----------------------------------------------------------------

    @property
    def from_analysis(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run_name is not None

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return fs.directory_of(fs.directory_of(fs.directory_of(self.maps_path)))
        else: return fs.directory_of(fs.directory_of(self.maps_path))

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

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        from ..analysis.run import AnalysisRun
        return AnalysisRun.from_name(self.modeling_path, self.analysis_run_name)

    # -----------------------------------------------------------------

    @property
    def maps_colours_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.colour_maps_path
        else: return self.environment.maps_colours_path

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.colour_maps_name
        else: return self.environment.maps_colours_name

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.ssfr_maps_path
        else: return self.environment.maps_ssfr_path

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.ssfr_maps_name
        else: return self.environment.maps_ssfr_name

    # -----------------------------------------------------------------

    @property
    def maps_tir_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.tir_maps_path
        else: return self.environment.maps_tir_path

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):

        """
        THis function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.tir_maps_name
        else: return self.environment.maps_tir_name

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.attenuation_maps_path
        else: return self.environment.maps_attenuation_path

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.attenuation_maps_name
        else: return self.environment.maps_attenuation_name

    # -----------------------------------------------------------------

    @property
    def maps_old_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.old_maps_path
        else: return self.environment.maps_old_path

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.old_maps_name
        else: return self.environment.maps_old_name

    # -----------------------------------------------------------------

    @property
    def maps_young_path(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.young_maps_path
        else: return self.environment.maps_young_path

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):

        """
        This function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.young_maps_name
        else: return self.environment.maps_young_name

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_path(self):

        """
        This fucntion ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.ionizing_maps_path
        else: return self.environment.maps_ionizing_path

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):

        """
        THis function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.ionizing_maps_name
        else: return self.environment.maps_ionizing_name

    # -----------------------------------------------------------------

    @property
    def maps_dust_path(self):

        """
        This fucntion ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.dust_maps_path
        else: return self.environment.maps_dust_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):

        """
        THis function ...
        :return:
        """

        if self.from_analysis: return self.analysis_run.dust_maps_name
        else: return self.environment.maps_dust_name

    # -----------------------------------------------------------------

    @property
    def maps_sub_paths(self):

        """
        This function ...
        :return:
        """

        return [self.maps_colours_path, self.maps_ssfr_path, self.maps_tir_path, self.maps_attenuation_path,
                self.maps_old_path, self.maps_young_path, self.maps_ionizing_path, self.maps_dust_path]

    # -----------------------------------------------------------------

    @property
    def maps_sub_names(self):

        """
        This function ...
        :return:
        """

        return [fs.name(path) for path in self.maps_sub_paths]

    # -----------------------------------------------------------------

    def get_origins_sub_name(self, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return: 
        """

        if self.from_analysis: return get_origins_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)
        else: return get_origins_sub_name(self.environment, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

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

    def get_all_filters_sub_name(self, name, method=None):

        """
        This function ...
        :param name:
        :param method:
        :return:
        """

        filters = set()
        origins = self.get_origins_sub_name(name, flatten=True, method=method)
        for name in origins: filters.update(origins[name])
        return list(filters)

    # -----------------------------------------------------------------

    def get_all_filters(self):

        """
        This function ...
        :return:
        """

        filters = set()
        for name in self.maps_sub_names: filters.update(self.get_all_filters_sub_name(name))
        return list(filters)

    # -----------------------------------------------------------------

    def get_methods_sub_name(self, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        if self.from_analysis: return get_methods_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)
        else: return get_methods_sub_name(self.environment, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_methods_sub_names(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        methods = dict()
        for name in self.maps_sub_names: methods[name] = self.get_methods_sub_name(name, flatten=flatten)
        return methods

    # -----------------------------------------------------------------

    # ORIGINS

    def get_colours_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return: 
        """

        return self.get_origins_sub_name(self.maps_colours_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def colour_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_colours_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def colour_origins_flat(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_colours_origins(flatten=True)

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

        return self.get_origins_sub_name(self.maps_ssfr_name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    @property
    def ssfr_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ssfr_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_tir_origins(self, flatten=False, methods=None):

        """
        This function ...
        :param flatten:
        :param methods:
        :return:
        """

        return self.get_origins_sub_name(self.maps_tir_name, flatten=flatten, methods=methods)

    # -----------------------------------------------------------------

    @property
    def tir_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def tir_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_attenuation_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_attenuation_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def attenuation_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def attenuation_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_old_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_old_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def old_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_old_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def old_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_old_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_young_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_young_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def young_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_young_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def young_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_young_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_ionizing_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_ionizing_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def ionizing_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_ionizing_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ionizing_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_ionizing_origins(flatten=True)

    # -----------------------------------------------------------------

    def get_dust_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(self.maps_dust_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def dust_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def dust_origins_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_origins(flatten=True)

    # -----------------------------------------------------------------

    # FILTERS

    def get_all_colours_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_colours_name)

    # -----------------------------------------------------------------

    def get_all_ssfr_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_ssfr_name)

    # -----------------------------------------------------------------

    def get_all_tir_filters(self):

        """
        Thisn function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_tir_name)

    # -----------------------------------------------------------------

    def get_all_attenuation_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_attenuation_name)

    # -----------------------------------------------------------------

    def get_all_old_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_old_name)

    # -----------------------------------------------------------------

    def get_all_young_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_young_name)

    # -----------------------------------------------------------------

    def get_all_ionizing_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_ionizing_name)

    # -----------------------------------------------------------------

    def get_all_dust_filters(self):

        """
        This function ...
        :return:
        """

        return self.get_all_filters_sub_name(self.maps_dust_name)

    # -----------------------------------------------------------------

    # METHODS

    def get_colours_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_colours_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def colour_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_colours_methods(flatten=False)

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

        return self.get_methods_sub_name(self.maps_ssfr_name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    @property
    def ssfr_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_tir_methods(self, flatten=False, methods=None):

        """
        This function ...
        :param flatten:
        :param methods:
        :return:
        """

        return self.get_methods_sub_name(self.maps_tir_name, flatten=flatten, methods=methods)

    # -----------------------------------------------------------------

    @property
    def tir_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_attenuation_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_attenuation_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def attenuation_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_old_methods(self, flatten=False):

        """
        Thisf unction ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_old_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def old_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_old_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_young_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_young_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def young_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_young_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_ionizing_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_ionizing_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def ionizing_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_ionizing_methods(flatten=False)

    # -----------------------------------------------------------------

    def get_dust_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(self.maps_dust_name, flatten=flatten)

    # -----------------------------------------------------------------

    @property
    def dust_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_methods(flatten=False)

    # -----------------------------------------------------------------

    # MAPS

    def get_colour_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_colours_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

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

    def get_colour_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_colours_name, flatten=flatten, framelist=framelist)

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

        return self.get_map_paths_sub_name(self.maps_ssfr_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

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

        return self.get_maps_sub_name(self.maps_ssfr_name, flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_ssfr_nans(self, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

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

        return self.get_nans_sub_name(self.maps_ssfr_name, flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

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

        return self.get_map_paths_sub_name(self.maps_tir_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_tir_maps(self, flatten=False, framelist=False, methods=None):

        """
        This function ...
        :param flatten:
        :param framelist:
        :param methods:
        :return:
        """

        return self.get_maps_sub_name(self.maps_tir_name, flatten=flatten, framelist=framelist, methods=methods)

    # -----------------------------------------------------------------

    def get_tir_nans(self, flatten=False, framelist=False, methods=None):

        """
        Thisfunction ...
        :param flatten:
        :param framelist:
        :param methods:
        :return:
        """

        return self.get_nans_sub_name(self.maps_tir_name, flatten=flatten, framelist=framelist, methods=methods)

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

        return self.get_map_paths_sub_name(self.maps_attenuation_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_attenuation_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist::
        :return:
        """

        return self.get_maps_sub_name(self.maps_attenuation_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_attenuation_extra_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_extra_maps_sub_name(self.maps_attenuation_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_attenuation_nans(self, flatten=False, framelist=False):

        """
        Thisfunction ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_attenuation_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_old_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_old_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_old_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_old_nans(self, flatten=False, framelist=False):

        """
        Thisn function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_old_name, flatten=flatten, framelist=framelist)

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

        return self.get_map_paths_sub_name(self.maps_young_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_young_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_young_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_young_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_young_name, flatten=flatten, framelist=framelist)

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

        return self.get_map_paths_sub_name(self.maps_ionizing_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_ionizing_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_ionizing_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_ionizing_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_ionizing_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_dust_map_paths(self, flatten=False, method=None, startswith=None, factors=None):

        """
        This function ...
        :param flatten:
        :param method:
        :param startswith:
        :param factors:
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_dust_name, flatten=flatten, method=method, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_dust_maps(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_dust_name, flatten=flatten, framelist=framelist)

    # -----------------------------------------------------------------

    def get_dust_negatives(self, flatten=False, framelist=False, method=None):

        """
        Thisf unction ...
        :param flatten:
        :param framelist:
        :param method:
        :return:
        """

        return self.get_negatives_sub_name(self.maps_dust_name, flatten=flatten, framelist=framelist, method=method)

    # -----------------------------------------------------------------

    def get_dust_nans(self, flatten=False, framelist=False):

        """
        This function ...
        :param flatten:
        :param framelist:
        :return:
        """

        return self.get_nans_sub_name(self.maps_dust_name, flatten=flatten, framelist=framelist)

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

    def get_fuv_attenuation_maps(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        if cortese: cortese = self.get_cortese_fuv_attenuation_maps()
        else: cortese = dict()
        if buat: buat = self.get_buat_fuv_attenuation_maps()
        else: buat = dict()

        if flatten:

            maps = dict()
            for name in cortese: maps["cortese__" + name] = cortese[name]
            for name in buat: maps["buat__" + name] = buat[name]
            return maps

        else:

            maps = dict()
            maps["cortese"] = cortese
            maps["buat"] = buat
            return maps

    # -----------------------------------------------------------------

    def get_fuv_attenuation_origins(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        if cortese: cortese = self.get_cortese_fuv_attenuation_origins()
        else: cortese = dict()
        if buat: buat = self.get_buat_fuv_attenuation_origins()
        else: buat = dict()

        if flatten:

            origins = dict()
            for name in cortese: origins["cortese__" + name] = cortese[name]
            for name in buat: origins["buat__" + name] = buat[name]
            return origins

        else:

            origins = dict()
            origins["cortese"] = cortese
            origins["buat"] = buat
            return origins

    # -----------------------------------------------------------------

    def get_fuv_attenuation_methods(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        if cortese: cortese = self.get_cortese_fuv_attenuation_methods()
        else: cortese = dict()
        if buat: buat = self.get_buat_fuv_attenuation_methods()
        else: buat = dict()

        if flatten:

            methods = dict()
            for name in cortese: methods["cortese__" + name] = cortese[name]
            for name in buat: methods["buat__" + name] = buat[name]
            return methods

        else:

            methods = dict()
            methods["cortese"] = cortese
            methods["buat"] = buat
            return methods

    # -----------------------------------------------------------------

    def get_fuv_attenuation_nans(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        if cortese: cortese = self.get_cortese_fuv_attenuation_nans()
        else: cortese = dict()
        if buat: buat = self.get_buat_fuv_attenuation_nans()
        else: buat = dict()

        if flatten:

            nans = dict()
            for name in cortese: nans["cortese__" + name] = cortese[name]
            for name in buat: nans["buat__" + name] = buat[name]
            return nans

        else:

            nans = dict()
            nans["cortese"] = cortese
            nans["buat"] = buat
            return nans

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_and_origins(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        maps = self.get_fuv_attenuation_maps(flatten=flatten, cortese=cortese, buat=buat)
        origins = self.get_fuv_attenuation_origins(flatten=flatten, cortese=cortese, buat=buat)

        # Check
        if not sequences.same_contents(maps.keys(), origins.keys()):

            log.error("Mismatch between FUV attenuation maps and their origins:")
            #raise ValueError("Mismatch between FUV attenuation maps names and their origins")

            sorted_keys_maps = sorted(maps.keys())
            sorted_keys_origins = sorted(origins.keys())

            if len(sorted_keys_maps) != len(sorted_keys_origins): log.error("Number of maps: " + str(len(sorted_keys_maps)) + " vs Number of origins: " + str(len(sorted_keys_origins)))

            indices = sequences.find_differences(sorted_keys_maps, sorted_keys_origins)

            log.error("Number of mismatches: " + str(len(indices)))

            #for index in indices: log.error(" - " + sorted_keys_maps[index] + " vs " + sorted_keys_origins[index])

            for index in range(min(len(sorted_keys_maps), len(sorted_keys_origins))):

                if sorted_keys_maps[index] == sorted_keys_origins[index]: log.success(" - " + sorted_keys_maps[index] + " = " + sorted_keys_origins[index])
                else: log.error(" - " + sorted_keys_maps[index] + " != " + sorted_keys_origins[index])

            # QUIT
            exit()

        #exit()

        # Return
        return maps, origins

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_origins_and_methods(self, flatten=False, cortese=True, buat=True):

        """
        This function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        # Already checked: maps vs. origins
        maps, origins = self.get_fuv_attenuation_maps_and_origins(flatten=flatten, cortese=cortese, buat=buat)
        methods = self.get_fuv_attenuation_methods(flatten=flatten, cortese=cortese, buat=buat)

        # Check
        if not sequences.same_contents(maps.keys(), methods.keys()):

            log.error("Mismatch between FUV attenuation maps and their methods:")

            sorted_keys_maps = sorted(maps.keys())
            sorted_keys_methods = sorted(methods.keys())

            if len(sorted_keys_maps) != len(sorted_keys_methods): log.error("Number of maps: " + str(len(sorted_keys_maps)) + " vs Number of methods: " + str(len(sorted_keys_methods)))

            indices = sequences.find_differences(sorted_keys_maps, sorted_keys_methods)

            log.error("Number of mismatches: " + str(len(indices)))

            for index in range(min(len(sorted_keys_maps), len(sorted_keys_methods))):

                if sorted_keys_maps[index] == sorted_keys_methods[index]: log.success(" - " + sorted_keys_maps[index] + " = " + sorted_keys_methods[index])
                else: log.error(" - " + sorted_keys_maps[index] + " != " + sorted_keys_methods[index])

            # QUIT
            exit()

        # Return
        return maps, origins, methods

    # -----------------------------------------------------------------

    def get_fuv_attenuation_maps_origins_methods_and_nans(self, flatten=False, cortese=False, buat=True):

        """
        Thisn function ...
        :param flatten:
        :param cortese:
        :param buat:
        :return:
        """

        # Already checked: maps, origins and methods
        maps, origins, methods = self.get_fuv_attenuation_maps_origins_and_methods(flatten=flatten, cortese=cortese, buat=buat)

        # Get nans
        nans = self.get_fuv_attenuation_nans(flatten=flatten, cortese=cortese, buat=buat)

        # Check
        if not sequences.same_contents(maps.keys(), nans.keys()):

            log.error("Mismatch between FUV attenuation maps and their nans:")

            sorted_keys_maps = sorted(maps.keys())
            sorted_keys_nans = sorted(nans.keys())

            if len(sorted_keys_maps) != len(sorted_keys_nans): log.error("Number of maps: " + str(len(sorted_keys_maps)) + " vs Number of nans: " + str(len(sorted_keys_nans)))

            indices = sequences.find_differences(sorted_keys_maps, sorted_keys_nans)

            log.error("Number of mismatches: " + str(len(indices)))

            for index in range(min(len(sorted_keys_maps), len(sorted_keys_nans))):

                if sorted_keys_maps[index] == sorted_keys_nans[index]: log.success(" - " + sorted_keys_maps[index] + " = " + sorted_keys_nans[index])
                else: log.error(" - " + sorted_keys_maps[index] + " != " + sorted_keys_nans[index])

            # QUIT
            exit()

        # Return
        return maps, origins, methods, nans

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
        origins = load_dict(table_path)

        # ONLY KEEP FUV
        for name in list(origins.keys()):
            if not name.startswith("FUV"): del origins[name]

        # Return
        return origins

    # -----------------------------------------------------------------

    def get_buat_nuv_attenuation_origins(self):

        """
        This function ...
        :return:
        """

        buat_path = fs.join(self.maps_attenuation_path, "buat")
        table_path = fs.join(buat_path, origins_filename)
        origins = load_dict(table_path)

        # ONLY KEEP NUV
        for name in list(origins.keys()):
            if not name.startswith("NUV"): del origins[name]

        # Return
        return origins

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

    def get_cortese_fuv_attenuation_methods(self):

        """
        Tihs function ...
        :return:
        """

        cortese_path = fs.join(self.maps_attenuation_path, "cortese")
        methods_path = fs.join(cortese_path, methods_filename)
        return load_dict(methods_path)

    # -----------------------------------------------------------------

    def get_cortese_fuv_attenuation_nans(self):

        """
        Thisn function ...
        :return:
        """

        cortese_path = fs.join(self.maps_attenuation_path, "cortese")

        nans = dict()

        for path, name in fs.files_in_path(cortese_path, returns=["path", "name"], extension="fits"):

            # Mask or None
            if "nans" in get_mask_names(path): nans[name] = Mask.from_file(path, plane="nans")
            else: nans[name] = None

        # Return the dictionary
        return nans

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_methods(self):

        """
        This function ...
        :return:
        """

        buat_path = fs.join(self.maps_attenuation_path, "buat")
        methods_path = fs.join(buat_path, methods_filename)
        methods = load_dict(methods_path)

        # ONLY KEEP FUV
        for name in list(methods.keys()):
            if not name.startswith("FUV"): del methods[name]

        # Return the methods
        return methods

    # -----------------------------------------------------------------

    def get_buat_fuv_attenuation_nans(self):

        """
        Thisnfunction ...
        :return:
        """

        buat_path = fs.join(self.maps_attenuation_path, "buat")

        # Initialize
        nans = dict()

        for path, name in fs.files_in_path(buat_path, returns=["path", "name"], extension="fits", contains="FUV"):

            # Mask or None
            if "nans" in get_mask_names(path): nans[name] = Mask.from_file(path, plane="nans")
            else: nans[name] = None

        # Return the dictionary
        return nans

    # -----------------------------------------------------------------

    @property
    def old_stellar_total_maps_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_old_path, "total")

    # -----------------------------------------------------------------

    @property
    def old_stellar_total_filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(name) for name in fs.files_in_path(self.old_stellar_total_maps_path, extension="fits", returns="name")]

    # -----------------------------------------------------------------

    def get_old_stellar_total_map(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        path = fs.join(self.old_stellar_total_maps_path, tostr(fltr, delimiter="_") + ".fits")
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    @property
    def old_stellar_bulge_maps_path(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.join(self.maps_old_path, "bulge")

    # -----------------------------------------------------------------

    @property
    def old_stellar_bulge_filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(name) for name in fs.files_in_path(self.old_stellar_bulge_maps_path, extension="fits", returns="name")]

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_map(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        path = fs.join(self.old_stellar_bulge_maps_path, tostr(fltr, delimiter="_") + ".fits")
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_maps(self, framelist=False):

        """
        This function ...
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_old_name, framelist=framelist, method="bulge")

    # -----------------------------------------------------------------

    @property
    def old_stellar_disk_maps_path(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.join(self.maps_old_path, "disk")

    # -----------------------------------------------------------------

    @property
    def old_stellar_disk_filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(name) for name in fs.files_in_path(self.old_stellar_disk_maps_path, extension="fits", returns="name")]

    # -----------------------------------------------------------------

    def get_old_stellar_disk_map(self, fltr):

        """
        This fucntion ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        path = fs.join(self.old_stellar_disk_maps_path, tostr(fltr, delimiter="_") + ".fits")
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    def get_old_stellar_disk_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_map_paths_sub_name(self.maps_old_name, method="disk")

    # -----------------------------------------------------------------

    def get_old_stellar_disk_maps(self, framelist=False):

        """
        This function ...
        :param framelist:
        :return:
        """

        return self.get_maps_sub_name(self.maps_old_name, framelist=framelist, method="disk")

    # -----------------------------------------------------------------

    def get_old_stellar_total_origins(self):

        """
        Thisn function ...
        :return:
        """

        return self.get_origins_sub_name(self.maps_old_name, method="total")

    # -----------------------------------------------------------------

    def get_old_stellar_disk_origins(self):

        """
        This function ...
        :return:
        """

        return self.get_origins_sub_name(self.maps_old_name, method="disk")

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_origins(self):

        """
        Thisn function ...
        :return:
        """

        return self.get_origins_sub_name(self.maps_old_name, method="bulge")

    # -----------------------------------------------------------------

    def get_old_stellar_total_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_methods_sub_name(self.maps_old_name, method="total")

    # -----------------------------------------------------------------

    def get_old_stellar_disk_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_methods_sub_name(self.maps_old_name, method="disk")

    # -----------------------------------------------------------------

    def get_old_stellar_bulge_methods(self):

        """
        This function ...
        :return:
        """

        return self.get_methods_sub_name(self.maps_old_name, method="bulge")

    # -----------------------------------------------------------------

    def get_hot_dust_maps(self):

        """
        This function ...
        :return:
        """

        hot_dust_path = fs.join(self.maps_dust_path, "hot")
        return NamedFrameList.from_directory(hot_dust_path).to_dictionary()

    # -----------------------------------------------------------------

    def get_hot_dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_map_paths_sub_name(dust_name, method="hot")

    # -----------------------------------------------------------------

    def get_not_hot_dust_maps(self, flatten=False):

        """
        This function ...
        :return:
        """

        return self.get_maps_sub_name(dust_name, not_method="hot", flatten=flatten)

    # -----------------------------------------------------------------

    def get_not_hot_dust_map_paths(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_map_paths_sub_name(dust_name, not_method="hot", flatten=flatten)

    # -----------------------------------------------------------------

    def get_not_hot_dust_origins(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_origins_sub_name(dust_name, not_method="hot", flatten=flatten)

    # -----------------------------------------------------------------

    def get_not_hot_dust_methods(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        return self.get_methods_sub_name(dust_name, not_method="hot", flatten=flatten)

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

    def get_hot_dust_methods(self):

        """
        This function ...
        :return:
        """

        hot_dust_path = fs.join(self.maps_dust_path, "hot")
        methods_path = fs.join(hot_dust_path, methods_filename)
        return load_dict(methods_path)

    # -----------------------------------------------------------------

    def get_map_paths_sub_name(self, name, flatten=False, method=None, not_method=None, not_methods=None, startswith=None, factors=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param method:
        :param not_method:
        :param not_methods:
        :param startswith:
        :param factors:
        :return:
        """

        if self.from_analysis: return get_map_paths_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, not_method=not_method, not_methods=not_methods, startswith=startswith, factors=factors)
        else: return get_map_paths_sub_name(self.environment, name, flatten=flatten, method=method, not_method=not_method, not_methods=not_methods, startswith=startswith, factors=factors)

    # -----------------------------------------------------------------

    def get_maps_sub_name(self, name, flatten=False, framelist=False, method=None, methods=None, not_method=None,
                          not_methods=None, factors=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :param factors:
        :return:
        """

        if self.from_analysis: return get_maps_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods,
                                                                 not_method=not_method, not_methods=not_methods, factors=factors)
        else: return get_maps_sub_name(self.environment, self.history, name, flatten=flatten, framelist=framelist, method=method,
                                       methods=methods, not_method=not_method, not_methods=not_methods, factors=factors)

    # -----------------------------------------------------------------

    def get_extra_maps_sub_name(self, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        Thisf unction ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        if self.from_analysis: return get_extra_maps_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)
        else: return get_extra_maps_sub_name(self.environment, name, flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_negatives_sub_name(self, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        if self.from_analysis: return get_negatives_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)
        else: return get_negatives_sub_name(self.environment, name, flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_nans_sub_name(self, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

        """
        This function ...
        :param name:
        :param flatten:
        :param framelist:
        :param method:
        :param methods:
        :param not_method:
        :param not_methods:
        :return:
        """

        if self.from_analysis: return get_nans_sub_name_analysis(self.analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)
        else: return get_nans_sub_name(self.environment, name, flatten=flatten, framelist=framelist, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # -----------------------------------------------------------------

    def get_maps_sub_names(self, flatten=False):

        """
        This function ...
        :param flatten:
        :return:
        """

        maps = dict()
        for name in self.maps_sub_names: maps[name] = self.get_maps_sub_name(name, flatten=flatten)
        return maps

    # -----------------------------------------------------------------

    def get_nmaps_sub_name(self, name, method=None):

        """
        This function ...
        :param name:
        :param method:
        :return:
        """

        paths = self.get_map_paths_sub_name(name, flatten=True, method=method)
        return len(paths)

    # -----------------------------------------------------------------

    def has_maps_sub_name(self, name, method=None):

        """
        This function ...
        :param name:
        :param method:
        :return:
        """

        return self.get_nmaps_sub_name(name, method=method) > 0

    # -----------------------------------------------------------------

    @property
    def ncolour_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(colours_name)
        #return len(self.colour_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.colour_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def nssfr_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(ssfr_name)
        #return len(self.ssfr_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.ssfr_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def ntir_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(tir_name)
        #return len(self.tir_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.tir_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def nattenuation_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(attenuation_name)
        #return len(self.attenuation_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.attenuation_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def nold_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(old_name)
        #return len(self.old_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.old_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def nyoung_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(young_name)
        #return len(self.young_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.young_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def nionizing_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(ionizing_name)
        #return len(self.ionizing_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.ionizing_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def ndust_maps(self):

        """
        This function ...
        :return:
        """

        #return self.get_nmaps_sub_name(dust_name)
        #return len(self.dust_maps_flat) # SO IT CAN GET LAZIFIED
        return len(self.dust_map_paths_flat)

    # -----------------------------------------------------------------

    @property
    def has_colour_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(colours_name)
        return self.ncolour_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(ssfr_name)
        return self.nssfr_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_tir_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(tir_name)
        return self.ntir_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(attenuation_name)
        return self.nattenuation_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_old_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(old_name)
        return self.nold_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_young_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(young_name)
        return self.nyoung_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(ionizing_name)
        return self.nionizing_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def has_dust_maps(self):

        """
        This function ...
        :return:
        """

        #return self.has_maps_sub_name(dust_name)
        return self.ndust_maps > 0 # SO IT CAN GET LAZIFIED

    # -----------------------------------------------------------------

    @property
    def colour_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_colour_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def colour_map_paths_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_colour_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def colour_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.colour_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def colour_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.colour_has_methods: raise ValueError("No methods")
        else: return self.colour_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def colour_map_names(self):

        """
        This ufnction ...
        :return:
        """

        if self.colour_has_methods: raise ValueError("Has different methods. Call colour_map_names_for_method.")
        else: return self.colour_map_paths.keys()

    # -----------------------------------------------------------------

    def colour_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.colour_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def colour_map_names_flat(self):

        """
        This function ...
        :return:
        """

        return self.colour_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def colour_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_colour_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def colour_maps_flat(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_colour_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ssfr_map_paths_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.ssfr_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.ssfr_has_methods: raise ValueError("No methods")
        else: return self.ssfr_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def ssfr_map_names(self):

        """
        This function ...
        :return:
        """

        if self.ssfr_has_methods: raise ValueError("Has different methods. Call ssfr_map_names_for_method.")
        else: return self.ssfr_map_paths.keys()

    # -----------------------------------------------------------------

    def ssfr_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.ssfr_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def ssfr_map_names_flat(self):

        """
        Thisfunction ...
        :return:
        """

        return self.ssfr_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ssfr_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_ssfr_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def tir_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def tir_map_paths_flat(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_tir_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def tir_has_methods(self):

        """
        Thisf unction ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.tir_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def tir_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.tir_has_methods: raise ValueError("No methods")
        else: return self.tir_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def tir_map_names(self):

        """
        This function ...
        :return:
        """

        if self.tir_has_methods: raise ValueError("Has methods. Call tir_map_names_for_method.")
        else: return self.tir_map_paths.keys()

    # -----------------------------------------------------------------

    def tir_map_names_for_method(self, method):

        """
        This function ....
        :param method:
        :return:
        """

        return self.tir_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def tir_map_names_flat(self):

        """
        This function ...
        :return:
        """

        return self.tir_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def tir_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def tir_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_tir_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def attenuation_map_paths_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.attenuation_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.attenuation_has_methods: raise ValueError("No methods")
        else: return self.attenuation_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def attenuation_map_names(self):

        """
        This function ...
        :return:
        """

        if self.attenuation_has_methods: raise ValueError("Has methods. Call attenuation_map_names_for_method.")
        else: return self.attenuation_map_paths.keys()

    # -----------------------------------------------------------------

    def attenuation_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.attenuation_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def attenuation_map_names_flat(self):

        """
        This funtion ...
        :return:
        """

        return self.attenuation_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def attenuation_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_attenuation_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def old_map_paths(self):

        """
        Thisfunction ...
        :return:
        """

        return self.get_old_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def old_map_paths_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_old_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def old_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.old_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def old_map_methods(self):

        """
        This function ....
        :return:
        """

        if not self.old_has_methods: raise ValueError("No methods")
        else: return self.old_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def old_map_names(self):

        """
        This ufnction ...
        :return:
        """

        if self.old_has_methods: raise ValueError("Has methods. Call old_map_names_for_method.")
        else: return self.old_map_paths.keys()

    # -----------------------------------------------------------------

    def old_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.old_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def old_map_names_flat(self):

        """
        Thisf unction ...
        :return:
        """

        return self.old_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def old_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_old_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def old_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_old_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def young_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_young_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def young_map_paths_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_young_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def young_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.young_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def young_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.young_has_methods: raise ValueError("No methods")
        else: return self.young_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def young_map_names(self):

        """
        This function ...
        :return:
        """

        if self.young_has_methods: raise ValueError("Has methods. Call young_map_names_for_method.")
        else: return self.young_map_paths.keys()

    # -----------------------------------------------------------------

    def young_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.young_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def young_map_names_flat(self):

        """
        Thisfunction ...
        :return:
        """

        return self.young_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def young_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_young_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def young_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_young_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_paths(self):

        """
        Thisfunction ...
        :return:
        """

        return self.get_ionizing_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_paths_flat(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_ionizing_map_paths(flatten=True)

    # ---------------------------------------unc--------------------------

    @property
    def ionizing_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.ionizing_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_methods(self):

        """
        Thisj function ....
        :return:
        """

        if not self.ionizing_has_methods: raise ValueError("No methods")
        else: return self.ionizing_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names(self):

        """
        This function ...
        :return:
        """

        if self.ionizing_has_methods: raise ValueError("Has methods. Call ionizing_map_names_for_method.")
        else: return self.ionizing_map_paths.keys()

    # -----------------------------------------------------------------

    def ionizing_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.ionizing_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names_flat(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def ionizing_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_ionizing_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ionizing_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_ionizing_maps(flatten=True)

    # -----------------------------------------------------------------

    @property
    def dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_map_paths(flatten=False)

    # -----------------------------------------------------------------

    @property
    def dust_map_paths_flat(self):

        """
        This funciton ...
        :return:
        """

        return self.get_dust_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def dust_has_methods(self):

        """
        This function ...
        :return:
        """

        return types.is_dictionary_of_dictionaries(self.dust_map_paths, passive=True)

    # -----------------------------------------------------------------

    @property
    def dust_map_methods(self):

        """
        This function ...
        :return:
        """

        if not self.dust_has_methods: raise ValueError("No methods")
        else: return self.dust_map_paths.keys()

    # -----------------------------------------------------------------

    @property
    def dust_map_names(self):

        """
        This function ...
        :return:
        """

        if self.dust_has_methods: raise ValueError("Has methods. Call dust_map_names_for_method")
        else: return self.dust_map_paths.keys()

    # -----------------------------------------------------------------

    def dust_map_names_for_method(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        return self.dust_map_paths[method].keys()

    # -----------------------------------------------------------------

    @property
    def dust_map_names_flat(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_paths_flat.keys()

    # -----------------------------------------------------------------

    @property
    def dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_maps(flatten=False)

    # -----------------------------------------------------------------

    @property
    def dust_maps_flat(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_maps(flatten=True)

# -----------------------------------------------------------------

StaticMapsCollection = create_lazified_class(MapsCollection, "StaticMapsCollection")

# -----------------------------------------------------------------

def get_map_paths_sub_name_analysis(analysis_run, name, flatten=False, method=None, methods=None, not_method=None,
                                    not_methods=None, factors=None, startswith=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :param factors:
    :param startswith:
    :return:
    """

    # Determine path
    sub_path = fs.join(analysis_run.maps_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Get map paths
    return get_map_paths_in_sub_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method,
                                     not_methods=not_methods, factors=factors, startswith=startswith)

# -----------------------------------------------------------------

def get_extra_map_paths_sub_name_analysis(analysis_run, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(analysis_run.maps_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Get map paths
    return get_extra_map_paths_in_sub_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_map_paths_sub_name(environment, name, flatten=False, method=None, methods=None, not_method=None,
                           not_methods=None, factors=None, startswith=None):

    """
    This function ...
    :param environment:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :param factors:
    :param startswith:
    :return:
    """

    # Determine path
    sub_path = fs.join(environment.maps_raw_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Get map paths
    return get_map_paths_in_sub_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method,
                                     not_methods=not_methods, factors=factors, startswith=startswith)

# -----------------------------------------------------------------

def get_extra_map_paths_sub_name(environment, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    Thisf unction ...
    :param environment:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(environment.maps_raw_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Get map paths
    return get_extra_map_paths_in_sub_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_map_paths_in_sub_path(sub_path, flatten=False, method=None, methods=None, not_method=None, not_methods=None, factors=None, startswith=None):

    """
    This function ...
    :param sub_path:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :param factors:
    :param startswith:
    :return:
    """

    if factors is not None: factor_strings = ["__" + repr(factor) for factor in factors]
    else: factor_strings = None

    # direct_origins_path = fs.join(sub_path, origins_filename)
    # No subdirectories
    # if fs.is_file(direct_origins_path): origins = load_dict(direct_origins_path)

    if not_method is not None:
        if not types.is_string_type(not_method): raise ValueError("'not_method' must be a string")
        if not_methods is not None: raise ValueError("Cannot specifify both 'not_method' and 'not_methods'")
        not_methods = [not_method]

    # Subdirectories
    if fs.contains_directories(sub_path, exact_not_name=directory_names_not_methods):

        # One method is specified
        if method is not None:

            # Check
            if not_methods is not None: raise ValueError("Cannot specify both 'method' and 'not_methods' simultaneously")

            # Check whether valid method
            method_path = fs.join(sub_path, method)
            if not fs.is_directory(method_path): raise ValueError("Directory not found for method '" + method + "'")

            # Return the file paths
            return fs.files_in_path(method_path, returns="dict", extension="fits", endswith=factor_strings, startswith=startswith)

        # Method not specified
        else:

            paths = dict()

            # Loop over the subdirectories
            for method_path, method_name in fs.directories_in_path(sub_path, returns=["path", "name"]):

                # Skip other method if method is defined
                if method is not None and method_name != method: continue

                # Skip method if not in list
                if methods is not None and method_name not in methods: continue

                # Skip methods in not_method
                if not_methods is not None and method_name in not_methods: continue

                # Get dictionary of file paths, but only FITS files
                files = fs.files_in_path(method_path, returns="dict", extension="fits", endswith=factor_strings, startswith=startswith)

                # Flatten into a one-level dict
                if flatten:
                    for map_name in files: paths[method_name + "_" + map_name] = files[map_name]

                # Don't flatten
                else: paths[method_name] = files

            # Return the paths
            return paths

    # Files present
    elif fs.contains_files(sub_path):

        # Method cannot be defined
        if method is not None: raise ValueError("Specified method '" + method + "', but all maps are in one directory")

        # Methods cannot be defined
        if methods is not None: raise ValueError("Specified methods '" + str(methods) + "', but all maps are in one directory")

        # Not method cannot be defined
        if not_methods is not None: raise ValueError("All maps are in one directory (no different methods)")

        # Return the file paths
        return fs.files_in_path(sub_path, returns="dict", extension="fits", endswith=factor_strings, startswith=startswith)

    # Nothing present
    else: return dict()

# -----------------------------------------------------------------

def get_extra_map_paths_in_sub_path(sub_path, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param sub_path:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    if not_method is not None:
        if not types.is_string_type(not_method): raise ValueError("'not_method' must be a string")
        if not_methods is not None: raise ValueError("Cannot specifify both 'not_method' and 'not_methods'")
        not_methods = [not_method]

    # Subdirectories
    if fs.contains_directories(sub_path, exact_not_name=directory_names_not_methods):

        # One method is specified
        if method is not None:

            # Check
            if not_methods is not None: raise ValueError("Cannot specify both 'method' and 'not_methods' simultaneously")

            # Check whether valid method
            method_path = fs.join(sub_path, method)
            if not fs.is_directory(method_path): raise ValueError("Directory not found for method '" + method + "'")

            # Find a directory that contains extra maps (a directorty that isn't called 'plot', 'contours', or 'profiles'
            subdirectory_names = fs.directories_in_path(method_path, exact_not_name=directory_names_not_methods, returns="name")

            # No extra maps
            if len(subdirectory_names) == 0: return dict()

            # One directory with extra maps
            elif len(subdirectory_names) == 1:

                # Get the path to the extra maps directory
                extra_maps_path = fs.join(method_path, subdirectory_names[0])

                # Return the FITS file paths
                return fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

            # More?
            else:

                # GIVE WARNING
                log.warning("Please remove directories other than 'plots', 'contours', and 'profiles'")

                # Get the name and path
                extra_maps_name = prompt_string("extra_directory_name", "directory with the extra maps in '" + fs.name(sub_path) + "/" + method, choices=subdirectory_names, required=False)
                if extra_maps_name is None: return dict() # none of the directories
                extra_maps_path = fs.join(method_path, extra_maps_name)

                # Return the FITS file paths
                return fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

        # Method not specified
        else:

            paths = dict()

            # Loop over the subdirectories
            for method_path, method_name in fs.directories_in_path(sub_path, returns=["path", "name"]):

                # Skip other method if method is defined
                if method is not None and method_name != method: continue

                # Skip method if not in list
                if methods is not None and method_name not in methods: continue

                # Skip methods in not_method
                if not_methods is not None and method_name in not_methods: continue

                # Find a directory that contains extra maps (a directorty that isn't called 'plot', 'contours', or 'profiles'
                subdirectory_names = fs.directories_in_path(method_path, exact_not_name=directory_names_not_methods, returns="name")

                # No extra maps
                if len(subdirectory_names) == 0: files = dict()

                # One directory with extra maps
                elif len(subdirectory_names) == 1:

                    # Get the path to the extra maps directory
                    extra_maps_path = fs.join(method_path, subdirectory_names[0])

                    # Get the FITS file paths
                    files = fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

                # More?
                else:

                    # GIVE WARNING
                    log.warning("Please remove directories other than 'plots', 'contours', and 'profiles'")

                    # Get the name and path
                    extra_maps_name = prompt_string("extra_directory_name", "directory with the extra maps in '" + fs.name(sub_path) + "/" + method_name, choices=subdirectory_names, required=False)
                    if extra_maps_name is None: files = dict()  # none of the directories
                    extra_maps_path = fs.join(method_path, extra_maps_name)

                    # Get the FITS file paths
                    files = fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

                # Flatten into a one-level dict
                if flatten:
                    for map_name in files: paths[method_name + "_" + map_name] = files[map_name]

                # Don't flatten
                else: paths[method_name] = files

            # Return the paths
            return paths

    # Files present
    elif fs.contains_files(sub_path):

        # Method cannot be defined
        if method is not None: raise ValueError("Specified method '" + method + "', but all maps are in one directory")

        # Methods cannot be defined
        if methods is not None: raise ValueError("Specified methods '" + str(methods) + "', but all maps are in one directory")

        # Not method cannot be defined
        if not_methods is not None: raise ValueError("All maps are in one directory (no different methods)")

        # Find a directory that contains extra maps (a directorty that isn't called 'plot', 'contours', or 'profiles'
        subdirectory_names = fs.directories_in_path(sub_path, exact_not_name=directory_names_not_methods, returns="name")

        # No extra maps
        if len(subdirectory_names) == 0: return dict()

        # One directory with extra maps
        elif len(subdirectory_names) == 1:

            # Get the path to the extra maps directory
            extra_maps_path = fs.join(sub_path, subdirectory_names[0])

            # Return the FITS file paths
            return fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

        # More?
        else:

            # GIVE WARNING
            log.warning("Please remove directories other than 'plots', 'contours', and 'profiles'")

            # Get the name and path
            extra_maps_name = prompt_string("extra_directory_name", "directory with the extra maps in '" + fs.name(sub_path), choices=subdirectory_names, required=False)
            if extra_maps_name is None: return dict()  # none of the directories
            extra_maps_path = fs.join(sub_path, extra_maps_name)

            # Return the FITS file paths
            return fs.files_in_path(extra_maps_path, returns="dict", extension="fits")

    # Nothing present
    else: return dict()

# -----------------------------------------------------------------

def get_maps_sub_name_analysis(analysis_run, name, flatten=False, framelist=False, method=None, methods=None,
                               not_method=None, not_methods=None, factors=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :param factors:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name_analysis(analysis_run, name, flatten=flatten, method=method, methods=methods,
                                            not_method=not_method, not_methods=not_methods, factors=factors)

    # Return the maps
    return get_maps_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_extra_maps_sub_name_analysis(analysis_run, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    Thisnf unction ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_extra_map_paths_sub_name_analysis(analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return the maps
    return get_extra_maps_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_negatives_sub_name_analysis(analysis_run, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name_analysis(analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return the map paths
    return get_negatives_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_nans_sub_name_analysis(analysis_run, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name_analysis(analysis_run, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return the map paths
    return get_nans_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_maps_sub_name(environment, history, name, flatten=False, framelist=False, method=None, methods=None,
                      not_method=None, not_methods=None, factors=None):

    """
    This function ...
    :param environment:
    :param history:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :param factors:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name(environment, name, flatten=flatten, method=method, methods=methods,
                                   not_method=not_method, not_methods=not_methods, factors=factors)

    # Return the maps
    return get_maps_sub_name_from_paths(paths, history=history, framelist=framelist)

# -----------------------------------------------------------------

def get_extra_maps_sub_name(environment, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param environment:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_extra_map_paths_sub_name(environment, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return
    return get_extra_maps_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_negatives_sub_name(environment, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This fnuction ...
    :param environment:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name(environment, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return the negative masks
    return get_negatives_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_nans_sub_name(environment, name, flatten=False, framelist=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param environment:
    :param name:
    :param flatten:
    :param framelist:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Get map paths
    paths = get_map_paths_sub_name(environment, name, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

    # Return the nans masks
    return get_nans_sub_name_from_paths(paths, framelist=framelist)

# -----------------------------------------------------------------

def get_maps_sub_name_from_paths(paths, history=None, framelist=False, images=True):

    """
    This function ...
    :param paths:
    :param history:
    :param framelist:
    :param images:
    :return:
    """

    # Initialize the maps dictionary
    maps = dict()

    # Loop over the entries
    for method_or_name in paths:

        # Methods
        if types.is_dictionary(paths[method_or_name]):

            method_name = method_or_name
            maps[method_name] = dict()

            # Loop over the paths, load the maps and add to dictionary
            for name in paths[method_or_name]:

                # Get the path
                map_path = paths[method_or_name][name]

                # Try to load
                try:
                    if images: maps[method_name][name] = Image.from_file(map_path, no_filter=True)
                    else: maps[method_name][name] = Frame.from_file(map_path, no_filter=True)
                except IOError:
                    command = command_for_sub_name(name)
                    log.warning("The " + method_name + "/" + name + " map is probably damaged. Run the '" + command + "' command again.")
                    log.warning("Removing the " + map_path + " map ...")
                    fs.remove_file(map_path)
                    history.remove_entries_and_save(command)

        # Just maps
        elif types.is_string_type(paths[method_or_name]):

            name = method_or_name

            # Get the path
            map_path = paths[method_or_name]

            # Try to load
            try:
                if images: maps[name] = Image.from_file(map_path, no_filter=True)
                else: maps[name] = Frame.from_file(map_path, no_filter=True)
            except IOError:
                command = command_for_sub_name(name)
                log.warning("The " + name + " map is probably damaged. Run the '" + command + "' command again.")
                log.warning("Removing the " + map_path + " map ...")
                fs.remove_file(map_path)
                history.remove_entries_and_save(command)

        # Something wrong
        else: raise RuntimeError("Something went wrong")

    # Return the maps
    if framelist: return NamedFrameList(**maps)
    else: return maps

# -----------------------------------------------------------------

def get_extra_maps_sub_name_from_paths(paths, framelist=False, images=True):

    """
    This function ...
    :param paths:
    :param framelist:
    :param images:
    :return:
    """

    # Initialize the maps dictionary
    maps = dict()

    # Loop over the entries
    for method_or_name in paths:

        # Methods
        if types.is_dictionary(paths[method_or_name]):

            method_name = method_or_name
            maps[method_name] = dict()

            # Loop over the paths, load the maps and add to dictionary
            for name in paths[method_or_name]:

                # Get the path
                map_path = paths[method_or_name][name]

                # Get the map
                if images: maps[method_name][name] = Image.from_file(map_path, no_filter=True)
                else: maps[method_name][name] = Frame.from_file(map_path, no_filter=True)

        # Just maps
        elif types.is_string_type(paths[method_or_name]):

            name = method_or_name

            # Get the path
            map_path = paths[method_or_name]

            # Get the map
            if images: maps[name] = Image.from_file(map_path, no_filter=True)
            else: maps[name] = Frame.from_file(map_path, no_filter=True)

        # Something wrong
        else: raise RuntimeError("Something went wrong")

    # Return the maps
    if framelist: return NamedFrameList(**maps)
    else: return maps

# -----------------------------------------------------------------

def get_negatives_sub_name_from_paths(paths, framelist=False):

    """
    Thisf unction ...
    :param paths:
    :param framelist:
    :return:
    """

    # Initialize the negatives dictionary
    negatives = dict()

    # Loop over the entries
    for method_or_name in paths:

        # Methods
        if types.is_dictionary(paths[method_or_name]):

            method_name = method_or_name
            negatives[method_name] = dict()

            # Loop over the paths
            for name in paths[method_or_name]:

                # Get the path
                map_path = paths[method_or_name][name]

                # Check mask planes
                if "negatives" in get_mask_names(map_path): negatives[method_name][name] = Mask.from_file(map_path, plane="negatives")
                else: negatives[method_name][name] = None

        # Just maps
        elif types.is_string_type(paths[method_or_name]):

            name = method_or_name

            # Get the path
            map_path = paths[method_or_name]

            # Check mask planes
            if "negatives" in get_mask_names(map_path): negatives[name] = Mask.from_file(map_path, plane="negatives")
            else: negatives[name] = None

    # Return the masks
    if framelist:
        raise NotImplementedError("Not implemented yet")
    else: return negatives

# -----------------------------------------------------------------

def get_nans_sub_name_from_paths(paths, framelist=False):

    """
    This function ...
    :param paths:
    :param framelist:
    :return:
    """

    # Initialize the NaN maps dictionary
    nans = dict()

    # Loop over the entries
    for method_or_name in paths:

        # Methods
        if types.is_dictionary(paths[method_or_name]):

            method_name = method_or_name
            nans[method_name] = dict()

            # Loop over the paths, load the maps and add to dictionary
            for name in paths[method_or_name]:

                # Get the path
                map_path = paths[method_or_name][name]

                # Check mask planes
                if "nans" in get_mask_names(map_path): nans[method_name][name] = Mask.from_file(map_path, plane="nans")
                else: nans[method_name][name] = None

        # Just maps
        elif types.is_string_type(paths[method_or_name]):

            name = method_or_name

            # Get the path
            map_path = paths[method_or_name]

            # Check mask planes
            if "nans" in get_mask_names(map_path): nans[name] = Mask.from_file(map_path, plane="nans")
            else: nans[name] = None

    # Return the masks
    if framelist:
        #return NamedFrameList(**maps)
        raise NotImplementedError("Not implemented yet")
    else: return nans

# -----------------------------------------------------------------

def get_origins_sub_name_analysis(analysis_run, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(analysis_run.maps_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Return
    return get_origins_sub_name_from_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_origins_sub_name(environment, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param environment:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(environment.maps_raw_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Return
    return get_origins_sub_name_from_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_origins_sub_name_from_path(sub_path, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param sub_path:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    direct_origins_path = fs.join(sub_path, origins_filename)

    if not_method is not None:
        if not types.is_string_type(not_method): raise ValueError("'not_method' must be a string")
        if not_methods is not None: raise ValueError("Cannot specifify both 'not_method' and 'not_methods'")
        not_methods = [not_method]

    # No subdirectories
    if fs.is_file(direct_origins_path):

        if method is not None: raise ValueError("No different methods")
        if methods is not None: raise ValueError("No different methods")
        if not_methods is not None: raise ValueError("No different methods")

        origins = load_dict(direct_origins_path)

    # Subdirectories
    else:

        if method is not None:

            if not_methods is not None: raise ValueError("Cannot specify both 'method' and 'not_methods' simultaneously")

            # Check whether valid method
            method_path = fs.join(sub_path, method)
            if not fs.is_directory(method_path): raise ValueError("Could not find a directory for the '" + method + "' method")
            origins_path = fs.join(method_path, origins_filename)
            if not fs.is_file(origins_path): raise ValueError("File '" + origins_path + "' is missing")

            # Load the origins
            origins = load_dict(origins_path)

        else:

            # Initialize
            origins = dict()

            # Loop over subdirectories
            for method_path in fs.directories_in_path(sub_path):

                origins_path = fs.join(method_path, origins_filename)
                if not fs.is_file(origins_path): raise ValueError("File '" + origins_path + "' is missing")

                # Determine method
                method_name = fs.name(method_path)

                # SKip if not in methods
                if methods is not None and method_name not in methods: continue

                # Skip methods
                if not_methods is not None and method_name in not_methods: continue

                # Load the origins for this method
                origins_method = load_dict(origins_path)

                # Flatten into a one-level dict
                if flatten:
                    for map_name in origins_method: origins[method_name + "_" + map_name] = origins_method[map_name]

                # Don't flatten: get nested dict
                else: origins[method_name] = origins_method

    # Return the origins
    return origins

# -----------------------------------------------------------------

def get_methods_sub_name_analysis(analysis_run, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param analysis_run:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(analysis_run.maps_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    # Return methods from sub path
    return get_methods_sub_name_from_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_methods_sub_name(environment, name, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param environment:
    :param name:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine path
    sub_path = fs.join(environment.maps_raw_path, name)
    if not fs.is_directory(sub_path): raise ValueError("Invalid name '" + name + "'")

    return get_methods_sub_name_from_path(sub_path, flatten=flatten, method=method, methods=methods, not_method=not_method, not_methods=not_methods)

# -----------------------------------------------------------------

def get_methods_sub_name_from_path(sub_path, flatten=False, method=None, methods=None, not_method=None, not_methods=None):

    """
    This function ...
    :param sub_path:
    :param flatten:
    :param method:
    :param methods:
    :param not_method:
    :param not_methods:
    :return:
    """

    # Determine methods path
    direct_methods_path = fs.join(sub_path, methods_filename)

    if not_method is not None:
        if not types.is_string_type(not_method): raise ValueError("'not_method' must be a string")
        if not_methods is not None: raise ValueError("Cannot specifify both 'not_method' and 'not_methods'")
        not_methods = [not_method]

    # No subdirectories
    if fs.is_file(direct_methods_path):

        if method is not None: raise ValueError("No different methods")
        if not_methods is not None: raise ValueError("No different methods")
        if methods is not None: raise ValueError("No different methods")

        the_methods = load_dict(direct_methods_path)

    # Subdirectories
    else:

        if method is not None:

            if not_methods is not None: raise ValueError("Cannot specify both 'method' and 'not_methods' simultaneously")

            # Check whether valid method
            method_path = fs.join(sub_path, method)
            if not fs.is_directory(method_path): raise ValueError("Could not find a directory for the '" + method + "' method")

            methods_filepath = fs.join(method_path, methods_filename)
            if not fs.is_file(methods_filepath): raise ValueError("File '" + methods_filepath + "' is missing")

            # Load the origins
            the_methods = load_dict(methods_filepath)

        else:

            # Initialize
            the_methods = dict()

            # Loop over subdirectories
            for method_path in fs.directories_in_path(sub_path):

                methods_filepath = fs.join(method_path, methods_filename)
                if not fs.is_file(methods_filepath): raise ValueError("File '" + methods_filepath + "' is missing")

                # Determine method
                method_name = fs.name(method_path)

                # Skip if not in methods
                if methods is not None and method_name not in methods: continue

                # Skip methods
                if not_methods is not None and method_name in not_methods: continue

                # Load the origins for this method
                methods_method = load_dict(methods_filepath)

                # Flatten into a one-level dict
                if flatten:
                    for map_name in methods_method: the_methods[method_name + "_" + map_name] = methods_method[map_name]

                # Don't flatten: get nested dict
                else: the_methods[method_name] = methods_method

    # Return the methods
    return the_methods

# -----------------------------------------------------------------
