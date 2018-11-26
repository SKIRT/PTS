#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.environment Contains the ModelingEnvironment class and derived classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from abc import ABCMeta
from collections import OrderedDict

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.configuration import Configuration
from ...core.data.sed import ObservedSED
from ...core.filter.filter import parse_filter
from ...magic.core.dataset import DataSet, StaticDataSet
from ...core.basics.range import QuantityRange
from pts.core.tools.utils import lazyproperty
from ...core.remote.host import load_host
from ...core.remote.remote import Remote
from ...magic.region.list import SkyRegionList
from ...magic.region.ellipse import SkyEllipseRegion
from ...magic.basics.stretch import SkyStretch
from ...core.tools import tables
from ..basics.properties import GalaxyProperties
from ...core.tools.serialization import load_dict
from .steps import galaxy_modeling, sed_modeling, images_modeling
from ...core.filter.filter import is_uv, is_optical, is_ir, is_nir, is_mir, is_fir, is_submm, is_fir_or_submm
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...magic.core.mask import Mask

# -----------------------------------------------------------------

config_basename = "modeling"
history_basename = "history"
commands_basename = "commands"

# -----------------------------------------------------------------

config_filename = config_basename + ".cfg"
history_filename = history_basename + ".dat"
commands_filename = commands_basename + ".txt"

# -----------------------------------------------------------------

fit_name = "fit"
analysis_name = "analysis"
reports_name = "reports"
visualisation_name = "visualisation"
plot_name = "plot"
log_name = "log"
config_name = "config"
show_name = "show"
build_name = "build"
in_name = "in"
html_name = "html"

# -----------------------------------------------------------------

def load_modeling_environment(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Determine the path to the modeling configuration
    config_path = fs.join(path, config_filename)

    # Check
    if not fs.is_file(config_path): raise ValueError("The path is not a modeling directory: config file missing")

    # Load the modeling configuration
    config = Configuration.from_file(config_path)

    # Check the type, and create the appropriate environment
    if config.modeling_type == galaxy_modeling: return GalaxyModelingEnvironment(path)
    elif config.modeling_type == sed_modeling: return SEDModelingEnvironment(path)
    elif config.modeling_type == images_modeling: return ImagesModelingEnvironment(path)
    else: raise ValueError("Invalid modeling configuration file: modeling type not found")

# -----------------------------------------------------------------

def load_modeling_environment_cwd():

    """
    This function ...
    :return:
    """

    # Get the modeling path (cwd)
    modeling_path = verify_modeling_cwd()

    # Load the modeling environment
    return load_modeling_environment(modeling_path)

# -----------------------------------------------------------------

def is_modeling_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    config_file_path = fs.join(path, config_filename)

    # Check for the presence of the configuration file
    if not fs.is_file(config_file_path): return False
    else: return True

# -----------------------------------------------------------------

def verify_modeling_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    config_file_path = fs.join(path, config_filename)

    # Check for the presence of the configuration file
    if not fs.is_file(config_file_path): raise ValueError("The path '" + path + "' is not a radiative transfer modeling directory")
    else: return path

# -----------------------------------------------------------------

def verify_modeling_cwd():

    """
    This function ...
    :return:
    """

    return verify_modeling_path(fs.cwd())

# -----------------------------------------------------------------

def find_modeling_path_up(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Check input
    if not fs.is_directory(path): raise ValueError("Specified path is not a directory")

    # Find a modeling configuration file in the directory
    filepath = fs.find_file_in_path(path, exact_name=config_basename, extension="cfg", return_none=True)

    # Not succesful: find in the directory upwards of the specified directory
    if filepath is None:

        # Stop at home
        if path == fs.home: raise ValueError("Could not find a modeling configuration file and home directory is reached")
        else: return find_modeling_path_up(fs.directory_of(path))

    # Return the directory path
    else: return path

# -----------------------------------------------------------------

def find_modeling_path_up_cwd():

    """
    This function ...
    :return:
    """

    return find_modeling_path_up(fs.cwd())

# -----------------------------------------------------------------

def find_modeling_environment_up(path):

    """
    This function ...
    :param path:
    :return:
    """

    return load_modeling_environment(find_modeling_path_up(path))

# -----------------------------------------------------------------

def find_modeling_environment_up_cwd():

    """
    This funciton ...
    :return:
    """

    return load_modeling_environment(find_modeling_path_up_cwd())

# -----------------------------------------------------------------

def find_modeling_path_down(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Check input
    if not fs.is_directory(path): raise ValueError("Specified path is not a directory")

    # Find a modeling configuration file in the directory
    filepath = fs.find_file_in_path(path, exact_name=config_basename, extension="cfg", return_none=True)

    # Not succesful
    if filepath is None:

        # Remember subdirectories with a modeling configuration file
        dirs_with_modeling = []

        # Loop over the subdirectories
        for dirpath in fs.directories_in_path(path):
            if fs.find_file_in_path(dirpath, exact_name=config_basename, extension="cfg", return_none=True): dirs_with_modeling.append(dirpath)

        # Check which subdirectory(ies)
        if len(dirs_with_modeling) == 0:
            #for dirpath in fs.directories_in_path(path):
            #    return find_modeling_path_down(dirpath)
            raise NotImplementedError("Not implemented in the right way")

        # Check
        if len(dirs_with_modeling) == 1: return dirs_with_modeling[0]

        # Multiple
        else: raise ValueError("There are multiple directories at the same level with a modeling configuration file")

    # Return the directory path
    else: return path

# -----------------------------------------------------------------

def find_modeling_path_down_cwd():

    """
    This function ...
    :return:
    """

    return find_modeling_path_down(fs.cwd())

# -----------------------------------------------------------------

def find_modeling_environment_down(path):

    """
    This function ...
    :param path:
    :return:
    """

    return load_modeling_environment(find_modeling_path_down(path))

# -----------------------------------------------------------------

def find_modeling_environment_down_cwd():

    """
    This function ...
    :return:
    """

    return load_modeling_environment(find_modeling_path_down_cwd())

# -----------------------------------------------------------------

class ModelingEnvironment(object):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path
        """

        # Set the modeling base path
        self.path = modeling_path

        # Determine the path to the modeling configuration file
        self.config_file_path = fs.join(self.path, config_filename)

        # Determine the path to the modeling history file
        self.history_file_path = fs.join(self.path, history_filename)

        # Initialize the history file
        if not fs.is_file(self.history_file_path):
            from .history import ModelingHistory
            history = ModelingHistory()
            history.saveto(self.history_file_path)

        # Determine the path to the commands file
        self.commands_file_path = fs.join(self.path, commands_filename)

        # Initialize the commands file
        if not fs.is_file(self.commands_file_path):
            from .commands import ModelingCommands
            commands = ModelingCommands()
            commands.saveto(self.commands_file_path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.fit_path = fs.create_directory_in(self.path, fit_name)
        self.analysis_path = fs.create_directory_in(self.path, analysis_name)
        self.reports_path = fs.create_directory_in(self.path, reports_name)
        self.visualisation_path = fs.create_directory_in(self.path, visualisation_name)
        self.plot_path = fs.create_directory_in(self.path, plot_name)
        self.log_path = fs.create_directory_in(self.path, log_name)
        self.config_path = fs.create_directory_in(self.path, config_name)
        self.in_path = fs.create_directory_in(self.path, in_name)
        self.show_path = fs.create_directory_in(self.path, show_name)
        self.build_path = fs.create_directory_in(self.path, build_name)
        self.html_path = fs.create_directory_in(self.path, html_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def modeling_configuration(self):

        """
        This function ...
        :return:
        """

        # Load the configuration
        config = Configuration.from_file(self.config_file_path)

        # Return the configuration
        return config

    # -----------------------------------------------------------------

    @lazyproperty
    def modeling_type(self):
        return self.modeling_configuration.modeling_type

    # -----------------------------------------------------------------

    @lazyproperty
    def history(self):

        """
        This function ...
        :return:
        """

        # Import the class
        from .history import ModelingHistory

        # Open the modeling history
        history = ModelingHistory.from_file(self.history_file_path)
        history.clean()
        return history

    # -----------------------------------------------------------------

    @lazyproperty
    def status(self):

        """
        This function ...
        :return:
        """

        # Import the class
        from .status import ModelingStatus

        # Open the modeling status
        status = ModelingStatus(self.path)
        return status

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_collection(self):
        from ..maps.collection import MapsCollection
        return MapsCollection.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_maps_collection(self):
        from ..maps.collection import StaticMapsCollection
        return StaticMapsCollection.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_selection(self):
        from ..maps.selection import ComponentMapsSelection
        return ComponentMapsSelection.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_maps_selection(self):
        from ..maps.selection import StaticComponentMapsSelection
        return StaticComponentMapsSelection.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_suite(self):
        from ..build.suite import ModelSuite
        return ModelSuite.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_model_suite(self):
        from ..build.suite import StaticModelSuite
        return StaticModelSuite.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_context(self):
        from ..fitting.context import FittingContext
        return FittingContext.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):
        return self.fitting_context.runs

# -----------------------------------------------------------------

data_name = "data"
prep_name = "prep"
inspect_name = "inspect"
truncated_name = "truncated"
phot_name = "phot"
maps_name = "maps"
components_name = "components"
deprojection_name = "deprojection"
playground_name = "playground"

# -----------------------------------------------------------------

maps_raw_name = "raw"
maps_components_name = "components"

# -----------------------------------------------------------------

colours_name = "colours"
ssfr_name = "ssfr"
tir_name = "tir"
attenuation_name = "attenuation"
old_name = "old"
young_name = "young"
ionizing_name = "ionizing"
dust_name = "dust"
rgb_name = "rgb"

# -----------------------------------------------------------------

map_sub_names = [colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name, rgb_name]

# -----------------------------------------------------------------

fluxes_name = "fluxes.dat"
asymptotic_sed_name = "asymptotic.dat"
truncated_sed_name = "truncated.dat"
properties_name = "properties.dat"
info_name = "info.dat"

# -----------------------------------------------------------------

index_filename = "index.html"
status_name = "status.html"

# -----------------------------------------------------------------

initial_dataset_name = "initial_dataset.dat"
dataset_name = "dataset.dat"
statistics_name = "statistics.dat"

# -----------------------------------------------------------------

disk_region_filename = "disk.reg"
truncation_ellipse_filename = "ellipse.reg"
significance_levels_filename = "levels.dat"

# -----------------------------------------------------------------

all_maps_filename = "all.html"
maps_summary_filename = "summary.html"
clip_maps_filename = "clip.html"
old_maps_filename = "old.html"
young_maps_filename = "young.html"
ionizing_maps_filename = "ionizing.html"
dust_maps_filename = "dust.html"
maps_selection_filename = "selection.html"

# -----------------------------------------------------------------

seds_dirname = "SEDs"
images_dirname = "images"

# -----------------------------------------------------------------

class GalaxyModelingEnvironment(ModelingEnvironment):

    """
    This function ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Call the constructor of the base class
        super(GalaxyModelingEnvironment, self).__init__(modeling_path)

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.data_path = fs.create_directory_in(self.path, data_name)
        self.prep_path = fs.create_directory_in(self.path, prep_name)
        self.inspect_path = fs.create_directory_in(self.path, inspect_name)
        self.truncation_path = fs.create_directory_in(self.path, truncated_name)
        self.phot_path = fs.create_directory_in(self.path, phot_name)
        self.maps_path = fs.create_directory_in(self.path, maps_name)
        self.components_path = fs.create_directory_in(self.path, components_name)
        self.playground_path = fs.create_directory_in(self.path, playground_name)

        # DISK REGION PATH
        self.disk_region_path = fs.join(self.components_path, disk_region_filename)

        # TRUNCATION ELLIPSE PATH
        self.truncation_ellipse_path = fs.join(self.truncation_path, truncation_ellipse_filename)

        # SIGNIFICANCE LEVELS PATH
        self.significance_levels_path = fs.join(self.truncation_path, significance_levels_filename)

        # DATA

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, fluxes_name)

        # Set the path to the galaxy properties file
        self.galaxy_properties_path = fs.join(self.data_path, properties_name)

        # Set the path to the galaxy info file
        self.galaxy_info_path = fs.join(self.data_path, info_name)

        # Set the ...
        self.data_seds_path = fs.create_directory_in(self.data_path, seds_dirname)

        # Set the ...
        self.data_images_path = fs.create_directory_in(self.data_path, images_dirname)

        # PHOTOMETRY
        self.observed_sed_path = fs.join(self.phot_path, fluxes_name)
        self.asymptotic_sed_path = fs.join(self.phot_path, asymptotic_sed_name)
        self.truncated_sed_path = fs.join(self.phot_path, truncated_sed_name)
        self.phot_images_path = fs.create_directory_in(self.phot_path, "images")

        # Truncation HTML path
        self.truncation_html_path = fs.create_directory_in(self.truncation_path, html_name)

        # DIFFERENT MAPS SUBDIRECTORIES:

        # The maps raw directory
        self.maps_raw_path = fs.create_directory_in(self.maps_path, maps_raw_name)

        # Set the paths to the raw maps subdirectories
        self.maps_colours_path = fs.create_directory_in(self.maps_raw_path, colours_name)
        self.maps_ssfr_path = fs.create_directory_in(self.maps_raw_path, ssfr_name)
        self.maps_tir_path = fs.create_directory_in(self.maps_raw_path, tir_name)
        self.maps_attenuation_path = fs.create_directory_in(self.maps_raw_path, attenuation_name)
        self.maps_old_path = fs.create_directory_in(self.maps_raw_path, old_name)
        self.maps_young_path = fs.create_directory_in(self.maps_raw_path, young_name)
        self.maps_ionizing_path = fs.create_directory_in(self.maps_raw_path, ionizing_name)
        self.maps_dust_path = fs.create_directory_in(self.maps_raw_path, dust_name)

        # NEW: Set the path to the maps/html directory
        self.maps_html_path = fs.create_directory_in(self.maps_path, html_name)

        # Set the paths to other pages
        self.all_maps_html_page_path = fs.join(self.maps_html_path, all_maps_filename)
        self.maps_summary_html_page_path = fs.join(self.maps_html_path, maps_summary_filename)
        self.clip_maps_html_page_path = fs.join(self.maps_html_path, clip_maps_filename)
        self.old_maps_html_page_path = fs.join(self.maps_html_path, old_maps_filename)
        self.young_maps_html_page_path = fs.join(self.maps_html_path, young_maps_filename)
        self.ionizing_maps_html_page_path = fs.join(self.maps_html_path, ionizing_maps_filename)
        self.dust_maps_html_page_path = fs.join(self.maps_html_path, dust_maps_filename)
        self.maps_selection_html_page_path = fs.join(self.maps_html_path, maps_selection_filename)

        # NEW: Set the path to the maps/components directory
        self.maps_components_path = fs.create_directory_in(self.maps_path, maps_components_name)

        # Set the paths to the component maps subdirectories
        self.maps_old_component_path = fs.create_directory_in(self.maps_components_path, old_name)
        self.maps_young_component_path = fs.create_directory_in(self.maps_components_path, young_name)
        self.maps_ionizing_component_path = fs.create_directory_in(self.maps_components_path, ionizing_name)
        self.maps_dust_component_path = fs.create_directory_in(self.maps_components_path, dust_name)

        # NEW

        self.html_index_path = fs.join(self.html_path, index_filename)
        self.html_status_path = fs.join(self.html_path, status_name)

        # Images directory
        self.html_images_path = fs.create_directory_in(self.html_path, "images")

        # Scripts directory
        self.html_scripts_path = fs.create_directory_in(self.html_path, "scripts")

        # NEW

        # Set the path to the initial dataset file
        self.initial_dataset_path = fs.join(self.prep_path, initial_dataset_name)

        # Set the path to the prepared dataset file
        self.prepared_dataset_path = fs.join(self.prep_path, dataset_name)

        # Set the path to the preparation statistics file
        self.preparation_statistics_path = fs.join(self.prep_path, statistics_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):
        return ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def truncated_sed(self):
        return ObservedSED.from_file(self.truncated_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def asymptotic_sed(self):
        return ObservedSED.from_file(self.asymptotic_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dustpedia_sed(self):
        return ObservedSED.from_file(self.observed_sed_dustpedia_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_dataset(self):
        return DataSet.from_directory(self.phot_images_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def static_photometry_dataset(self):
        return StaticDataSet.from_directory(self.phot_images_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_names(self):
        return fs.files_in_path(self.phot_images_path, extension="fits", returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_filters(self):
        return [parse_filter(filter_name) for filter_name in self.photometry_image_names]

    # -----------------------------------------------------------------

    def has_photometry_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        from ...core.tools import types
        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.photometry_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_paths(self):
        return fs.files_in_path(self.phot_images_path, extension="fits", returns="path")

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_paths_for_filters(self):

        """
        This function ...
        :return:
        """

        paths = dict()

        # Loop over the image file names
        for name in self.photometry_image_names:

            # Get the filter
            fltr = parse_filter(name)

            # Get the image path
            path = fs.join(self.phot_images_path, name + ".fits")

            # Add to dictionary
            paths[fltr] = path

        # Return the dictionary
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_paths_for_filter_names(self):
        return {str(fltr): path for fltr, path in self.photometry_image_paths_for_filters.items()}

    # -----------------------------------------------------------------

    def get_photometry_image_path(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return fs.join(self.phot_images_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_photometry_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_photometry_image_path(name))

    # -----------------------------------------------------------------

    def get_photometry_image_path_for_filter(self, fltr):

        """
        Thisfunction ...
        :param fltr:
        :return:
        """

        return self.get_photometry_image_path(str(fltr))

    # -----------------------------------------------------------------

    def has_photometry_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_photometry_image_path_for_filter(fltr))

    # -----------------------------------------------------------------

    def get_photometry_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_photometry_image_path(name)
        return Image.from_file(path)

    # -----------------------------------------------------------------

    def get_photometry_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.get_photometry_image(str(fltr))

    # -----------------------------------------------------------------

    def get_photometry_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_photometry_image_path(name)
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    def get_photometry_frame_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.get_photometry_frame(str(fltr))

    # -----------------------------------------------------------------

    def get_photometry_background(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_photometry_image_path(name)
        return Frame.from_file(path, plane="background")

    # -----------------------------------------------------------------

    def get_photometry_background_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.get_photometry_background(str(fltr))

    # -----------------------------------------------------------------

    def get_photometry_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_photometry_image_path(name)
        return Mask.nans_from_file(path)

    # -----------------------------------------------------------------

    def get_photometry_mask_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.get_photometry_mask(str(fltr))

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):
        return fs.name(self.maps_colours_path)

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):
        return fs.name(self.maps_ssfr_path)

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):
        return fs.name(self.maps_tir_path)

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):
        return fs.name(self.maps_attenuation_path)

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):
        return fs.name(self.maps_old_path)

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):
        return fs.name(self.maps_young_path)

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):
        return fs.name(self.maps_ionizing_path)

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):
        return fs.name(self.maps_dust_path)

    # -----------------------------------------------------------------

    @property
    def cache_host_id(self):
        return self.modeling_configuration.cache_host_id

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_host(self):
        return load_host(self.cache_host_id)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_remote(self):
        return Remote(host_id=self.cache_host_id)

    # -----------------------------------------------------------------

    @property
    def has_initial_dataset(self):
        return fs.is_file(self.initial_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_dataset(self):
        return DataSet.from_file(self.initial_dataset_path, check=False)  # don't check whether the file are actually present (caching on remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_names(self):
        return sorted(self.initial_dataset.names, key=lambda filter_name: parse_filter(filter_name).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_paths(self):
        return sorted(self.initial_dataset.path_list, key=lambda path: parse_filter(fs.name(fs.directory_of(path))).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @property
    def nimages(self):
        return len(self.preparation_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):
        return [parse_filter(name) for name in self.preparation_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def uv_filters(self):
        return [fltr for fltr in self.filters if is_uv(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def optical_filters(self):
        return [fltr for fltr in self.filters if is_optical(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def ir_filters(self):
        return [fltr for fltr in self.filters if is_ir(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def nir_filters(self):
        return [fltr for fltr in self.filters if is_nir(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def mir_filters(self):
        return [fltr for fltr in self.filters if is_mir(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def fir_filters(self):
        return [fltr for fltr in self.filters if is_fir(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def submm_filters(self):
        return [fltr for fltr in self.filters if is_submm(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def fir_submm_filters(self):
        return [fltr for fltr in self.filters if is_fir_or_submm(fltr)]

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):
        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def uv_filter_names(self):
        return [str(fltr) for fltr in self.uv_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def optical_filter_names(self):
        return [str(fltr) for fltr in self.optical_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def ir_filter_names(self):
        return [str(fltr) for fltr in self.ir_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def nir_filter_names(self):
        return [str(fltr) for fltr in self.nir_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def mir_filter_names(self):
        return [str(fltr) for fltr in self.mir_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def fir_filter_names(self):
        return [str(fltr) for fltr in self.fir_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def submm_filter_names(self):
        return [str(fltr) for fltr in self.submm_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):
        return [fltr.wavelength for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def min_pixelscale(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.min_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def max_pixelscale(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.max_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def pixelscale_range(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.pixelscale_range

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        #return self.initial_dataset.min_wavelength
        return min(fltr.wavelength for fltr in self.filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        #return self.initial_dataset.max_wavelength
        return max(fltr.wavelength for fltr in self.filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_range(self):
        # WARNING: INITIALIZED FILES CAN BE CACHED!
        #return self.initial_dataset.wavelength_range
        return QuantityRange(self.min_wavelength, self.max_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def plotting_colours(self):

        """
        This function ...
        :return:
        """

        # NOTE: CAN ONLY BE USED AFTER PREPARATION

        # Get the colour map
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("jet")

        min_log_wavelength = np.log10(self.min_wavelength.to("micron").value)
        max_log_wavelength = np.log10(self.max_wavelength.to("micron").value)
        log_wavelength_range = max_log_wavelength - min_log_wavelength

        normalized_wavelengths = [(np.log10(wavelength.to("micron").value) - min_log_wavelength)/log_wavelength_range for wavelength in self.wavelengths]

        # Get colours in RGBA
        colours = cmap(normalized_wavelengths)

        # Convert to HEX
        from ...core.basics.colour import rgb_to_hex
        colours = [rgb_to_hex(np.array(colour[0:3])*255) for colour in colours]

        # Return the colours
        return colours

    # -----------------------------------------------------------------

    @lazyproperty
    def plotting_colours_dict(self):

        """
        This function ...
        :return:
        """

        # NOTE: CAN ONLY BE USED AFTER PREPARATION

        colours = dict()
        for name, colour in zip(self.preparation_names, self.plotting_colours):
             colours[name] = colour
        return colours

    # -----------------------------------------------------------------

    @property
    def has_dataset(self):
        return fs.is_file(self.prepared_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dataset(self):
        return DataSet.from_file(self.prepared_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_list(self):
        return self.dataset.get_framelist(named=False)  # on filter

    # -----------------------------------------------------------------

    def get_names_for_filters(self, filters):

        """
        Thisf unction ...
        :param filters:
        :return:
        """

        return self.dataset.get_names_for_filters(filters)

    # -----------------------------------------------------------------

    def get_frames_for_filters(self, filters, framelist=False, named=True):

        """
        Thisf unction ...
        :param filters:
        :param framelist:
        :param named:
        :return:
        """

        if framelist: return self.dataset.get_framelist_for_filters(filters, named=named)
        else: return self.dataset.get_frames_for_filters(filters)

    # -----------------------------------------------------------------

    def get_frame_path(self, name):
        return self.dataset.get_frame_path(name)

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr):
        return self.dataset.get_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_frame_path_for_filter(self, fltr):
        return self.dataset.get_frame_path_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_frame(self, name):
        return self.dataset.get_frame(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def named_frame_list(self):
        return self.dataset.get_framelist(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def errormap_list(self):
        return self.dataset.get_errormap_list(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def named_errormap_list(self):
        return self.dataset.get_errormap_list(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_path_list(self):
        return self.dataset.get_frame_path_list(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def named_frame_path_list(self):
        return self.dataset.get_frame_path_list(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_context(self):
        from ..analysis.context import AnalysisContext
        return AnalysisContext.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @property
    def analysis_runs(self):
        return self.analysis_context.runs

    # -----------------------------------------------------------------

    @property
    def cached_analysis_runs(self):
        return self.analysis_context.cached_runs

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.galaxy_properties_path): raise IOError("The galaxy properties file is not present. Perform 'fetch_properties' to create this file'")

        # Load the properties
        properties = GalaxyProperties.from_file(self.galaxy_properties_path)

        # Return the property map
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_info(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.galaxy_info_path): raise IOError("The galaxy info file is not (yet) present.")

        # Load the info table
        table = tables.from_file(self.galaxy_info_path)

        # To ordered dict
        info = OrderedDict()
        for name in table.colnames: info[name] = table[name][0]

        # Return the info
        return info

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_type(self):
        return self.galaxy_info["Hubble Type"]

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_stage(self):
        return self.galaxy_info["Hubble Stage"]

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_ellipse(self):  # from properties

        """
        This function ...
        :return:
        """

        # Get properties
        center = self.galaxy_properties.center
        major = self.galaxy_properties.major_arcsec
        position_angle = self.galaxy_properties.position_angle
        ellipticity = self.galaxy_properties.ellipticity

        # 1 / axial_ratio = 1 - ellipticity
        axial_ratio = 1. / (1. - ellipticity)

        # Set radius
        minor = major * axial_ratio
        radius = SkyStretch(major, minor)

        # Create and return the region
        region = SkyEllipseRegion(center, radius, position_angle)
        return region

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_distance(self):
        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_center(self):
        return self.galaxy_properties.center

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_inclination(self):
         return self.galaxy_properties.inclination

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_position_angle(self):
        return self.galaxy_properties.position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_redshift(self):
        return self.galaxy_properties.redshift

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Open the region
        region = SkyRegionList.from_file(self.disk_region_path)

        # Return the first and only shape
        return region[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_position_angle(self):
        return self.disk_ellipse.angle

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_axial_ratio(self):
        return self.disk_ellipse.axial_ratio

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_ellipticity(self):
        return self.disk_ellipse.ellipticity

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_inclination(self):

        """
        This function ...
        :return:
        """

        #from ..decomposition.decomposition import axial_ratio_to_inclination_mosenkov
        #return axial_ratio_to_inclination_mosenkov(self.disk_axial_ratio, self.hubble_stage)
        from ..decomposition.decomposition import axial_ratio_to_inclination_with_intrinsic, q_degeyter
        return axial_ratio_to_inclination_with_intrinsic(self.disk_axial_ratio, q_degeyter)

    # -----------------------------------------------------------------

    @property
    def has_truncation_ellipse(self):
        return fs.is_file(self.truncation_ellipse_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        # Load the ellipse
        region = SkyRegionList.from_file(self.truncation_ellipse_path)
        ellipse = region[0]

        # Return the (sky) ellipse
        return ellipse

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_factor(self):

        """
        Thisf unction ...
        :return:
        """

        major_factor = self.truncation_ellipse.major / self.disk_ellipse.major
        minor_factor = self.truncation_ellipse.minor / self.disk_ellipse.minor

        # Check
        if not np.isclose(major_factor, minor_factor): raise ValueError("Inconsistent major and minor axes truncation factor")

        # Return the factor
        return major_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def physical_truncation_ellipse(self):
        return self.truncation_ellipse.to_physical(distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    def physical_truncation_ellipse_for_frame(self, frame):
        return self.physical_truncation_ellipse_for_wcs(frame.wcs)

    # -----------------------------------------------------------------

    def physical_truncation_ellipse_for_wcs(self, wcs):
        return self.truncation_ellipse.to_physical(distance=self.galaxy_distance, wcs=wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_radius(self):
        return self.truncation_ellipse.semimajor

    # -----------------------------------------------------------------

    @lazyproperty
    def physical_truncation_radius(self):
        return (self.truncation_radius * self.galaxy_distance).to("kpc", equivalencies=dimensionless_angles())

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_area(self):

        """
        This function ...
        :return:
        """

        # Convert the semi minor and semi major axis lengths from angular to physical sizes
        semimajor = (self.truncation_ellipse.semimajor * self.galaxy_distance).to("kpc", equivalencies=dimensionless_angles())
        semiminor = (self.truncation_ellipse.semiminor * self.galaxy_distance).to("kpc", equivalencies=dimensionless_angles())

        # Calculate the area in kpc^2
        # A = pi * a * b
        area = math.pi * semimajor * semiminor

        # Return the area
        return area

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_box(self):
        return self.truncation_ellipse.bounding_box

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_box_axial_ratio(self):
        return self.truncation_box.axial_ratio

    # -----------------------------------------------------------------

    @property
    def has_significance_levels(self):
        return fs.is_file(self.significance_levels_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def significance_levels(self):
        return load_dict(self.significance_levels_path)

# -----------------------------------------------------------------

input_name = "input"
sed_filename = "sed.dat"
template_ski_filename = "template.ski"

# -----------------------------------------------------------------

class SEDModelingEnvironment(ModelingEnvironment):

    """
    This class ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Call the constructor of the base class
        super(SEDModelingEnvironment, self).__init__(modeling_path)

        # Set the SED path
        self.sed_path = fs.join(self.path, sed_filename)

        # Set the SED plot path
        self.sed_plot_path = fs.join(self.path, "sed.pdf")

        # Set the ski template path
        self.ski_path = fs.join(self.path, template_ski_filename)

        # Set the ski input path
        self.ski_input_path = fs.create_directory_in(self.path, input_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return ObservedSED.from_file(self.sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_filter_names(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.filter_names()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_filters(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.filters()

# -----------------------------------------------------------------

images_name = "images"
header_filename = "header.txt"

# -----------------------------------------------------------------

class ImagesModelingEnvironment(ModelingEnvironment):

    """
    This class ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Call the constructor of the base class
        super(ImagesModelingEnvironment, self).__init__(modeling_path)

        # Set images path
        self.images_path = fs.create_directory_in(self.path, images_name)

        # Set images header path
        self.images_header_path = fs.join(self.images_path, header_filename)

        # Set the ski template path
        self.ski_path = fs.join(self.path, template_ski_filename)

        # Set the ski input path
        self.ski_input_path = fs.create_directory_in(self.path, input_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def image_filter_names(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.images_path, extension="fits", returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def image_filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(name) for name in self.image_filter_names]

# -----------------------------------------------------------------
