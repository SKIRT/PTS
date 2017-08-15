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
from ...magic.core.dataset import DataSet
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

# -----------------------------------------------------------------

config_filename = "modeling.cfg"
history_filename = "history.dat"

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        from ..maps.collection import MapsCollection
        return MapsCollection.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_suite(self):

        """
        This function ...
        :return:
        """

        from ..build.suite import ModelSuite
        return ModelSuite.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_context(self):

        """
        This function ...
        :return:
        """

        from ..fitting.context import FittingContext
        return FittingContext.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):

        """
        Thi function ...
        :return:
        """

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

# -----------------------------------------------------------------

colours_name = "colours"
ssfr_name = "ssfr"
tir_name = "tir"
attenuation_name = "attenuation"
old_name = "old"
young_name = "young"
ionizing_name = "ionizing"
dust_name = "dust"

# -----------------------------------------------------------------

map_sub_names = [colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name]

# -----------------------------------------------------------------

fluxes_name = "fluxes.dat"
properties_name = "properties.dat"
info_name = "info.dat"
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
significance_maps_filename = "significance.html"

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
        self.deprojection_path = fs.create_directory_in(self.path, deprojection_name)
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
        self.data_seds_path = fs.create_directory_in(self.data_path, "SEDs")

        # Set the ...
        self.data_images_path = fs.create_directory_in(self.data_path, "images")

        # DIFFERENT MAPS SUBDIRECTORIES:

        # The maps raw directory
        self.maps_raw_path = fs.create_directory_in(self.maps_path, maps_raw_name)

        # The path to the maps/colours directory
        self.maps_colours_path = fs.create_directory_in(self.maps_raw_path, colours_name)

        # The path to the maps/ssfr directory
        self.maps_ssfr_path = fs.create_directory_in(self.maps_raw_path, ssfr_name)

        # The path to the maps/TIR directory
        self.maps_tir_path = fs.create_directory_in(self.maps_raw_path, tir_name)

        # The path to the maps/attenuation directory
        self.maps_attenuation_path = fs.create_directory_in(self.maps_raw_path, attenuation_name)

        # Set the path to the maps/old directory
        self.maps_old_path = fs.create_directory_in(self.maps_raw_path, old_name)

        # Set the path to the maps/young directory
        self.maps_young_path = fs.create_directory_in(self.maps_raw_path, young_name)

        # Set the path to the maps/ionizing directory
        self.maps_ionizing_path = fs.create_directory_in(self.maps_raw_path, ionizing_name)

        # Set the path to the maps/dust directory
        self.maps_dust_path = fs.create_directory_in(self.maps_raw_path, dust_name)

        # NEW: Set the path to the maps/html directory
        self.maps_html_path = fs.create_directory_in(self.maps_path, html_name)
        self.all_maps_html_page_path = fs.join(self.maps_html_path, all_maps_filename)
        self.maps_summary_html_page_path = fs.join(self.maps_html_path, maps_summary_filename)
        self.significance_maps_html_page_path = fs.join(self.maps_html_path, significance_maps_filename)

        # NEW

        self.html_status_path = fs.join(self.html_path, status_name)
        self.html_images_path = fs.create_directory_in(self.html_path, "images")

        # NEW

        # Set the path to the initial dataset file
        self.initial_dataset_path = fs.join(self.prep_path, initial_dataset_name)

        # Set the path to the prepared dataset file
        self.prepared_dataset_path = fs.join(self.prep_path, dataset_name)

        # Set the path to the preparation statistics file
        self.preparation_statistics_path = fs.join(self.prep_path, statistics_name)

    # -----------------------------------------------------------------

    @property
    def maps_colours_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_colours_path)

    # -----------------------------------------------------------------

    @property
    def maps_ssfr_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_ssfr_path)

    # -----------------------------------------------------------------

    @property
    def maps_tir_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_tir_path)

    # -----------------------------------------------------------------

    @property
    def maps_attenuation_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_attenuation_path)

    # -----------------------------------------------------------------

    @property
    def maps_old_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_old_path)

    # -----------------------------------------------------------------

    @property
    def maps_young_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_young_path)

    # -----------------------------------------------------------------

    @property
    def maps_ionizing_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_ionizing_path)

    # -----------------------------------------------------------------

    @property
    def maps_dust_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_dust_path)

    # -----------------------------------------------------------------

    @property
    def cache_host_id(self):

        """
        This function ...
        :return: 
        """

        return self.modeling_configuration.cache_host_id

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_host(self):

        """
        This function ...
        :return:
        """

        return load_host(self.cache_host_id)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_remote(self):

        """
        This function ...
        :return:
        """

        return Remote(host_id=self.cache_host_id)

    # -----------------------------------------------------------------

    @property
    def has_initial_dataset(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.initial_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_dataset(self):

        """
        This function ...
        :return:
        """

        return DataSet.from_file(self.initial_dataset_path, check=False)  # don't check whether the file are actually present (caching on remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_names(self):

        """
        This function ...
        :return:
        """

        return sorted(self.initial_dataset.names, key=lambda filter_name: parse_filter(filter_name).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_paths(self):

        """
        This function ...
        :return:
        """

        return sorted(self.initial_dataset.path_list, key=lambda path: parse_filter(fs.name(fs.directory_of(path))).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @property
    def nimages(self):

        """
        This function ...
        :return:
        """

        return len(self.preparation_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(name) for name in self.preparation_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return [fltr.wavelength for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.min_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.max_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def pixelscale_range(self):

        """
        This function ...
        :return:
        """

        # WARNING: INITIALIZED FILES CAN BE CACHED!
        return self.initial_dataset.pixelscale_range

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        # WARNING: INITIALIZED FILES CAN BE CACHED!
        #return self.initial_dataset.min_wavelength

        return min(fltr.wavelength for fltr in self.filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        # WARNING: INITIALIZED FILES CAN BE CACHED!
        #return self.initial_dataset.max_wavelength

        return max(fltr.wavelength for fltr in self.filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_range(self):

        """
        This function ...
        :return:
        """

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

        colours = dict()
        for name, colour in zip(self.preparation_names, self.plotting_colours):
            colours[name] = colour
        return colours

    # -----------------------------------------------------------------

    @property
    def has_dataset(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.prepared_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dataset(self):

        """
        This funtion ...
        :return:
        """

        return DataSet.from_file(self.prepared_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_framelist(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def named_frame_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_framelist(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def errormap_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_errormap_list(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def named_errormap_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_errormap_list(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_path_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_frame_path_list(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def named_frame_path_list(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_frame_path_list(named=True)  # on name

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_context(self):

        """
        This function ...
        :return:
        """

        from ..analysis.context import AnalysisContext
        return AnalysisContext.from_modeling_path(self.path)

    # -----------------------------------------------------------------

    @property
    def analysis_runs(self):

        """
        This function ...
        :return:
        """

        return self.analysis_context.runs

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

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_center(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.center

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_inclination(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.inclination

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_redshift(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.disk_ellipse.angle

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
    def physical_truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")
        # return self.truncation_area

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

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse.bounding_box

    # -----------------------------------------------------------------

    @lazyproperty
    def significance_levels(self):

        """
        This function ...
        :return:
        """

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
