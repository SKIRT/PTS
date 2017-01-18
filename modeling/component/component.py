#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.component.component Contains the ModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.basics.filter import Filter
from ...core.basics.configuration import Configuration
from ..core.history import ModelingHistory

# -----------------------------------------------------------------

class ModelingComponent(Configurable):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelingComponent, self).__init__(config)

        # Attributes
        #self.galaxy_name = None

        # The modeling configuration file
        self.config_file_path = None

        # The modeling history file
        self.history_file_path = None

        # Modeling directories
        #self.data_path = None
        #self.prep_path = None
        #self.truncation_path = None
        #self.phot_path = None
        #self.maps_path = None
        #self.components_path = None
        #self.deprojection_path = None
        self.fit_path = None
        self.analysis_path = None
        self.reports_path = None
        self.visualisation_path = None
        self.plot_path = None
        self.log_path = None
        self.config_path = None
        self.show_path = None

        # PTS directories
        self.kernels_path = None

        # The path to the observed SEDs
        #self.observed_sed_path = None
        #self.observed_sed_dustpedia_path = None

        # The path to the fitting configuration file
        self.fitting_configuration_path = None

        # The path to the maps
        #self.old_stellar_map_path = None
        #self.young_stellar_map_path = None
        #self.ionizing_stellar_map_path = None
        #self.dust_map_path = None

        # The paths to the significance maps
        #self.old_stellar_significance_path = None
        #self.young_stellar_significance_path = None
        #self.ionizing_stellar_significance_path = None
        #self.dust_significance_path = None

        # The paths to the cutoff masks
        #self.old_stellar_cutoff_path = None
        #self.young_stellar_cutoff_path = None
        #self.ionizing_stellar_cutoff_path = None
        #self.dust_cutoff_path = None

        # The path to the galaxy properties file
        #self.galaxy_properties_path = None

        # The path to the components/models directory
        #self.components_models_path = None

        # The path to the components/projections directory
        #self.components_projections_path = None

        # The paths to the bulge and disk models
        #self.bulge_model_path = None
        #self.disk_model_path = None

        # The paths to the different galaxy projections
        #self.earth_projection_path = None
        #self.edgeon_projection_path = None
        #self.faceon_projection_path = None

        # The path to the components/images directory
        #self.components_images_path = None

        # The paths to the final bulge, disk and model images
        #self.bulge_image_path = None
        #self.disk_image_path = None
        #self.model_image_path = None

        # The path to the truncation/masks directory
        #self.truncation_masks_path = None

        # The path to the truncation mask of the reference image (and rebinned images in the dataset)
        #self.reference_mask_path = None

        # The path to the data/seds directory
        #self.data_seds_path = None

        # The path to the DustPedia observed SED for the galaxy
        #self.dustpedia_sed_path = None

        # The path to the data/images directory
        #self.data_images_path = None

        # The path to the disk region file
        #self.disk_region_path = None

        # The path to the initial data set file
        #self.initial_dataset_path = None

        # The path to the prepared data set file
        #self.prepared_dataset_path = None

        # The path to the preparation statistics file
        #self.preparation_statistics_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Determine the path to the modeling configuration file
        self.config_file_path = fs.join(self.config.path, "modeling.cfg")

        # Check for the presence of the configuration file
        if not fs.is_file(self.config_file_path): raise ValueError("The current working directory is not a radiative transfer modeling directory (the configuration file is missing)")

        # Determine the path to the modeling history file
        self.history_file_path = fs.join(self.config.path, "history.dat")

        # Initialize the history file
        if not fs.is_file(self.history_file_path):
            history = ModelingHistory()
            history.saveto(self.history_file_path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        #self.data_path = fs.create_directory_in(self.config.path, "data")
        #self.prep_path = fs.create_directory_in(self.config.path, "prep")
        #self.truncation_path = fs.create_directory_in(self.config.path, "truncated")
        #self.phot_path = fs.create_directory_in(self.config.path, "phot")
        #self.maps_path = fs.create_directory_in(self.config.path, "maps")
        #self.components_path = fs.create_directory_in(self.config.path, "components")
        #self.deprojection_path = fs.create_directory_in(self.config.path, "deprojection")
        self.fit_path = fs.create_directory_in(self.config.path, "fit")
        self.analysis_path = fs.create_directory_in(self.config.path, "analysis")
        self.reports_path = fs.create_directory_in(self.config.path, "reports")
        self.visualisation_path = fs.create_directory_in(self.config.path, "visualisation")
        self.plot_path = fs.create_directory_in(self.config.path, "plot")
        self.log_path = fs.create_directory_in(self.config.path, "log")
        self.config_path = fs.create_directory_in(self.config.path, "config")
        self.show_path = fs.create_directory_in(self.config.path, "show")

        # Determine the path to the kernels user directory
        self.kernels_path = fs.join(introspection.pts_user_dir, "kernels")

        # Set the path to the observed SED
        #self.observed_sed_path = fs.join(self.phot_path, "fluxes.dat")

        # Set the path to the DustPedia observed SED
        #self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

        # Set the path to the fitting configuration file
        #self.fitting_configuration_path = fs.join(self.fit_path, "configuration.cfg")

        # Set the paths to the input maps
        #self.old_stellar_map_path = fs.join(self.maps_path, "old_stars.fits")
        #self.young_stellar_map_path = fs.join(self.maps_path, "young_stars.fits")
        #self.ionizing_stellar_map_path = fs.join(self.maps_path, "ionizing_stars.fits")
        #self.dust_map_path = fs.join(self.maps_path, "dust.fits")

        # The paths to the significance masks
        #self.old_stellar_significance_path = fs.join(self.maps_path, "old_stars_significance.fits")
        #self.young_stellar_significance_path = fs.join(self.maps_path, "young_stars_significance.fits")
        #self.ionizing_stellar_significance_path = fs.join(self.maps_path, "ionizing_stars_significance.fits")
        #self.dust_significance_path = fs.join(self.maps_path, "dust_significance.fits")

        # The paths to the significance masks
        #self.old_stellar_cutoff_path = fs.join(self.maps_path, "old_stars_cutoff.fits")
        #self.young_stellar_cutoff_path = fs.join(self.maps_path, "young_stars_cutoff.fits")
        #self.ionizing_stellar_cutoff_path = fs.join(self.maps_path, "ionizing_stars_cutoff.fits")
        #self.dust_cutoff_path = fs.join(self.maps_path, "dust_cutoff.fits")

        # Set the path to the galaxy properties file
        #self.galaxy_properties_path = fs.join(self.data_path, "properties.dat")

        # Set the path to the components/models directory
        #self.components_models_path = fs.create_directory_in(self.components_path, "models")

        # Set the path to the components/projections directory
        #self.components_projections_path = fs.create_directory_in(self.components_path, "projections")

        # Set the paths to the bulge and disk models
        #self.bulge_model_path = fs.join(self.components_models_path, "bulge.mod")
        #self.disk_model_path = fs.join(self.components_models_path, "disk.mod")

        # The paths to the different galaxy projections
        #self.earth_projection_path = fs.join(self.components_projections_path, "earth.proj")
        #self.edgeon_projection_path = fs.join(self.components_projections_path, "edgeon.proj")
        #self.faceon_projection_path = fs.join(self.components_projections_path, "faceon.proj")

        # Set the path to the components/images directory
        #self.components_images_path = fs.create_directory_in(self.components_path, "images")

        # Set the path to the final bulge, disk and model images
        #self.bulge_image_path = fs.join(self.components_images_path, "bulge.fits")
        #self.disk_image_path = fs.join(self.components_images_path, "disk.fits")
        #self.model_image_path = fs.join(self.components_images_path, "model.fits")

        # Set the path to the truncation/masks directory
        #self.truncation_masks_path = fs.create_directory_in(self.truncation_path, "masks")

        # The path to the truncation mask of the reference image (and rebinned images in the dataset)
        #self.reference_mask_path = fs.join(self.truncation_masks_path, "reference.fits")

        # Set ...
        #self.data_seds_path = fs.create_directory_in(self.data_path, "SEDs")

        # The DustPedia SED path
        #self.dustpedia_sed_path = fs.join(self.data_seds_path, "DustPedia.dat")

        # Set ...
        #self.data_images_path = fs.create_directory_in(self.data_path, "images")

        # Set the path to the disk region file
        #self.disk_region_path = fs.join(self.components_path, "disk.reg")

        # Set the path to the initial dataset file
        #self.initial_dataset_path = fs.join(self.prep_path, "initial_dataset.dat")

        # Set the path to the prepared dataset file
        #self.prepared_dataset_path = fs.join(self.prep_path, "dataset.dat")

        # Set the path to the preparation statistics file
        #self.preparation_statistics_path = fs.join(self.prep_path, "statistics.dat")

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
    def history(self):

        """
        This function ...
        :return:
        """

        # Open the modeling history
        return ModelingHistory.from_file(self.history_file_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return map(Filter.from_string, self.fitting_filter_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_configuration(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.fitting_configuration_path) if fs.is_file(self.fitting_configuration_path) else None

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.filters if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_labels(self):

        """
        THIS FUNCTION GUARANTEES THAT THE LABELS ARE ALWAYS ORDERED ALPHABETICALLY !!
        :return:
        """

        return sorted(self.fitting_configuration.free_parameters) if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        ranges = dict()
        for label in self.free_parameter_labels:
            parameter_range = self.fitting_configuration[label + "_range"]
            ranges[label] = parameter_range
        return ranges

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_descriptions(self):

        """
        This function ...
        :return:
        """

        from ..config.parameters import choices
        return choices

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        from ..config.parameters import units
        return units

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):

        """
        This function ...
        :return:
        """

        return Filter.from_string("GALEX FUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):

        """
        This function ...
        :return:
        """

        return Filter.from_string("IRAC I1")

    # -----------------------------------------------------------------

    @lazyproperty
    def pacs_red_filter(self):

        """
        This function ...
        :return:
        """

        return Filter.from_string("Pacs 160")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_psw_filter(self):

        """
        This function ...
        :return:
        """

        return Filter.from_string("SPIRE PSW")

    # -----------------------------------------------------------------

    @lazyproperty
    def v_band_wavelength(self):

        """
        This function ...
        :return:
        """

        return 0.55 * Unit("micron")

# -----------------------------------------------------------------

def load_fitting_configuration(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the path to the fitting configuration file
    fitting_configuration_path = fs.join(modeling_path, "fit", "configuration.cfg")

    # Open the configuration and return it
    return Configuration.from_file(fitting_configuration_path)

# -----------------------------------------------------------------

def get_config_file_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the path to the configuration file
    path = fs.join(modeling_path, "modeling.cfg")

    # Return the path
    return path

# -----------------------------------------------------------------

def load_modeling_configuration(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the path to the modeling configuration file
    path = get_config_file_path(modeling_path)

    # Open the configuration and return it
    return Configuration.from_file(path)

# -----------------------------------------------------------------

def load_modeling_history(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine history file path
    history_file_path = fs.join(modeling_path, "history.dat")

    # Create new history
    if not fs.is_file(history_file_path):

        history = ModelingHistory()
        history.saveto(history_file_path)

    else: history = ModelingHistory.from_file(history_file_path)

    # Return the history
    return history

# -----------------------------------------------------------------
