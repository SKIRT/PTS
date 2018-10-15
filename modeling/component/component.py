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
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.filter.broad import BroadBandFilter
from ...core.basics.configuration import Configuration
from ..core.history import ModelingHistory
from ..core.commands import ModelingCommands
from ..core.environment import GalaxyModelingEnvironment, SEDModelingEnvironment, ImagesModelingEnvironment
from ...core.tools.utils import lazyproperty
from ...core.tools import parsing

# -----------------------------------------------------------------

class ModelingComponent(Configurable):
    
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
        super(ModelingComponent, self).__init__(*args, **kwargs)

        # The modeling configuration file
        self.config_file_path = None

        # The modeling environemnt
        self.environment = None

        # PTS directories
        self.kernels_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup(**kwargs)

        # Determine the path to the modeling configuration file
        self.config_file_path = fs.join(self.config.path, "modeling.cfg")

        # Check for the presence of the configuration file
        if not fs.is_file(self.config_file_path): raise ValueError("The current working directory (" + self.config.path + ") is not a radiative transfer modeling directory (the configuration file is missing)")

        # Determine the path to the kernels user directory
        self.kernels_path = fs.join(introspection.pts_user_dir, "kernels")

        # Create the modeling environment
        if self.is_galaxy_modeling: self.environment = GalaxyModelingEnvironment(self.config.path)
        elif self.is_sed_modeling: self.environment = SEDModelingEnvironment(self.config.path)
        elif self.is_images_modeling: self.environment = ImagesModelingEnvironment(self.config.path)

    # -----------------------------------------------------------------

    @property
    def history_file_path(self):
        return self.environment.history_file_path

    # -----------------------------------------------------------------

    @property
    def commands_file_path(self):
        return self.environment.commands_file_path

    # -----------------------------------------------------------------

    @property
    def fit_path(self):
        return self.environment.fit_path

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):
        return self.environment.analysis_path

    # -----------------------------------------------------------------

    @property
    def reports_path(self):
        return self.environment.reports_path

    # -----------------------------------------------------------------

    @property
    def visualisation_path(self):
        return self.environment.visualisation_path

    # -----------------------------------------------------------------

    @property
    def plot_path(self):
        return self.environment.plot_path

    # -----------------------------------------------------------------

    @property
    def log_path(self):
        return self.environment.log_path

    # -----------------------------------------------------------------

    @property
    def config_path(self):
        return self.environment.config_path

    # -----------------------------------------------------------------

    @property
    def show_path(self):
        return self.environment.show_path

    # -----------------------------------------------------------------

    @property
    def build_path(self):
        return self.environment.build_path

    # -----------------------------------------------------------------

    @property
    def html_path(self):
        return self.environment.html_path

    # -----------------------------------------------------------------

    @property
    def object_name(self):
        return self.modeling_configuration.name

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):
        # Return the observed SED
        if self.is_galaxy_modeling: return self.environment.observed_sed
        elif self.is_sed_modeling: return self.environment.observed_sed
        else: raise ValueError("Observed SED is not defined for modeling types other than 'galaxy' or 'sed'")

    # -----------------------------------------------------------------

    @lazyproperty
    def truncated_sed(self):
        if not self.is_galaxy_modeling: raise RuntimeError("Something went wrong")
        return self.environment.truncated_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_path(self):
        # Return the correct path
        if self.is_galaxy_modeling: return self.environment.observed_sed_path
        elif self.is_sed_modeling: return self.environment.sed_path
        else: raise ValueError("Observed SED not defined for modeling types other than 'galaxy' or 'sed'")

    # -----------------------------------------------------------------

    @lazyproperty
    def truncated_sed_path(self):
        if not self.is_galaxy_modeling: raise RuntimeError("Something went wrong")
        return self.environment.truncated_sed_path

    # -----------------------------------------------------------------

    def observed_flux(self, fltr, unit=None, add_unit=True):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.observed_sed.photometry_for_filter(fltr, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters(self):
        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names(self):
        return [str(fltr) for fltr in self.observed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_wavelengths(self):
        return [fltr.wavelength for fltr in self.observed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filters(self):
        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filter_names(self):
        return [str(fltr) for fltr in self.sed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filter_wavelengths(self):
        return [fltr.pivot for fltr in self.sed_filters]

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

    @property
    def is_galaxy_modeling(self):
        return self.modeling_configuration.modeling_type == "galaxy"

    # -----------------------------------------------------------------

    @property
    def is_sed_modeling(self):
        return self.modeling_configuration.modeling_type == "sed"

    # -----------------------------------------------------------------

    @property
    def is_images_modeling(self):
        return self.modeling_configuration.modeling_type == "images"

    # -----------------------------------------------------------------

    @lazyproperty
    def history(self):
        return self.environment.history

    # -----------------------------------------------------------------

    @lazyproperty
    def status(self):
        return self.environment.status

    # -----------------------------------------------------------------

    @lazyproperty
    def nuv_filter(self):
        return BroadBandFilter("GALEX NUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def nuv_wavelength(self):
        return self.nuv_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):
        return BroadBandFilter("GALEX FUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_wavelength(self):
        return self.fuv_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):
        return BroadBandFilter("IRAC I1")

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_wavelength(self):
        return self.i1_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def i2_filter(self):
        return BroadBandFilter("IRAC I2")

    # -----------------------------------------------------------------

    @lazyproperty
    def jhk_filters(self):
        return parsing.lazy_broad_band_filter_list("2MASS")

    # -----------------------------------------------------------------

    @lazyproperty
    def pacs_red_filter(self):
        return BroadBandFilter("Pacs 160")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_psw_filter(self):
        return BroadBandFilter("SPIRE PSW")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_pmw_filter(self):
        return BroadBandFilter("SPIRE PMW")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_plw_filter(self):
        return BroadBandFilter("SPIRE PLW")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_filters(self):
        return parsing.lazy_broad_band_filter_list("SPIRE")

    # -----------------------------------------------------------------

    @lazyproperty
    def planck_filters(self):
        return parsing.lazy_broad_band_filter_list("Planck")

    # -----------------------------------------------------------------

    @lazyproperty
    def hfi_filters(self):
        return parsing.lazy_broad_band_filter_list("HFI")

    # -----------------------------------------------------------------

    @lazyproperty
    def lfi_filters(self):
        return parsing.lazy_broad_band_filter_list("LFI")

    # -----------------------------------------------------------------

    @lazyproperty
    def iras_filters(self):
        return parsing.lazy_broad_band_filter_list("IRAS")

    # -----------------------------------------------------------------

    @lazyproperty
    def iras_and_planck_filters(self):
        return self.planck_filters + self.iras_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def ignore_sed_plot_filters(self):
        return self.planck_filters + self.iras_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_no_iras(self):
        return [fltr for fltr in self.observed_filters if fltr not in self.iras_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_no_iras(self):
        return [str(fltr) for fltr in self.observed_filters_no_iras]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_no_iras_planck(self):
        # Get the filters
        return [fltr for fltr in self.observed_filters if fltr not in self.iras_and_planck_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_no_iras_planck(self):
        return [str(fltr) for fltr in self.observed_filters_no_iras_planck]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_wavelengths_no_iras_planck(self):
        return [fltr.wavelength for fltr in self.observed_filters_no_iras_planck]

    # -----------------------------------------------------------------

    @lazyproperty
    def v_band_wavelength(self):
        return 0.55 * Unit("micron")

    # -----------------------------------------------------------------

    @property
    def maps_collection(self):
        return self.environment.maps_collection

    # -----------------------------------------------------------------

    @property
    def static_maps_collection(self):
        return self.environment.static_maps_collection

    # -----------------------------------------------------------------

    @property
    def maps_selection(self):
        return self.environment.maps_selection

    # -----------------------------------------------------------------

    @property
    def static_maps_selection(self):
        return self.environment.static_maps_selection

    # -----------------------------------------------------------------

    @property
    def model_suite(self):
        return self.environment.model_suite

    # -----------------------------------------------------------------

    @property
    def static_model_suite(self):
        return self.environment.static_model_suite

    # -----------------------------------------------------------------

    @property
    def fitting_context(self):
        return self.environment.fitting_context

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):
        return self.environment.fitting_runs

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

def get_modeling_type(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    configuration = load_modeling_configuration(modeling_path)
    return configuration.modeling_type

# -----------------------------------------------------------------

def get_default_fitting_method(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    configuration = load_modeling_configuration(modeling_path)
    return configuration.fitting_method

# -----------------------------------------------------------------

def get_cache_host_id(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    configuration = load_modeling_configuration(modeling_path)
    return configuration.cache_host_id

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

    else:

        history = ModelingHistory.from_file(history_file_path)
        history.clean()

    # Return the history
    return history

# -----------------------------------------------------------------

def load_modeling_commands(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the commands file path
    commands_file_path = fs.join(modeling_path, "commands.txt")

    # Create new commands file
    if not fs.is_file(commands_file_path):

        commands = ModelingCommands()
        commands.saveto(commands_file_path)

    else: commands = ModelingCommands.from_file(commands_file_path)

    # Return the commands
    return commands

# -----------------------------------------------------------------

def get_configuration_file_paths(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the config path
    config_path = fs.join(modeling_path, "config")

    # Return the file paths
    return fs.files_in_path(config_path, extension="cfg")

# -----------------------------------------------------------------

def get_log_file_paths(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the log path
    log_path = fs.join(modeling_path, "log")

    # Return the file paths
    return fs.files_in_path(log_path, extension="txt")

# -----------------------------------------------------------------

def load_modeling_status(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    from ..core.status import ModelingStatus
    return ModelingStatus(modeling_path)

# -----------------------------------------------------------------