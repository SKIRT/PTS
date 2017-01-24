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

        # The modeling configuration file
        self.config_file_path = None

        # The modeling history file
        self.history_file_path = None

        # Modeling directories
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

        # The path to the fitting configuration file
        self.fitting_configuration_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

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

        # Set the path to the fitting configuration file
        self.fitting_configuration_path = fs.join(self.fit_path, "configuration.cfg")

    # -----------------------------------------------------------------

    @property
    def object_name(self):

        """
        This function ...
        :return:
        """

        return self.modeling_configuration.name

    # -----------------------------------------------------------------

    @property
    def input_map_paths(self):

        """
        This function ...
        :return:
        """

        return []

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        from .sed import get_observed_sed as get_sed_other
        from .galaxy import get_observed_sed as get_sed_galaxy

        if self.modeling_type == "galaxy": return get_sed_galaxy(self.config.path)
        else: return get_sed_other(self.config.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_path(self):

        """
        This function ...
        :return:
        """

        from .sed import get_observed_sed_file_path as get_path_other
        from .galaxy import get_observed_sed_file_path as get_path_galaxy

        if self.modeling_type == "galaxy": return get_path_galaxy(self.config.path)
        else: return get_path_other(self.config.path)

    # -----------------------------------------------------------------

    def observed_flux(self, fltr, unit=None, add_unit=True):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.observed_sed.flux_for_filter(fltr, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filters(self):

        """
        This property ...
        :return:
        """

        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.sed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filter_wavelengths(self):

        """
        This function ...
        :return:
        """

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

        # Open the modeling history
        return ModelingHistory.from_file(self.history_file_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return map(Filter, self.fitting_filter_names)

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

        return self.fitting_configuration.descriptions

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.units

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):

        """
        This function ...
        :return:
        """

        return Filter("GALEX FUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):

        """
        This function ...
        :return:
        """

        return Filter("IRAC I1")

    # -----------------------------------------------------------------

    @lazyproperty
    def pacs_red_filter(self):

        """
        This function ...
        :return:
        """

        return Filter("Pacs 160")

    # -----------------------------------------------------------------

    @lazyproperty
    def spire_psw_filter(self):

        """
        This function ...
        :return:
        """

        return Filter("SPIRE PSW")

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
    if not fs.is_file(fitting_configuration_path): return None

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
