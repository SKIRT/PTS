#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.environment Contains the ModelingEnvironment class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .history import ModelingHistory
from ...core.basics.configuration import Configuration
from ...core.data.sed import ObservedSED
from ...core.filter.filter import parse_filter

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
        self.config_file_path = fs.join(self.path, "modeling.cfg")

        # Determine the path to the modeling history file
        self.history_file_path = fs.join(self.path, "history.dat")

        # Initialize the history file
        if not fs.is_file(self.history_file_path):
            history = ModelingHistory()
            history.saveto(self.history_file_path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.fit_path = fs.create_directory_in(self.path, "fit")
        self.analysis_path = fs.create_directory_in(self.path, "analysis")
        self.reports_path = fs.create_directory_in(self.path, "reports")
        self.visualisation_path = fs.create_directory_in(self.path, "visualisation")
        self.plot_path = fs.create_directory_in(self.path, "plot")
        self.log_path = fs.create_directory_in(self.path, "log")
        self.config_path = fs.create_directory_in(self.path, "config")
        self.show_path = fs.create_directory_in(self.path, "show")
        self.build_path = fs.create_directory_in(self.path, "build")

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
        self.sed_path = fs.join(self.path, "sed.dat")

        # Set the SED plot path
        self.sed_plot_path = fs.join(self.path, "sed.pdf")

        # Set the ski template path
        self.ski_path = fs.join(self.path, "template.ski")

        # Set the ski input path
        self.ski_input_path = fs.create_directory_in(self.path, "input")

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
        self.images_path = fs.create_directory_in(self.path, "images")

        # Set images header path
        self.images_header_path = fs.join(self.images_path, "header.txt")

        # Set the ski template path
        self.ski_path = fs.join(self.path, "template.ski")

        # Set the ski input path
        self.ski_input_path = fs.create_directory_in(self.path, "input")

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
