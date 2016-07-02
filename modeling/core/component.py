#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.component Contains the ModelingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.basics.configurable import NewConfigurable
from ...core.tools import inspection
from ...core.tools import filesystem as fs
from ..core.sed import ObservedSED

# -----------------------------------------------------------------

class ModelingComponent(NewConfigurable):
    
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
        self.galaxy_name = None

        # Modeling directories
        self.data_path = None
        self.prep_path = None
        self.truncation_path = None
        self.phot_path = None
        self.maps_path = None
        self.components_path = None
        self.fit_path = None
        self.analysis_path = None
        self.reports_path = None
        self.visualisation_path = None
        self.plot_path = None
        self.log_path = None

        # PTS directories
        self.kernels_path = None

        # Reference image
        self.reference_image = "Pacs red"

        # The path to the observed SEDs
        self.observed_sed_path = None
        self.observed_sed_dustpedia_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Get the full paths to the necessary subdirectories
        self.data_path = fs.join(self.config.path, "data")
        self.prep_path = fs.join(self.config.path, "prep")
        self.truncation_path = fs.join(self.config.path, "truncated")
        self.phot_path = fs.join(self.config.path, "phot")
        self.maps_path = fs.join(self.config.path, "maps")
        self.components_path = fs.join(self.config.path, "components")
        self.fit_path = fs.join(self.config.path, "fit")
        self.analysis_path = fs.join(self.config.path, "analysis")
        self.reports_path = fs.join(self.config.path, "reports")
        self.visualisation_path = fs.join(self.config.path, "visualisation")
        self.plot_path = fs.join(self.config.path, "plot")
        self.log_path = fs.join(self.config.path, "log")

        # Determine the path to the kernels user directory
        self.kernels_path = fs.join(inspection.pts_user_dir, "kernels")

        # Check whether the 'data' directory exists, otherwise exit with an error
        if fs.is_directory(self.data_path):

            # Create the prep path if it does not exist yet
            fs.create_directories([self.prep_path, self.truncation_path, self.maps_path, self.phot_path,
                                   self.maps_path, self.components_path, self.fit_path, self.analysis_path,
                                   self.reports_path, self.visualisation_path, self.plot_path, self.log_path])

        # Exit with an error
        else: raise ValueError("The current working directory is not a radiative transfer modeling directory (the data directory is missing)")

        # Set the path to the observed SED
        self.observed_sed_path = fs.join(self.phot_path, "fluxes.dat")

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return ObservedSED.from_file(self.observed_sed_path)

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
