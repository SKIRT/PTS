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
from ....core.tools import filesystem as fs
from ....core.basics.log import log

# -----------------------------------------------------------------

class MapsAnalysisComponent(AnalysisComponent):
    
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
        super(MapsAnalysisComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

        # Set the paths to the raw maps subdirectories
        self.maps_colours_path = fs.create_directory_in(self.maps_raw_path, colours_name)
        self.maps_ssfr_path = fs.create_directory_in(self.maps_raw_path, ssfr_name)
        self.maps_tir_path = fs.create_directory_in(self.maps_raw_path, tir_name)
        self.maps_attenuation_path = fs.create_directory_in(self.maps_raw_path, attenuation_name)
        self.maps_old_path = fs.create_directory_in(self.maps_raw_path, old_name)
        self.maps_young_path = fs.create_directory_in(self.maps_raw_path, young_name)
        self.maps_ionizing_path = fs.create_directory_in(self.maps_raw_path, ionizing_name)
        self.maps_dust_path = fs.create_directory_in(self.maps_raw_path, dust_name)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MapsAnalysisComponent, self).setup(**kwargs)

        # Load the analysis run
        self.analysis_run = None

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
