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
from ...maps.collection import MapsCollection, StaticMapsCollection
from ....core.filter.filter import parse_filter

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

    def get_frame(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse filter
        fltr = parse_filter(fltr)

        # Get frame
        return self.get_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    @property
    def frame_list(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulated_frame_list

    #  -----------------------------------------------------------------

    @property
    def named_frame_list(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulated_named_frame_list

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
