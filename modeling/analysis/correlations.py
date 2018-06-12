#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.correlations Contains the CorrelationsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log

# -----------------------------------------------------------------

class CorrelationsAnalyser(DustHeatingAnalysisComponent):
    
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
        super(CorrelationsAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the data points
        self.get_data()

        # 4. Writing
        self.write()

        # 5. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CorrelationsAnalyser, self).setup()

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        #
        self.get_ssfr_funev()

    # -----------------------------------------------------------------

    def get_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        ssfr = self.model.total_ssfr_map_earth
        funev =

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        self.plot_ssfr_funev()

# -----------------------------------------------------------------
