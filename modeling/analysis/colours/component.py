#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.colours.component Contains the ColourAnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log

# -----------------------------------------------------------------

simulated_name = "simulated"
observed_name = "observed"
residuals_name = "residuals"
distributions_name = "distributions"

# -----------------------------------------------------------------

class ColourAnalysisComponent(AnalysisComponent):
    
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
        super(ColourAnalysisComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    @property
    def colours_observed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.colours_path, observed_name)

    # -----------------------------------------------------------------

    @property
    def colours_simulated_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.colours_path, simulated_name)

    # -----------------------------------------------------------------

    @property
    def colours_residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.colours_path, residuals_name)

    # -----------------------------------------------------------------

    @property
    def colours_distributions_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.colours_path, distributions_name)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColourAnalysisComponent, self).setup()

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
