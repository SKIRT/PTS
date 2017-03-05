#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.component Contains the DustHeatingAnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.tools.logging import log

# -----------------------------------------------------------------

contributions = ["total", "old", "young", "ionizing", "unevolved"]

# -----------------------------------------------------------------

class DustHeatingAnalysisComponent(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(DustHeatingAnalysisComponent, self).__init__(config, interactive)

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    @property
    def cell_heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.heating_path, "cell")

    # -----------------------------------------------------------------

    @property
    def projected_heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.heating_path, "projected")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustHeatingAnalysisComponent, self).setup(**kwargs)

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
