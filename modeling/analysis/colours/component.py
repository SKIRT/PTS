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

# -----------------------------------------------------------------

class ColourAnalysisComponent(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ColourAnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The path to the colours/observed directory
        self.colours_observed_path = None

        # The path to the colours/simulated directory
        self.colours_simulated_path = None

        # The path to the colours/residuals directory
        self.colours_residuals_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ColourAnalysisComponent, self).setup()

        # Set the path to the colours/observed directory
        self.colours_observed_path = fs.join(self.analysis_colours_path, "observed")

        # Set the path to the colours/simulated directory
        self.colours_simulated_path = fs.join(self.analysis_colours_path, "simulated")

        # Set the path to the colours/residuals directory
        self.colours_residuals_path = fs.join(self.analysis_colours_path, "residuals")

        # Create the directories
        fs.create_directories(self.colours_observed_path, self.colours_simulated_path, self.colours_residuals_path)

# -----------------------------------------------------------------
