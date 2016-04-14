#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.attenuation Contains the AttenuationAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent

# -----------------------------------------------------------------

# IR-through-UV extinction curve (Fitzpatrick+, 2007):
# http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJ/663/320
#  - J/ApJ/663/320/stars: Survey stars (tables 1, 3 and 4 of paper) (328 rows)

# UVOT imaging of M81 and Holmberg IX (Hoversten+, 2011):
#  - J/AJ/141/205/table3: Source photometry and fitted parameters
#    http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/AJ/141/205/table3
#  -> has "AV": Best fitting internal extinction A_V

# -----------------------------------------------------------------

class AttenuationAnalyser(AnalysisComponent):
    
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
        super(AttenuationAnalyser, self).__init__(config)

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new AttenuationAnalyser instance
        analyser = cls()

        # Set the modeling path
        analyser.config.path = arguments.path

        # Return the new instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

# -----------------------------------------------------------------
