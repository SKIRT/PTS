#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.reporting.component Contains the ReportingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ReportingComponent(GalaxyModelingComponent):
    
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
        super(ReportingComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The path to the data report
        self.data_report_path = None

        # The path to the preparation initialization report
        self.preparation_initialization_report_path = None

        # The path to the preparation report
        self.preparation_report_path = None

        # The path to the decomposition report
        self.decomposition_report_path = None

        # The path to the photometry report
        self.photometry_report_path = None

        # The path to the map making report
        self.map_making_report_path = None

        # The path to the input initialization report
        self.input_initialization_report_path = None

        # The path to the exploration report
        self.exploration_report_path = None

        # The path to the fitting report
        self.fitting_report_path = None

        # The path to the analysis report
        self.analysis_report_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ReportingComponent, self).setup(**kwargs)

        # Set the report paths
        self.data_report_path = fs.join(self.reports_path, "data.txt")
        self.preparation_initialization_report_path = fs.join(self.reports_path, "preparation_initialization.txt")
        self.preparation_report_path = fs.join(self.reports_path, "preparation.txt")
        self.decomposition_report_path = fs.join(self.reports_path, "decomposition.txt")
        self.photometry_report_path = fs.join(self.reports_path, "photometry.txt")
        self.map_making_report_path = fs.join(self.reports_path, "map_making.txt")
        self.input_initialization_report_path = fs.join(self.reports_path, "input_initialization.txt")
        self.exploration_report_path = fs.join(self.reports_path, "exploration.txt")
        self.fitting_report_path = fs.join(self.reports_path, "fitting.txt")
        self.analysis_report_path = fs.join(self.reports_path, "analysis.txt")

# -----------------------------------------------------------------
