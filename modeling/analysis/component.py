#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from .context import AnalysisContext

# -----------------------------------------------------------------

class AnalysisComponent(GalaxyModelingComponent):
    
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
        super(AnalysisComponent, self).__init__(*args, **kwargs)

        # The analysis context
        self.context = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        # Create the analysis context
        self.context = AnalysisContext(self.analysis_path)

    # -----------------------------------------------------------------

    @property
    def timing_table_path(self):

        """
        This function ...
        :return:
        """

        return self.context.timing_table_path

    # -----------------------------------------------------------------

    @property
    def memory_table_path(self):

        """
        This function ...
        :return:
        """

        return self.context.memory_table_path

    # -----------------------------------------------------------------

    @property
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return self.context.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return self.context.memory_table

    # -----------------------------------------------------------------

    @property
    def cached_table(self):

        """
        This fucntion ...
        :return:
        """

        return self.context.cached_table

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return self.context.analysis_run_names

    # -----------------------------------------------------------------

    def get_run_path(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run_path(run_name)

    # -----------------------------------------------------------------

    def get_run_info(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run_info(run_name)

    # -----------------------------------------------------------------

    def get_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run(run_name)

# -----------------------------------------------------------------
