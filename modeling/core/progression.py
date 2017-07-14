#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.progression Contains the GalaxyModelingProgression class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .environment import GalaxyModelingEnvironment
from ..fitting.run import FittingRuns, FittingRun
from ..analysis.run import AnalysisRuns, AnalysisRun
from ...core.basics.configuration import prompt_string

# -----------------------------------------------------------------

def create_modeling_progression(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get analysis runs
    analysis_runs = AnalysisRuns(modeling_path)

    # Choose analysis run
    if analysis_runs.empty: analysis_run_name = None
    elif analysis_runs.has_single: analysis_run_name = analysis_runs.single_name
    else: analysis_run_name = prompt_string("analysis_run", "analysis run", choices=analysis_runs.names, default=analysis_runs.last_name, required=True)

    # Check whether fitting run has to be choosen
    analysis_run = AnalysisRun.from_name(modeling_path, analysis_run_name)
    if analysis_run_name is None or not analysis_run.from_fitting:

        # Get fitting runs
        fitting_runs = FittingRuns(modeling_path)

        # Choose fitting run
        if fitting_runs.empty: fitting_run_name = None
        elif fitting_runs.has_single: fitting_run_name = fitting_runs.single_name
        else: fitting_run_name =

    # Just take fitting run name for given analysis run
    else: fitting_run_name = analysis_run.fitting_run_name

    # Check whether model representation has to be choosen
    fitting_run = FittingRun.from_name(modeling_path, fitting_run_name)
    if fitting_run_name is None:


    else: representation_name = fitting_run.initial_representation_name

    # Check whether the model name has to be choosen
    if representation_name is None:

    else: model_name =



    # Return
    return GalaxyModelingProgression(modeling_path, model_name, representation_name, fitting_run_name, analysis_run_name)

# -----------------------------------------------------------------

class GalaxyModelingProgression(object):

    """
    This function ...
    """

    def __init__(self, modeling_path, model_name=None, representation_name=None, fitting_run_name=None, analysis_run_name=None):

        """
        This function ...
        :param modeling_path:
        :param model_name:
        :param representation_name:
        :param fitting_run_name:
        :param analysis_run_name:
        """

        # Set attributes
        self.modeling_path = modeling_path
        self.model_name = model_name
        self.representation_name = representation_name
        self.fitting_run_name = fitting_run_name
        self.analysis_run_name = analysis_run_name

    # -----------------------------------------------------------------

    @lazyproperty
    def environment(self):

        """
        This function ...
        :return:
        """

        return GalaxyModelingEnvironment(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def representation(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ..
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
