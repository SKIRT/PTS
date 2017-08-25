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

# Import the relevant PTS classes and modules
from .environment import GalaxyModelingEnvironment
from ..fitting.run import FittingRuns, FittingRun
from ..analysis.run import AnalysisRuns, AnalysisRun
from ..build.suite import ModelSuite
from ...core.basics.configuration import prompt_string
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

def create_modeling_progression(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Check modeling path
    if modeling_path is None: raise ValueError("Modeling path is not specified")

    # Get analysis runs
    analysis_runs = AnalysisRuns(modeling_path)

    # Choose analysis run
    if analysis_runs.empty: analysis_run_name = None
    elif analysis_runs.has_single: analysis_run_name = analysis_runs.single_name
    else: analysis_run_name = prompt_string("analysis_run", "name of the analysis run", choices=analysis_runs.names, default=analysis_runs.last_name, required=True)

    # Check whether fitting run has to be choosen
    analysis_run = AnalysisRun.from_name(modeling_path, analysis_run_name) if analysis_run_name is not None else None
    if analysis_run_name is None or not analysis_run.from_fitting:

        # Get fitting runs
        fitting_runs = FittingRuns(modeling_path)

        # Choose fitting run
        if fitting_runs.empty: fitting_run_name = None
        elif fitting_runs.has_single: fitting_run_name = fitting_runs.single_name
        else: fitting_run_name = prompt_string("fitting_run", "name of the fitting run", choices=fitting_runs.names, required=True)

    # Just take fitting run name for given analysis run
    elif analysis_run is not None: fitting_run_name = analysis_run.fitting_run_name
    else: raise RuntimeError("We shouldn't get here")

    # Get the model suite
    suite = ModelSuite.from_modeling_path(modeling_path)

    # Check whether model representation has to be choosen
    fitting_run = FittingRun.from_name(modeling_path, fitting_run_name) if fitting_run_name is not None else None
    if fitting_run_name is None:

        # Choose representation
        if suite.no_representations: representation_name = None
        elif suite.has_single_representation: representation_name = suite.single_representation_name
        else: representation_name = prompt_string("representation", "name of the model representation", choices=suite.representation_names, required=True)

    # Else
    elif fitting_run is not None: representation_name = fitting_run.initial_representation_name
    else: raise RuntimeError("We shouldn't get here")

    # Check whether the model name has to be choosen
    if representation_name is None:

        # No models
        if suite.no_models: model_name = None

        # Choose from models
        else: model_name = prompt_string("model", "name of the model", choices=suite.model_names, required=True)

    # Get the model name for the chosen representation
    else: model_name = suite.get_model_name_for_representation(representation_name)

    #print("MODEL NAME", model_name)
    #print("REPRESENTATION NAME", representation_name)
    #print("FITTING RUN NAME", fitting_run_name)
    #print("ANALYSIS RUN NAME", analysis_run_name)

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
    def model_suite(self):

        """
        This function ...
        :return:
        """

        if self.model_name is None: return None # no models yet or none choosen: don't load the model suite
        return ModelSuite.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_context(self):

        """
        This function ...
        :return:
        """

        from ..fitting.context import FittingContext
        if self.fitting_run_name is None: return None # no fitting runs yet or none choosen: don't load the fitting context
        return FittingContext.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):

        """
        This function ...
        :return:
        """

        return self.fitting_context.runs

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_context(self):

        """
        This function ...
        :return:
        """

        from ..analysis.context import AnalysisContext
        if self.analysis_run_name is None: return None
        return AnalysisContext.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def analysis_runs(self):

        """
        This function ...
        :return:
        """

        return self.analysis_context.runs

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

        if self.model_name is None: return None
        return self.model_suite.get_model_definition(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def representation(self):

        """
        This function ...
        :return:
        """

        if self.representation_name is None: return None
        return self.model_suite.get_model_representation(self.representation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ..
        :return:
        """

        if self.fitting_run_name is None: return None
        return self.fitting_runs.load(self.fitting_run_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """



# -----------------------------------------------------------------
