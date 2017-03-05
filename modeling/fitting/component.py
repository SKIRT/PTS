#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from .tables import RunsTable
from .run import FittingRun, get_generation_names, get_finished_generations

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config, interactive)

        # -- Attributes --

        # Runs table path
        self.runs_table_path = None

        # The database path
        self.database_path = None

        # The statistics path
        self.statistics_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup(**kwargs)

        # Set the path to the runs table
        self.runs_table_path = fs.join(self.fit_path, "runs.dat")

        # Create the runs table if it doesn't exist yet
        if not fs.is_file(self.runs_table_path):
            table = RunsTable()
            table.saveto(self.runs_table_path)

        # Set the path to the database
        self.database_path = fs.join(self.fit_path, "database.db")

        # Set the path to the statistics file
        self.statistics_path = fs.join(self.fit_path, "statistics.csv")

    # -----------------------------------------------------------------

    def load_fitting_run(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        model_name = self.model_for_run(name)
        return FittingRun(self.config.path, name, model_name)

    # -----------------------------------------------------------------

    @property
    def runs_table(self):

        """
        This function ...
        :return:
        """

        return RunsTable.from_file(self.runs_table_path)

    # -----------------------------------------------------------------

    @property
    def run_names(self):

        """
        This function ...
        :return:
        """

        return self.runs_table.run_names

    # -----------------------------------------------------------------

    def model_for_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.runs_table.model_for_run(run_name)

    # -----------------------------------------------------------------

    @property
    def statistics(self):

        """
        This function ...
        :return:
        """

        return Table.read(self.statistics_path)

# -----------------------------------------------------------------

def get_run_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fit_path = fs.join(modeling_path, "fit")
    return fs.directories_in_path(fit_path, returns="name")

# -----------------------------------------------------------------

def get_run_generation_combinations(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    combinations = []

    for run_name in get_run_names(modeling_path):
        for generation in get_generation_names(modeling_path, run_name):
            combinations.append((run_name, generation))

    # Return the combinations
    return combinations

# -----------------------------------------------------------------

def get_run_finished_generation_combinations(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    combinations = []

    for run_name in get_run_names(modeling_path):
        for generation in get_finished_generations(modeling_path, run_name):
            combinations.append((run_name, generation))

    # Return the combinations
    return combinations

# -----------------------------------------------------------------

def load_fitting_run(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    model_name = get_model_for_run(modeling_path, name)
    return FittingRun(modeling_path, name, model_name)

# -----------------------------------------------------------------

def get_runs_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fit_path = fs.join(modeling_path, "fit")
    return fs.join(fit_path, "runs.dat")

# -----------------------------------------------------------------

def get_runs_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_runs_table_path(modeling_path)
    return RunsTable.from_file(path)

# -----------------------------------------------------------------

def get_model_for_run(modeling_path, name):

    """
    This function ...
    :return:
    """

    table = get_runs_table(modeling_path)
    return table.model_for_run(name)

# -----------------------------------------------------------------
