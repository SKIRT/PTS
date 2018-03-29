#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.context Contains the FittingContext class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .tables import RunsTable
from ...evolve.optimize.stepwise import load_populations
from .run import get_generation_names, get_finished_generations
from ...core.tools import filesystem as fs
from .run import FittingRuns
from ...core.tools import tables

# -----------------------------------------------------------------

runs_filename = "runs.dat"
database_filename = "database.db"
statistics_filename = "statistics.csv"
populations_filename = "populations.dat"
earth_instrument_name = "earth"

# -----------------------------------------------------------------

class FittingContext(object):

    """
    This class ...
    """

    def __init__(self, fit_path):

        """
        This function ...
        :param fit_path:
        """

        # Set the fitting path
        self.path = fit_path

        # Set the path to the runs table
        self.runs_table_path = fs.join(self.path, runs_filename)

        # Create the runs table if it doesn't exist yet
        if not fs.is_file(self.runs_table_path):

            table = RunsTable()
            table.saveto(self.runs_table_path)

        # Set the path to the database
        self.database_path = fs.join(self.path, database_filename)

        # Set the path to the statistics file
        self.statistics_path = fs.join(self.path, statistics_filename)

        # Set the path to the populations file
        self.populations_path = fs.join(self.path, populations_filename)

        # The name of the instrument
        self.earth_instrument_name = earth_instrument_name

        # The fitting runs
        self.runs = FittingRuns(self.modeling_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_modeling_path(cls, modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        fit_path = fs.join(modeling_path, "fit")
        return cls(fit_path)

    # -----------------------------------------------------------------

    @property
    def fit_path(self):

        """
        This function ...
        :return:
        """

        return self.path

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def generation_combinations(self):

        """
        This function ...
        :return:
        """

        combinations = []

        for run_name in self.runs.names:

            for generation in get_generation_names(self.modeling_path, run_name):
                combinations.append((run_name, generation))

        # Return the combinations
        return combinations

    # -----------------------------------------------------------------

    @property
    def finished_generation_combinations(self):

        """
        This function ...
        :return:
        """

        combinations = []

        for run_name in self.runs.names:
            for generation in get_finished_generations(self.modeling_path, run_name):
                combinations.append((run_name, generation))

        # Return the combinations
        return combinations

    # -----------------------------------------------------------------

    @property
    def fitting_run_names(self):

        """
        This function ...
        :return:
        """

        return self.runs.names

    # -----------------------------------------------------------------

    def load_fitting_run(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.runs.load(name)

    # -----------------------------------------------------------------

    @property
    def runs_table(self):

        """
        This function ...
        :param self:
        :return:
        """

        return RunsTable.from_file(self.runs_table_path)

    # -----------------------------------------------------------------

    def get_model_for_run(self, name):

        """
        This function ...
        :return:
        """

        return self.runs_table.model_for_run(name)

    # -----------------------------------------------------------------

    @property
    def statistics(self):

        """
        Thisf ucntion ...
        :param self:
        :return:
        """

        from ...evolve.analyse.statistics import load_statistics
        return load_statistics(self.statistics_path)

    # -----------------------------------------------------------------

    @property
    def database(self):

        """
        This function ...
        :return:
        """

        from ...evolve.analyse.database import load_database
        return load_database(self.database_path)

    # -----------------------------------------------------------------

    @property
    def populations(self):

        """
        This function ...
        :return:
        """

        return load_populations(self.populations_path)

    # -----------------------------------------------------------------

    def get_populations_for_run(self, run_id):

        """
        This function ...
        :param run_id:
        :return:
        """

        return self.populations[run_id]

    # -----------------------------------------------------------------

    def get_statistics_for_run(self, run_id):

        """
        This function ...
        :param run_id:
        :return:
        """

        return tables.filtered(self.statistics, "Identifier", run_id)

# -----------------------------------------------------------------

def get_runs_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fit_path = fs.join(modeling_path, "fit")
    filepath = fs.join(fit_path, runs_filename)
    return RunsTable.from_file(filepath)

# -----------------------------------------------------------------

def get_model_name_for_run(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    runs_table = get_runs_table(modeling_path)
    return runs_table.model_for_run(name)

# -----------------------------------------------------------------
