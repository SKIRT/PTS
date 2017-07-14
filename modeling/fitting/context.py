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
from .run import FittingRun, get_generation_names, get_finished_generations
from ...core.tools import filesystem as fs

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

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

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

    def get_fit_path(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        return fs.join(modeling_path, "fit")

    # -----------------------------------------------------------------

    def get_statistics_path(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        return fs.join(get_fit_path(modeling_path), "statistics.csv")

    # -----------------------------------------------------------------

    def get_statistics(modeling_path):

        """
        Thisf ucntion ...
        :param modeling_path:
        :return:
        """

        from ...evolve.analyse.statistics import load_statistics
        return load_statistics(get_statistics_path(modeling_path))

    # -----------------------------------------------------------------

    def get_database_path(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        return fs.join(get_fit_path(modeling_path), "database.db")

    # -----------------------------------------------------------------

    def get_database(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        from ...evolve.analyse.database import load_database
        return load_database(get_database_path(modeling_path))

    # -----------------------------------------------------------------

    def get_populations_path(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        return fs.join(get_fit_path(modeling_path), "populations.dat")

    # -----------------------------------------------------------------

    def get_populations(modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        return load_populations(get_populations_path(modeling_path))

# -----------------------------------------------------------------
