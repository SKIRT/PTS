#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem, tables

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
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
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        # The path to the fit/in directory
        self.fit_in_path = None

        # The path to the fit/out directory
        self.fit_out_path = None

        # The path to the fit/res directory
        self.fit_res_path = None

        # The path to the fit/plot directory
        self.fit_plot_path = None

        # The path to the fit/best directory
        self.fit_best_path = None

        # The path to the fit/prob directory
        self.fit_prob_path = None

        # The path to the ski file
        self.fit_ski_path = None

        # The path to the parameter table
        self.parameter_table_path = None

        # The path to the chi squared table
        self.chi_squared_table_path = None

        # The path to the weights table
        self.weights_table_path = None

        # The path to the runtime table
        self.runtime_table_path = None

        # The path to the scripts directory
        self.fit_scripts_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup()

        # Set the output path
        self.config.output_path = self.fit_path

        # Set the path to the fit/in directory
        self.fit_in_path = filesystem.join(self.fit_path, "in")

        # Set the path to the fit/out directory
        self.fit_out_path = filesystem.join(self.fit_path, "out")

        # Set the path to the fit/res directory
        self.fit_res_path = filesystem.join(self.fit_path, "res")

        # Set the path to the fit/plot directory
        self.fit_plot_path = filesystem.join(self.fit_path, "plot")

        # Set the path to the fit/best directory
        self.fit_best_path = filesystem.join(self.fit_path, "best")

        # Set the path to the fit/prob directory
        self.fit_prob_path = filesystem.join(self.fit_path, "prob")

        # Create the fit/in, fit/out, fit/res, fit/plot, fit/best and fit/prob directories
        filesystem.create_directories([self.fit_in_path, self.fit_out_path, self.fit_res_path, self.fit_plot_path, self.fit_best_path, self.fit_prob_path])

        # Set the path to the fit/scripts directory
        self.fit_scripts_path = filesystem.join(self.fit_path, "scripts")

        # Create the fit/scripts directory
        if not filesystem.is_directory(self.fit_scripts_path): filesystem.create_directory(self.fit_scripts_path)

        # Set the path to the parameter file
        self.parameter_table_path = filesystem.join(self.fit_path, "parameters.dat")

        # Initialize the parameter file if that hasn't been done yet
        if not filesystem.is_file(self.parameter_table_path):

            # Create an empty table
            names = ["Simulation name", "FUV young", "FUV ionizing", "Dust mass"]
            data = [[], [], [], []]
            dtypes = ["S24", "float64", "float64", "float64"]
            table = tables.new(data, names, dtypes=dtypes)
            tables.write(table, self.parameter_table_path, format="ascii.ecsv")

        # Determine the path to the ski file
        self.fit_ski_path = filesystem.join(self.fit_path, self.galaxy_name + ".ski")

        # Set the path to the chi squared table file
        self.chi_squared_table_path = filesystem.join(self.fit_path, "chi_squared.dat")

        # Initialize the chi squared file if that hasn't been done yet
        if not filesystem.is_file(self.chi_squared_table_path):

            # Initialize
            names = ["Simulation name", "Chi squared"]
            data = [[], []]
            dtypes = ["S24", "float64"]
            table = tables.new(data, names, dtypes=dtypes)
            tables.write(table, self.chi_squared_table_path, format="ascii.ecsv")

        # Set the path to the weights table file
        self.weights_table_path = filesystem.join(self.fit_path, "weights.dat")

        # Set the path to the runtime table file
        self.runtime_table_path = filesystem.join(self.fit_path, "runtimes.dat")

        # Initialize the runtime file if that hasn't been done yet
        if not filesystem.is_file(self.runtime_table_path):

            # Initialize
            names = ["Simulation name", "Host id", "Cluster name", "Cores", "Hyperthreads per core", "Processes", "Packages", "Runtime"]
            data = [[], [], [], [], [], [], [], []]
            dtypes = ["S24", "S15", "S15", "int64", "int64", "int64", "int64", "float64"]
            table = tables.new(data, names, dtypes=dtypes)
            table["Runtime"] = "s" # runtimes are expressed in seconds
            tables.write(table, self.runtime_table_path, format="ascii.ecsv")

# -----------------------------------------------------------------
