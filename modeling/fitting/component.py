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
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...core.launch.timing import TimingTable

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

        # The names of the fit parameters
        self.parameter_names = ["FUV young", "FUV ionizing", "Dust mass"]

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

        # The path to the timing table
        self.timing_table_path = None

        # The path to the scripts directory
        self.fit_scripts_path = None

        # The paths to the probability distribution tables
        self.distribution_table_paths = dict()

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
        self.fit_in_path = fs.join(self.fit_path, "in")

        # Set the path to the fit/out directory
        self.fit_out_path = fs.join(self.fit_path, "out")

        # Set the path to the fit/res directory
        self.fit_res_path = fs.join(self.fit_path, "res")

        # Set the path to the fit/plot directory
        self.fit_plot_path = fs.join(self.fit_path, "plot")

        # Set the path to the fit/best directory
        self.fit_best_path = fs.join(self.fit_path, "best")

        # Set the path to the fit/prob directory
        self.fit_prob_path = fs.join(self.fit_path, "prob")

        # Create the fit/in, fit/out, fit/res, fit/plot, fit/best and fit/prob directories
        fs.create_directories([self.fit_in_path, self.fit_out_path, self.fit_res_path, self.fit_plot_path, self.fit_best_path, self.fit_prob_path])

        # Set the path to the fit/scripts directory
        self.fit_scripts_path = fs.join(self.fit_path, "scripts")

        # Create the fit/scripts directory
        if not fs.is_directory(self.fit_scripts_path): fs.create_directory(self.fit_scripts_path)

        # Set the path to the parameter file
        self.parameter_table_path = fs.join(self.fit_path, "parameters.dat")

        # Initialize the parameter file if that hasn't been done yet
        if not fs.is_file(self.parameter_table_path):

            # Create the table
            names = ["Simulation name", "FUV young", "FUV ionizing", "Dust mass"]
            data = [[], [], [], []]
            dtypes = ["S24", "float64", "float64", "float64"]
            table = tables.new(data, names, dtypes=dtypes)

            # Create the (empty) table
            tables.write(table, self.parameter_table_path, format="ascii.ecsv")

        # Determine the path to the ski file
        self.fit_ski_path = fs.join(self.fit_path, self.galaxy_name + ".ski")

        # Set the path to the chi squared table file
        self.chi_squared_table_path = fs.join(self.fit_path, "chi_squared.dat")

        # Initialize the chi squared file if that hasn't been done yet
        if not fs.is_file(self.chi_squared_table_path):

            # Create the table
            names = ["Simulation name", "Chi squared"]
            data = [[], []]
            dtypes = ["S24", "float64"]
            table = tables.new(data, names, dtypes=dtypes)

            # Write the (empty) table
            tables.write(table, self.chi_squared_table_path, format="ascii.ecsv")

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.fit_path, "weights.dat")

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.fit_path, "timing.dat")

        # Initialize the timing table
        timing_table = TimingTable(self.timing_table_path)

        # Set the paths to the probability distribution tables
        for parameter_name in self.parameter_names:

            # Determine the path to the table
            path = fs.join(self.fit_prob_path, parameter_name.lower() + ".dat")

            # Set the path
            self.distribution_table_paths[parameter_name] = path

# -----------------------------------------------------------------
