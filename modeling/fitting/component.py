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

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .generations import GenerationsTable

# -----------------------------------------------------------------

contributions = ["old", "young", "ionizing"]

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

        # The path to the fit/evolution directory
        self.fit_generations_path = None

        # The path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = None

        # The path to the fit/dust grids directory
        self.fit_dust_grids_path = None

        # The path to the fit/best directory
        self.fit_best_path = None

        # The path to the fit/prob directory
        self.fit_prob_path = None

        # The path to the fit/grid directory
        self.fit_grid_path = None

        # The path to the ski file
        self.fit_ski_path = None

        # The path to the generations table
        self.generations_table_path = None

        # The path to the parameter table
        self.parameter_table_path = None

        # The path to the chi squared table
        self.chi_squared_table_path = None

        # The path to the weights table
        self.weights_table_path = None

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

        # The path to the scripts directory
        self.fit_scripts_path = None

        # The path to the geometries directory
        self.fit_geometries_path = None

        # The paths of the directories of the simulations that calculate the contributions of the various stellar populations
        self.fit_best_contribution_paths = dict()

        # The path of the directory to generate simulated images for the best model
        self.fit_best_images_path = None

        # The paths to the probability distribution tables
        self.distribution_table_paths = dict()

        # The path to the reference image
        self.reference_path = None

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
        self.fit_in_path = fs.create_directory_in(self.fit_path, "in")

        # Set the path to the fit/generations directory
        self.fit_generations_path = fs.create_directory_in(self.fit_path, "generations")

        # Set the path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = fs.create_directory_in(self.fit_path, "wavelength grids")

        # Set the path to the fit/dust grids directory
        self.fit_dust_grids_path = fs.create_directory_in(self.fit_path, "dust grids")

        # Set the path to the fit/best directory
        self.fit_best_path = fs.create_directory_in(self.fit_path, "best")

        # Set the path to the fit/prob directory
        self.fit_prob_path = fs.create_directory_in(self.fit_path, "prob")

        # Set the path to the fit/grid directory
        self.fit_grid_path = fs.create_directory_in(self.fit_path, "grid")

        # Set the path to the fit/scripts directory
        self.fit_scripts_path = fs.create_directory_in(self.fit_path, "scripts")

        # Set the path to the fit/geometries directory
        self.fit_geometries_path = fs.create_directory_in(self.fit_path, "geometries")

        # Set and create the paths to the fit/best/ contribution directories
        for contribution in contributions:

            path = fs.join(self.fit_best_path, contribution)
            fs.create_directory(path)
            self.fit_best_contribution_paths[contribution] = path

        # Set the path to the fit/best/images directory
        self.fit_best_images_path = fs.create_directory_in(self.fit_best_path, "images")

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.fit_path, "generations.dat")

        # Initialize the generations table
        generations_table = GenerationsTable(self.generations_table_path)

        # -----------------------------------------------------------------

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

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.fit_path, "memory.dat")

        # Initialize the memory table
        memory_table = MemoryTable(self.memory_table_path)

        # Set the paths to the probability distribution tables
        for parameter_name in self.parameter_names:

            # Determine the path to the table
            path = fs.join(self.fit_prob_path, parameter_name.lower() + ".dat")

            # Set the path
            self.distribution_table_paths[parameter_name] = path

        # Set the path to the reference image
        self.reference_path = fs.join(self.truncation_path, self.reference_image + ".fits")

# -----------------------------------------------------------------
