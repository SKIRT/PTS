#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.run Contains the FittingRun class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class FittingRun(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # The path to the template ski file
        self.template_ski_path = None

        # The path to the fixed parameter values file
        self.fixed_parameters_path = None

        # The path to the fit/evolution directory
        self.fit_generations_path = None

        # The path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = None

        # The path to the wavelength grids table
        self.wavelength_grids_table_path = None

        # The path to the fit/dust grids directory
        self.fit_dust_grids_path = None

        # The path to the dust grids table
        self.dust_grids_table_path = None

        # The path to the fit/best directory
        self.fit_best_path = None

        # The path to the fit/instruments directory
        self.fit_instruments_path = None

        # The paths to the SED and frame instrument files
        self.sed_instrument_path = None
        self.frame_instrument_path = None
        self.simple_instrument_path = None

        # The path to the fit/prob directory
        self.fit_prob_path = None

        # The path to the generations table
        self.generations_table_path = None

        # The path to the best parameters table
        self.best_parameters_table_path = None

        # The path to the weights table
        self.weights_table_path = None

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

        # The path to the geometries directory
        self.fit_geometries_path = None

        # The directory with the probability distributions for the different free parameters
        self.prob_distributions_path = None

        # The directory for the posterior / parameter probability distribution tables
        self.prob_parameters_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Set the path to the template ski file
        self.template_ski_path = fs.join(self.fit_path, "template.ski")

        # Set the path to the fixed parameters file
        self.fixed_parameters_path = fs.join(self.fit_path, "fixed.dat")

        # Set the path to the fit/generations directory
        self.fit_generations_path = fs.create_directory_in(self.fit_path, "generations")

        # Set the path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = fs.create_directory_in(self.fit_path, "wavelength grids")

        # Set the path to the wavelength grids table
        self.wavelength_grids_table_path = fs.join(self.fit_wavelength_grids_path, "grids.dat")

        # Set the path to the fit/dust grids directory
        self.fit_dust_grids_path = fs.create_directory_in(self.fit_path, "dust grids")

        # Set the path to the dust grids table
        self.dust_grids_table_path = fs.join(self.fit_dust_grids_path, "grids.dat")

        # Set the path to the fit/best directory
        self.fit_best_path = fs.create_directory_in(self.fit_path, "best")

        # Set the path to the fit/prob directory
        self.fit_prob_path = fs.create_directory_in(self.fit_path, "prob")

        # Set the path to the fit/instruments directory
        self.fit_instruments_path = fs.create_directory_in(self.fit_path, "instruments")

        # Set the path to the SED and frame instrument
        self.sed_instrument_path = fs.join(self.fit_instruments_path, "sed.instr")
        self.frame_instrument_path = fs.join(self.fit_instruments_path, "frame.instr")
        self.simple_instrument_path = fs.join(self.fit_instruments_path, "simple.instr")

        # Set the path to the fit/geometries directory
        self.fit_geometries_path = fs.create_directory_in(self.fit_path, "geometries")

        # -----------------------------------------------------------------

        ## WEIGHTS TABLE

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.fit_path, "weights.dat")

        ## TIMING TABLE

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.fit_path, "timing.dat")

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.fit_path, "memory.dat")

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

        ## GENERATIONS TABLE

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.fit_path, "generations.dat")

        # Initialize the generations table if necessary
        if not fs.is_file(self.generations_table_path) and self.free_parameter_labels is not None:
            generations_table = GenerationsTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            generations_table.saveto(self.generations_table_path)

        ## PROBABILITY DISTRIBUTION TABLES

        # The directory with the probability distributions for the different free parameters
        self.prob_distributions_path = fs.create_directory_in(self.fit_prob_path, "distributions")

        # The directory with the combined probability tables for the different free parameters
        self.prob_parameters_path = fs.create_directory_in(self.fit_prob_path, "parameters")

        ## BEST PARAMETERS TABLE

        # Set the path to the best parameters table
        self.best_parameters_table_path = fs.join(self.fit_path, "best_parameters.dat")

        # Initialize the best parameters table if necessary
        if not fs.is_file(self.best_parameters_table_path) and self.free_parameter_labels is not None:
            best_parameters_table = BestParametersTable(parameters=self.free_parameter_labels,
                                                        units=self.parameter_units)
            best_parameters_table.saveto(self.best_parameters_table_path)

        
# -----------------------------------------------------------------
