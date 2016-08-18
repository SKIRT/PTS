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

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .tables import GenerationsTable
from ...core.simulation.skifile import SkiFile, LabeledSkiFile
from ...core.basics.distribution import Distribution
from ..basics.instruments import load_instrument

# -----------------------------------------------------------------

contributions = ["old", "young", "ionizing"]

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        # The path to the template ski file
        self.template_ski_path = None

        # The path to the fit/evolution directory
        self.fit_generations_path = None

        # The path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = None

        # The path to the fit/dust grids directory
        self.fit_dust_grids_path = None

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

        # The path to the weights table
        self.weights_table_path = None

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

        # The path to the geometries directory
        self.fit_geometries_path = None

        # The paths of the directories of the simulations that calculate the contributions of the various stellar populations
        self.fit_best_contribution_paths = dict()

        # The path of the directory to generate simulated images for the best model
        self.fit_best_total_path = None

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

        # Set the path to the template ski file
        self.template_ski_path = fs.join(self.fit_path, "template.ski")

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

        # Set the path to the fit/instruments directory
        self.fit_instruments_path = fs.create_directory_in(self.fit_path, "instruments")

        # Set the path to the SED and frame instrument
        self.sed_instrument_path = fs.join(self.fit_instruments_path, "sed.instr")
        self.frame_instrument_path = fs.join(self.fit_instruments_path, "frame.instr")
        self.simple_instrument_path = fs.join(self.fit_instruments_path, "simple.instr")

        # Set the path to the fit/geometries directory
        self.fit_geometries_path = fs.create_directory_in(self.fit_path, "geometries")

        # Set and create the paths to the fit/best/ contribution directories
        for contribution in contributions:

            path = fs.join(self.fit_best_path, contribution)
            fs.create_directory(path)
            self.fit_best_contribution_paths[contribution] = path

        # Set the path to the fit/best/total directory
        self.fit_best_total_path = fs.create_directory_in(self.fit_best_path, "total")

        # -----------------------------------------------------------------

        ## WEIGHTS TABLE

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.fit_path, "weights.dat")

        ## TIMING TABLE

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.fit_path, "timing.dat")

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):
            timing_table = TimingTable.initialize()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.fit_path, "memory.dat")

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable.initialize()
            memory_table.saveto(self.memory_table_path)

        ## GENERATIONS TABLE

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.fit_path, "generations.dat")

        # Initialize the generations table if necessary
        if not fs.is_file(self.generations_table_path) and self.free_parameter_labels is not None:
            generations_table = GenerationsTable.initialize(self.free_parameter_labels)
            generations_table.saveto(self.generations_table_path)

        # Set the paths to the probability distribution tables
        if self.free_parameter_labels is not None:
            for label in self.free_parameter_labels:

                # Determine the path to the table
                path = fs.join(self.fit_prob_path, label + ".dat")

                # Set the path
                self.distribution_table_paths[label] = path

    # -----------------------------------------------------------------

    @lazyproperty
    def generations_table(self):

        """
        This function ...
        :return:
        """

        return GenerationsTable.from_file(self.generations_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.generation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        return LabeledSkiFile(self.template_ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values(self):

        """
        This function ...
        :return:
        """

        # TODO: get parameter values of best fitting model

        # Get the current values in the ski file prepared by InputInitializer
        young_luminosity_guess, young_filter = self.ski_template.get_stellar_component_luminosity("Young stars")
        ionizing_luminosity_guess, ionizing_filter = self.ski_template.get_stellar_component_luminosity("Ionizing stars")
        dust_mass_guess = self.ski_template.get_dust_component_mass(0)

    # -----------------------------------------------------------------

    def get_parameter_distribution(self, label, normalized=True):

        """
        This function ...
        :param label:
        :param normalized:
        :return:
        """

        # Load the probability distribution
        distribution = Distribution.from_file(self.distribution_table_paths[label])

        # Normalize the distribution
        if normalized: distribution.normalize(value=1.0, method="max")

        # Return the distribution
        return distribution

    # -----------------------------------------------------------------

    def has_distribution(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return fs.is_file(self.distribution_table_paths[label])

    # -----------------------------------------------------------------

    @lazyproperty
    def last_generation_index(self):

        """
        This function ...
        :return:
        """

        highest_index = -1

        # Find the highest index
        for i in range(len(self.generations_table)):
            if not self.generations_table["Generation index"].mask[i]:
                index = self.generations_table["Generation index"][i]
                if index  > highest_index: highest_index = index

        # Return the highest generation index
        return highest_index

    # -----------------------------------------------------------------

    @lazyproperty
    def ngenerations(self):

        """
        This function ...
        :return:
        """

        return len(self.generations_table)

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_wavelength_grid_level(self):

        """
        This function ...
        :return:
        """

        # Return the last filename, sorted as integers
        return int(fs.files_in_path(self.fit_wavelength_grids_path, not_contains="grids", extension="txt", returns="name", sort=int)[-1])

    # -----------------------------------------------------------------

    @lazyproperty
    def current_wavelength_grid_level(self):

        """
        This function ...
        :return:
        """

        if len(self.generations_table) > 0: return self.generations_table["Wavelength grid level"][-1]
        else: return 0

    # -----------------------------------------------------------------

    def wavelength_grid_path_for_level(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        return fs.join(self.fit_wavelength_grids_path, str(level) + ".txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_dust_grid_level(self):

        """
        This function ...
        :return:
        """

        # Return the last filename, sorted as integers
        return int(fs.files_in_path(self.fit_dust_grids_path, not_contains="grids", extension="dg", returns="name", sort=int)[-1])

    # -----------------------------------------------------------------

    @lazyproperty
    def current_dust_grid_level(self):

        """
        This function ...
        :return:
        """

        if len(self.generations_table) > 0: return self.generations_table["Dust grid level"][-1]
        else: return 0

    # -----------------------------------------------------------------

    def dust_grid_path_for_level(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        return fs.join(self.fit_dust_grids_path, str(level) + ".dg")

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.sed_instrument_path): raise RuntimeError("The SED instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.sed_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.frame_instrument_path): raise RuntimeError("The frame instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.frame_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.simple_instrument_path): raise RuntimeError("The simple instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.simple_instrument_path)

# -----------------------------------------------------------------
