#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization.base Contains the FittingInitializerBase class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import FittingComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.prep.wavelengthgrids import WavelengthGridGenerator, WavelengthGridsTable
from ..tables import WeightsTable
from ....core.tools.serialization import write_dict
from ....core.prep.smile import SKIRTSmileSchema
from ...build.suite import ModelSuite
from ..weights import WeightsCalculator

# -----------------------------------------------------------------

class FittingInitializerBase(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FittingInitializerBase, self).__init__(*args, **kwargs)

        # The INITIAL model representation
        self.representation = None

        # The fitting run
        self.fitting_run = None

        # The ski file
        self.ski = None

        # The table of weights for each band
        self.weights = None

        # The wavelength grid generator
        self.wg_generator = None

        # The wavelength grid table
        self.wg_table = None

        # The dictionary of input map paths
        self.input_map_paths = dict()

        # The models suite
        self.suite = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        FittingComponent.setup(self, **kwargs)
        
        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Create the table to contain the weights
        self.weights = WeightsTable()

        # Set the models suite
        self.suite = ModelSuite.from_modeling_path(self.config.path)

        # Clear
        if self.has_weights_table: fs.remove_file(self.weights_table_path)
        if self.has_input_maps_file: fs.remove_file(self.input_maps_path)
        if self.config.regenerate_wavelength_grids:
            fs.clear_directory(self.wavelength_grids_path)
            if self.has_wavelength_grids_table: fs.remove_file(self.wavelength_grids_table_path)

        # Initialize the wavelength grids table
        if self.has_wavelength_grids_table: self.wg_table = WavelengthGridsTable.from_file(self.wavelength_grids_table_path)
        else: self.wg_table = WavelengthGridsTable()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_name

    # -----------------------------------------------------------------

    def load_representation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the model representation ...")

        # Load the initial representation
        self.representation = self.fitting_run.initial_representation

    # -----------------------------------------------------------------

    def set_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust grid ...")

        # Check whether we can use the file tree dust grid
        smile = SKIRTSmileSchema()
        if smile.supports_file_tree_grids and self.representation.has_dust_grid_tree:

            # Create file tree dust grid
            dust_grid = self.representation.create_file_tree_dust_grid(write=False)

            # Make sure it is only the file name, not a complete path
            dust_grid.filename = fs.name(dust_grid.filename)

        # Just take the real dust grid object
        else: dust_grid = self.representation.dust_grid

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(dust_grid)

    # -----------------------------------------------------------------

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Fixed wavelengths (always in the grid)
        fixed = [self.fuv_wavelength, self.i1_wavelength]

        # Set options
        self.wg_generator.config.show = False
        self.wg_generator.config.write = False

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range=self.config.wg.npoints_range, ngrids=self.config.wg.ngrids,
                              fixed=fixed, add_emission_lines=self.config.wg.add_emission_lines,
                              min_wavelength=self.config.wg.min_wavelength, max_wavelength=self.config.wg.max_wavelength,
                              filters=self.fitting_run.fitting_filters, table=self.wg_table)

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Calculate weights
        calculator = WeightsCalculator(self.config.weighing)
        calculator.config.write = False
        calculator.config.show = True
        calculator.run(filters=self.fitting_run.fitting_filters)

        # Set weights
        self.weights = calculator.table

    # -----------------------------------------------------------------

    @property
    def weights_table_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.weights_table_path

    # -----------------------------------------------------------------

    @property
    def has_weights_table(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.weights_table_path)

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights to " + self.weights_table_path + " ...")

        # Write the table with weights
        self.weights.saveto(self.weights_table_path)

    # -----------------------------------------------------------------

    @property
    def wavelength_grids_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.wavelength_grids_path

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Loop over the grids
        for index, grid in enumerate(self.wg_generator.grids):

            # Debugging
            log.debug("Writing wavelength grid " + str(index) + " ...")

            # Determine the path to the grid
            path = fs.join(self.wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

    # -----------------------------------------------------------------

    @property
    def wavelength_grids_table_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.wavelength_grids_table_path

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grids_table(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.wavelength_grids_table_path)

    # -----------------------------------------------------------------

    def write_wavelength_grid_table(self):

        """
        Tihs function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid table ...")

        # Write the wavelength grids table
        self.wg_table.saveto(self.wavelength_grids_table_path)

    # -----------------------------------------------------------------

    @property
    def input_maps_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.input_maps_file_path

    # -----------------------------------------------------------------

    @property
    def has_input_maps_file(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.input_maps_path)

    # -----------------------------------------------------------------

    def write_input_map_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dictionary of input map paths ...")

        # Write
        write_dict(self.input_map_paths, self.input_maps_path)

# -----------------------------------------------------------------
