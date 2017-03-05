#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization.sed Contains the SEDFittingInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import FittingComponent
from ....core.tools import tables
from ....core.tools import filesystem as fs
from ....magic.tools import wavelengths
from ....core.tools.logging import log
from ....core.prep.wavelengthgrids import WavelengthGridGenerator
from ...component.sed import SEDModelingComponent
from ....core.simulation.skifile import LabeledSkiFile
from ..tables import WeightsTable
from ....core.basics.unit import parse_unit as u

# -----------------------------------------------------------------

class SEDFittingInitializer(FittingComponent, SEDModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructors of the base classes
        FittingComponent.__init__(self, config, interactive)
        SEDModelingComponent.__init__(self, config, interactive)

        # The fitting run
        self.fitting_run = None

        # The table of weights for each band
        self.weights = None

        # The wavelength grid generator
        self.wg_generator = None

        # The ski file
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski template
        self.load_ski_template()

        # 3. Create the wavelength grids
        self.create_wavelength_grids()

        # 5. Adjust the ski template
        self.adjust_ski()

        # 6. Calculate the weight factor to give to each band
        self.calculate_weights()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup functions of the base classes
        FittingComponent.setup(self, **kwargs)
        SEDModelingComponent.setup(self, **kwargs)

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the table to contain the weights
        self.weights = WeightsTable()

    # -----------------------------------------------------------------

    def load_ski_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the ski file
        self.ski = LabeledSkiFile(self.ski_path)

    # -----------------------------------------------------------------

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Get fixed wavelengths
        fixed = self.sed_filter_wavelengths

        # Set options
        self.wg_generator.config.show = False
        self.wg_generator.config.write = False

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range=self.config.wg.npoints_range, ngrids=self.config.wg.ngrids,
                              add_emission_lines=self.config.wg.add_emission_lines,
                              min_wavelength=self.config.wg.min_wavelength,
                              max_wavelength=self.config.wg.max_wavelength, fixed=fixed, filters=self.fitting_filters)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # 1. Set the instrument
        self.set_instrument()

        # 2. Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # 3. Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # 4. Set parameter values
        self.set_initial_parameter_values()

        # 5. Set dust emissivity
        self.set_dust_emissivity()

        # 6. Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # 7. Dust self-absorption
        self.set_selfabsorption()

        # 8. Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the instrument ...")

        instrument_names = self.ski.get_instrument_names()
        if len(instrument_names) > 1: raise ValueError("The ski file can only have one instrument")
        old_name = instrument_names[0]
        self.ski.set_instrument_name(old_name, "earth")

    # -----------------------------------------------------------------

    def set_initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the initial parameter values ...")

        # Loop over the free parameter labels
        for label in self.free_parameter_labels:

            # Determine value as the mean of the given range
            value = self.free_parameter_ranges[label].mean

            # Set the value for the parameter
            self.ski.set_labeled_value(label, value)

    # -----------------------------------------------------------------

    def set_dust_emissivity(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust emissivity ...")

        # Set dust emissivity if present
        if self.ski.has_dust_emissivity:

            # Enable or disable
            if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
            else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust self-absorption ...")

        # Set
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Initialize lists to contain the filters of the different wavelength ranges
        uv_bands = []
        optical_bands = []
        nir_bands = []
        mir_bands = []
        fir_bands = []
        submm_bands = []

        # Loop over the observed SED filters
        for fltr in self.fitting_filters:

            # Get the central wavelength
            wavelength = fltr.center

            # Get a string identifying which portion of the wavelength spectrum this wavelength belongs to
            spectrum = wavelengths.name_in_spectrum(wavelength)

            # Determine to which group
            if spectrum[0] == "UV": uv_bands.append(fltr)
            elif spectrum[0] == "Optical": optical_bands.append(fltr)
            elif spectrum[0] == "Optical/IR": optical_bands.append(fltr)
            elif spectrum[0] == "IR":
                if spectrum[1] == "NIR": nir_bands.append(fltr)
                elif spectrum[1] == "MIR": mir_bands.append(fltr)
                elif spectrum[1] == "FIR": fir_bands.append(fltr)
                else: raise RuntimeError("Unknown IR range")
            elif spectrum[0] == "Submm": submm_bands.append(fltr)
            else: raise RuntimeError("Unknown wavelength range")

        # Set the number of groups
        number_of_groups = 0

        # Check wihch groups are present
        has_uv = len(uv_bands) > 0
        has_optical = len(optical_bands) > 0
        has_nir = len(nir_bands) > 0
        has_mir = len(mir_bands) > 0
        has_fir = len(fir_bands) > 0
        has_submm = len(submm_bands) > 0

        if has_uv: number_of_groups += 1
        if has_optical: number_of_groups += 1
        if has_nir: number_of_groups += 1
        if has_mir: number_of_groups += 1
        if has_fir: number_of_groups += 1
        if has_submm: number_of_groups += 1

        # Detemrine total number of data points
        number_of_data_points = len(self.fitting_filters)

        # Determine the weight for each group of filters
        uv_weight = 1. / (len(uv_bands) * number_of_groups) * number_of_data_points if has_uv else 0.0
        optical_weight = 1. / (len(optical_bands) * number_of_groups) * number_of_data_points if has_optical else 0.0
        nir_weight = 1. / (len(nir_bands) * number_of_groups) * number_of_data_points if has_nir else 0.0
        mir_weight = 1. / (len(mir_bands) * number_of_groups) * number_of_data_points if has_mir else 0.0
        fir_weight = 1. / (len(fir_bands) * number_of_groups) * number_of_data_points if has_fir else 0.0
        submm_weight = 1. / (len(submm_bands) * number_of_groups) * number_of_data_points if has_submm else 0.0

        # Debugging
        if has_uv: log.debug("UV: number of bands = " + str(len(uv_bands)) + ", weight = " + str(uv_weight))
        if has_optical: log.debug("Optical: number of bands = " + str(len(optical_bands)) + ", weight = " + str(optical_weight))
        if has_nir: log.debug("NIR: number of bands = " + str(len(nir_bands)) + ", weight = " + str(nir_weight))
        if has_mir: log.debug("MIR: number of bands = " + str(len(mir_bands)) + ", weight = " + str(mir_weight))
        if has_fir: log.debug("FIR: number of bands = " + str(len(fir_bands)) + ", weight = " + str(fir_weight))
        if has_submm: log.debug("Submm: number of bands = " + str(len(submm_bands)) + ", weight = " + str(submm_weight))

        # Loop over the bands in each group and set the weight in the weights table
        for fltr in uv_bands: self.weights.add_point(fltr, uv_weight)
        for fltr in optical_bands: self.weights.add_point(fltr, optical_weight)
        for fltr in nir_bands: self.weights.add_point(fltr, nir_weight)
        for fltr in mir_bands: self.weights.add_point(fltr, mir_weight)
        for fltr in fir_bands: self.weights.add_point(fltr, fir_weight)
        for fltr in submm_bands: self.weights.add_point(fltr, submm_weight)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the ski file
        self.write_ski()

        # 2. Write the weights table
        self.write_weights()

        # 3. Write the wavelength grids
        self.write_wavelength_grids()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.template_ski_path + " ...")

        # Save the ski template file
        self.ski.saveto(self.template_ski_path)

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights to " + self.fitting_run.weights_table_path + " ...")

        # Write the table with weights
        self.weights.saveto(self.fitting_run.weights_table_path)

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Loop over the grids
        index = 0
        for grid in self.wg_generator.grids:

            # Determine the path to the grid
            path = fs.join(self.fitting_run.wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

            # Increment the index
            index += 1

        # Write the wavelength grids table
        tables.write(self.wg_generator.table, self.fitting_run.wavelength_grids_table_path)

# -----------------------------------------------------------------
