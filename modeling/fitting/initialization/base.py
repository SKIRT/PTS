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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..component import FittingComponent
from ....core.tools import filesystem as fs
from ....magic.tools import wavelengths
from ....core.basics.log import log
from ....core.prep.wavelengthgrids import WavelengthGridGenerator, WavelengthGridsTable
from ..tables import WeightsTable
from ....core.tools.serialization import write_dict
from ....core.prep.smile import SKIRTSmileSchema
from ...build.suite import ModelSuite
from ....core.tools.utils import lazyproperty
from ....core.tools import sequences

# -----------------------------------------------------------------

wavelength_regimes = ["uv", "optical", "nir", "mir", "fir", "submm-microwave"]

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

    @lazyproperty
    def regimes(self):

        """
        This function ...
        :return:
        """

        # Only UV
        if self.config.only_uv:

            if self.config.no_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = ["uv"]

        # Only optical
        elif self.config.only_optical:

            if self.config.no_optical: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = ["optical"]

        # Only NIR
        elif self.config.only_nir:

            if self.config.no_nir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = ["nir"]

        # Only MIR
        elif self.config.only_mir:

            if self.config.no_mir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = ["mir"]

        # Only FIR
        elif self.config.only_fir:

            if self.config.no_fir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = ["fir"]

        # Only submm/microwave
        elif self.config.only_submm_microwave:

            if self.config.no_submm_microwave: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            regimes = ["submm-microwave"]

        # Regimes
        else: regimes = self.config.regimes[:]

        # Ignore certain regimes?
        if self.config.no_uv: regimes = sequences.removed_item(regimes, "uv")
        if self.config.no_optical: regimes = sequences.removed_item(regimes, "optical")
        if self.config.no_nir: regimes = sequences.removed_item(regimes, "nir")
        if self.config.no_mir: regimes = sequences.removed_item(regimes, "mir")
        if self.config.no_fir: regimes = sequences.removed_item(regimes, "fir")
        if self.config.no_submm_microwave: regimes = sequences.removed_item(regimes, "submm-microwave")

        # Check number of regimes
        if len(regimes) == 0: raise ValueError("No regimes")

        # Return the regimes
        return regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def uv_weight(self):

        """
        This function ...
        :return:
        """

        if "uv" in self.regimes: return self.config.uv
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def optical_weight(self):

        """
        This function ...
        :return:
        """

        if "optical" in self.regimes: return self.config.optical
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def nir_weight(self):

        """
        This function ...
        :return:
        """

        if "nir" in self.regimes: return self.config.nir
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def mir_weight(self):

        """
        This function ...
        :return:
        """

        if "mir" in self.regimes: return self.config.mir
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def fir_weight(self):

        """
        This function ...
        :return:
        """

        if "fir" in self.regimes: return self.config.fir
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def submm_microwave_weight(self):

        """
        This function ...
        :return:
        """

        if "submm-microwave" in self.regimes: return self.config.submm_microwave
        else: return 0.

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Get the weights
        weights = calculate_weights_filters(self.fitting_run.fitting_filters, uv=self.uv_weight, optical=self.optical_weight, nir=self.nir_weight, mir=self.mir_weight, fir=self.fir_weight, submm_microwave=self.submm_microwave_weight)

        # Add to weights table
        for fltr in weights: self.weights.add_point(fltr, weights[fltr])

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

def calculate_weights_filters_old(filters, uv=1, optical=1, nir=1, mir=1, fir=1, submm_microwave=1):

    """
    This function ...
    :param filters:
    :param uv:
    :param optical:
    :param nir:
    :param mir:
    :param fir:
    :param submm_microwave:
    :return:
    """

    # Get bands per regime
    uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands = split_filters_regimes_old(filters)

    # Get nbands per regime
    nuv = len(uv_bands)
    noptical = len(optical_bands)
    nnir = len(nir_bands)
    nmir = len(mir_bands)
    nfir = len(fir_bands)
    nsubmm_microwave = len(submm_microwave_bands)

    # Determine regime weights
    uv_weight, optical_weight, nir_weight, mir_weight, fir_weight, submm_radio_weight = calculate_weights(nuv, noptical, nnir, nmir, nfir, nsubmm_microwave, uv=uv, optical=optical, nir=nir, mir=mir, fir=fir, submm_microwave=submm_microwave)

    # Initialize dictionary for the weights per filter
    weights = OrderedDict()

    # Loop over the bands in each group and set the weight in the weights table
    for fltr in uv_bands: weights[fltr] = uv_weight
    for fltr in optical_bands: weights[fltr] = optical_weight
    for fltr in nir_bands: weights[fltr] = nir_weight
    for fltr in mir_bands: weights[fltr] = mir_weight
    for fltr in fir_bands: weights[fltr] = fir_weight
    for fltr in submm_microwave_bands: weights[fltr] = submm_radio_weight

    # Return the weights
    return weights

# -----------------------------------------------------------------

def calculate_weights_filters(filters, uv=1, optical=1, nir=1, mir=1, fir=1, submm_microwave=1):

    """
    This function ...
    :param filters:
    :param uv:
    :param optical:
    :param nir:
    :param mir:
    :param fir:
    :param submm_microwave:
    :return:
    """

    # Get bands per regime
    uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands = split_filters_regimes(filters)

    # Get nbands per regime
    nuv = len(uv_bands)
    noptical = len(optical_bands)
    nnir = len(nir_bands)
    nmir = len(mir_bands)
    nfir = len(fir_bands)
    nsubmm_microwave = len(submm_microwave_bands)

    # Determine regime weights
    uv_weight, optical_weight, nir_weight, mir_weight, fir_weight, submm_radio_weight = calculate_weights(nuv, noptical, nnir, nmir, nfir, nsubmm_microwave, uv=uv, optical=optical, nir=nir, mir=mir, fir=fir, submm_microwave=submm_microwave)

    # Initialize dictionary for the weights per filter
    weights = OrderedDict()

    # Loop over the bands in each group and set the weight in the weights table
    for fltr in uv_bands: weights[fltr] = uv_weight
    for fltr in optical_bands: weights[fltr] = optical_weight
    for fltr in nir_bands: weights[fltr] = nir_weight
    for fltr in mir_bands: weights[fltr] = mir_weight
    for fltr in fir_bands: weights[fltr] = fir_weight
    for fltr in submm_microwave_bands: weights[fltr] = submm_radio_weight

    # Return the weights
    return weights

# -----------------------------------------------------------------

def split_filters_regimes_old(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    # Initialize lists to contain the filters of the different wavelength ranges
    uv_bands = []
    optical_bands = []
    nir_bands = []
    mir_bands = []
    fir_bands = []
    submm_microwave_bands = []

    # Loop over the observed SED filters
    for fltr in filters:

        # Get the central wavelength
        wavelength = fltr.wavelength

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
        elif spectrum[0] == "Submm": submm_microwave_bands.append(fltr)
        elif spectrum[0] == "Radio" and spectrum[1] == "Microwave": submm_microwave_bands.append(fltr)
        else: raise RuntimeError("Unknown wavelength range: " + str(spectrum))

    # Return
    return uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands

# -----------------------------------------------------------------

def get_nbands_per_regime_old(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands = split_filters_regimes_old(filters)
    return len(uv_bands), len(optical_bands), len(nir_bands), len(mir_bands), len(fir_bands), len(submm_microwave_bands)

# -----------------------------------------------------------------

def calculate_weights_old(nuv, noptical, nnir, nmir, nfir, nsubmm_microwave, uv=1, optical=1, nir=1, mir=1, fir=1, submm_microwave=1):

    """
    This function ...
    :param nuv:
    :param noptical:
    :param nnir:
    :param nmir:
    :param nfir:
    :param nsubmm_microwave:
    :param uv:
    :param optical:
    :param nir:
    :param mir:
    :param fir:
    :param submm_microwave:
    :return:
    """

    # Set the number of groups
    number_of_groups = 0

    # Check which groups are present
    has_uv = nuv > 0
    has_optical = noptical > 0
    has_nir = nnir > 0
    has_mir = nmir > 0
    has_fir = nfir > 0
    has_submm_microwave = nsubmm_microwave > 0

    if has_uv: number_of_groups += 1
    if has_optical: number_of_groups += 1
    if has_nir: number_of_groups += 1
    if has_mir: number_of_groups += 1
    if has_fir: number_of_groups += 1
    if has_submm_microwave: number_of_groups += 1

    nall_groups = 6

    # Determine total number of data points
    number_of_data_points = nuv + noptical + nnir + nmir + nfir + nsubmm_microwave

    # Determine normalizations
    total_normalization = uv + optical + nir + mir + fir + submm_microwave
    uv = float(uv) / total_normalization * nall_groups
    optical = float(optical) / total_normalization * nall_groups
    nir = float(nir) / total_normalization * nall_groups
    mir = float(mir) / total_normalization * nall_groups
    fir = float(fir) / total_normalization * nall_groups
    submm = float(submm_microwave) / total_normalization * nall_groups

    # Determine the weight for each group of filters
    uv_weight = uv / (nuv * number_of_groups) * number_of_data_points if has_uv else 0.0
    optical_weight = optical / (noptical * number_of_groups) * number_of_data_points if has_optical else 0.0
    nir_weight = nir / (nnir * number_of_groups) * number_of_data_points if has_nir else 0.0
    mir_weight = mir / (nmir * number_of_groups) * number_of_data_points if has_mir else 0.0
    fir_weight = fir / (nfir * number_of_groups) * number_of_data_points if has_fir else 0.0
    submm_microwave_weight = submm / (nsubmm_microwave * number_of_groups) * number_of_data_points if has_submm_microwave else 0.0

    # Debugging
    if has_uv: log.debug("UV: number of bands = " + str(nuv) + ", weight = " + str(uv_weight))
    if has_optical: log.debug("Optical: number of bands = " + str(noptical) + ", weight = " + str(optical_weight))
    if has_nir: log.debug("NIR: number of bands = " + str(nnir) + ", weight = " + str(nir_weight))
    if has_mir: log.debug("MIR: number of bands = " + str(nmir) + ", weight = " + str(mir_weight))
    if has_fir: log.debug("FIR: number of bands = " + str(nfir) + ", weight = " + str(fir_weight))
    if has_submm_microwave: log.debug("Submm/microwave: number of bands = " + str(nsubmm_microwave) + ", weight = " + str(submm_microwave_weight))

    # Return the weights
    return uv_weight, optical_weight, nir_weight, mir_weight, fir_weight, submm_microwave_weight

# -----------------------------------------------------------------
