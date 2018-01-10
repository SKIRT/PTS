#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.refitter Contains the Refitter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from .initialization.base import calculate_weights_filters
from ...core.tools.utils import lazyproperty
from .tables import WeightsTable
from ...core.tools import filesystem as fs
from ...core.tools import tables
from .modelanalyser import FluxDifferencesTable
from ...core.tools import sequences

# -----------------------------------------------------------------

earth_instrument_name = "earth"

# -----------------------------------------------------------------

class Refitter(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(Refitter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The table of weights for each band
        self.weights = None

        # The flux differences
        self.differences = defaultdict(lambda: defaultdict)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the weights
        self.get_weights()

        # 2. Get the differences
        self.get_differences()

        # 3. Calculate the chi squared values
        self.calculate_chi_squared()

        # 2. Get the parameters of the best models for each generation
        self.get_best_parameters()

        # 3. Calculate the probabilities
        self.calculate_probabilities()

        # 4. Calculate the probability distributions
        self.create_distributions()

        # 3. Writing
        self.write()

        # Show
        if self.config.show: self.show()

        # Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    @lazyproperty
    def path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.fitting_run.refitting_path, self.config.name)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Refitter, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.fitting_run)

        # Create the table to contain the weights
        self.weights = WeightsTable()

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
            if self.config.only_submm: raise ValueError("Error")
            regimes = ["uv"]

        # Only optical
        elif self.config.only_optical:

            if self.config.no_optical: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm: raise ValueError("Error")
            regimes = ["optical"]

        # Only NIR
        elif self.config.only_nir:

            if self.config.no_nir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm: raise ValueError("Error")
            regimes = ["nir"]

        # Only MIR
        elif self.config.only_mir:

            if self.config.no_mir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm: raise ValueError("Error")
            regimes = ["mir"]

        # Only FIR
        elif self.config.only_fir:

            if self.config.no_fir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_submm: raise ValueError("Error")
            regimes = ["fir"]

        # Only submm
        elif self.config.only_submm:

            if self.config.no_submm: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            regimes = ["submm"]

        # Regimes
        else: regimes = self.config.regimes[:]

        # Ignore certain regimes?
        if self.config.no_uv: regimes = sequences.removed_item(regimes, "uv")
        if self.config.no_optical: regimes = sequences.removed_item(regimes, "optical")
        if self.config.no_nir: regimes = sequences.removed_item(regimes, "nir")
        if self.config.no_mir: regimes = sequences.removed_item(regimes, "mir")
        if self.config.no_fir: regimes = sequences.removed_item(regimes, "fir")
        if self.config.no_submm: regimes = sequences.removed_item(regimes, "submm")

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
    def submm_weight(self):

        """
        This function ...
        :return:
        """

        if "submm" in self.regimes: return self.config.submm
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        if self.config.generations is not None: return self.config.generations
        else: return self.fitting_run.generation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def generations(self):

        """
        This function ...
        :return:
        """

        gens = OrderedDict()
        for name in self.generation_names: gens[name] = self.fitting_run.get_generation(name)
        return gens

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        if self.config.filters is not None: return self.config.filters
        else: return self.fitting_run.fitting_filters

    # -----------------------------------------------------------------

    @property
    def different_filters(self):

        """
        This function ...
        :return:
        """

        if self.config.filters is not None: return sequences.same_contents(self.config.filters, self.fitting_run.fitting_filters)
        else: return False

    # -----------------------------------------------------------------

    @property
    def reweigh(self):

        """
        This function ...
        :return:
        """

        if self.config.reweigh is not None: return self.config.reweigh
        else: return self.different_filters

    # -----------------------------------------------------------------

    def get_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting weights ...")

        # Recalculate the weights
        if self.reweigh: self.calculate_weights()

        # Load the original weights
        else: self.load_weights()

    # -----------------------------------------------------------------

    def load_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the weights to give to each band ...")

        # Load the weights
        self.weights = self.fitting_run.weights

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Get the weights
        weights = calculate_weights_filters(self.filters, uv=self.uv_weight, optical=self.optical_weight, nir=self.nir_weight, mir=self.mir_weight, fir=self.fir_weight, submm=self.submm_weight)

        # Add to weights table
        for fltr in weights: self.weights.add_point(fltr, weights[fltr])

    # -----------------------------------------------------------------

    @property
    def rediff(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.rediff is not None: return self.config.rediff
        else: return self.different_filters

    # -----------------------------------------------------------------

    def get_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting differences ...")

        # Recalculate the differences
        if self.rediff: self.calculate_differences()

        # Load the original differences
        else: self.load_differences()

    # -----------------------------------------------------------------

    def load_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the differences between observed and simulated fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Loading the differences for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulation names
            for simulation_name in generation.simulation_names:

                # Get the differences filepath
                if not generation.has_sed_differences(simulation_name): raise IOError("Differences file is not found for simulation '" + simulation_name + "'")
                differences = generation.get_simulation_sed_differences(simulation_name)

                # Set table
                self.differences[generation_name][simulation_name] = differences

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between observed and simulated fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Calculating differences for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulations
            for simulation in generation.analysed_simulations_basic:

                # Get simulation name
                simulation_name = simulation.name

                # Initialize the differences table
                differences = FluxDifferencesTable()

                # Get mock SED
                if not generation.has_mock_sed(simulation_name): raise IOError("Mock SED is not found for simulation '" + simulation_name + "'")
                mock_sed = generation.get_simulation_mock_sed(simulation_name)

                # Loop over the entries in the fluxdensity table (SED) derived from the simulation
                for i in range(len(mock_sed)):

                    # Get instrument, band and flux density
                    instrument = mock_sed["Instrument"][i]
                    band = mock_sed["Band"][i]
                    fluxdensity = mock_sed["Photometry"][i]

                    # Find the corresponding flux in the SED derived from observation
                    observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy").value

                    # Find the corresponding flux error in the SED derived from observation
                    observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value

                    # If no match with (instrument, band) is found in the observed SED
                    if observed_fluxdensity is None:
                        log.warning("The observed flux density could not be found for the " + instrument + " " + band + " band")
                        continue

                    difference = fluxdensity - observed_fluxdensity
                    relative_difference = difference / observed_fluxdensity

                    # Find the index of the current band in the weights table
                    index = tables.find_index(self.weights, key=[instrument, band], column_name=["Instrument", "Band"])
                    if index is None: continue  # Skip this band if a weight is not found
                    weight = self.weights["Weight"][index]

                    # Calculate the chi squared term
                    chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

                    # Add entry to the table
                    differences.add_entry(instrument, band, difference, relative_difference, chi_squared_term)

                # Set table
                self.differences[generation_name][simulation_name] = differences

    # -----------------------------------------------------------------

    @property
    def nfree_parameters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.nfree_parameters

    # -----------------------------------------------------------------

    # @property
    # def ndof(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return self.ndifferences - self.nfree_parameters - 1 # number of data points - number of fitted parameters - 1

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating chi squared values ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Loop over the simulations
            for simulation_name in self.differences[generation_name]:

                # Get the differences table
                differences = self.differences[generation_name][simulation_name]

                # Inform the user
                log.info("Calculating the chi squared value for this model ...")

                # The (reduced) chi squared value is the sum of all the terms (for each band),
                # divided by the number of degrees of freedom
                chi_squared = np.sum(differences["Chi squared term"]) / self.ndof

                # Debugging
                log.debug("Found a (reduced) chi squared value of " + str(chi_squared))

    # -----------------------------------------------------------------

    def get_best_parameters(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def calculate_probabilities(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Weights
        self.write_weights()

    # -----------------------------------------------------------------

    @property
    def weights_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "weights.dat")

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

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def plot(self):


        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
