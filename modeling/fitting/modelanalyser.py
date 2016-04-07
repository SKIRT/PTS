#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelanalyser Contains the ModelAnalyser class, used for analysing the goodness
#  of the radiative transfer model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem, tables
from ..core.sed import ObservedSED

# -----------------------------------------------------------------

class ModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelAnalyser, self).__init__(config)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The flux calculator
        self.flux_calculator = None

        # The observed fluxes
        self.fluxes = None

        # The flux differences table
        self.differences = None

        # The calculated chi squared value
        self.chi_squared = None

    # -----------------------------------------------------------------

    def run(self, simulation, flux_calculator):

        """
        This function ...
        :param simulation:
        :param flux_calculator:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, flux_calculator)

        # 2. Load the observed SED
        self.load_observed_sed()

        # 3. Calculate the differences
        self.calculate_differences()

        # 4. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 5. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the model analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.flux_calculator = None
        self.fluxes = None
        self.differences = None
        self.chi_squared = None

    # -----------------------------------------------------------------

    def setup(self, simulation, flux_calculator):

        """
        This function ...
        :param simulation:
        :param flux_calculator:
        :return:
        """

        # Call the setup function of the base class
        super(ModelAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Make a local reference to the flux calculator
        if flux_calculator is None:
            raise RuntimeError("No ObservedFluxCalculator found; the calculate_observed_fluxes flag must be enabled on "
                               "each simulation that is part of the radiative transfer modeling")
        self.flux_calculator = flux_calculator

        # Initialize the differences table
        names = ["Instrument", "Band", "Flux difference", "Relative difference"]
        data = [[], [], [], []]
        dtypes = ["S5", "S7", "float64", "float64"]
        self.differences = tables.new(data, names, dtypes=dtypes)

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the fluxes table
        fluxes_path = filesystem.join(self.phot_path, "fluxes.dat")

        # Load the observed SED
        self.fluxes = ObservedSED.from_file(fluxes_path)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # In the flux-density tables derived from the simulation (created by the ObservedFluxCalculator object),
        # search the one corresponding to the "earth" instrument
        table_name = self.galaxy_name + "_earth"
        if table_name not in self.flux_calculator.tables: raise RuntimeError("Could not find a flux-density table for the 'earth' instrument")

        # Get the table
        table = self.flux_calculator.tables[table_name]

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(table)):

            #observatory = table["Observatory"][i]
            instrument = table["Instrument"][i]
            band = table["Band"][i]
            wavelength = table["Wavelength"][i]
            fluxdensity = table["Flux"][i]

            # Find the corresponding flux in the SED derived from observation
            observed_fluxdensity = self.fluxes.flux_for_band(instrument, band, unit="Jy").value

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + instrument + " " + band + " band")
                continue

            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Add an entry to the differences table
            self.differences.add_row([instrument, band, difference, relative_difference])

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        self.chi_squared = 1.

        return

        # From Sébastien for M31 ..

        refSED   = "Files/M31skirtSED_weight.dat"

        D = 2.4222569e22 # 0.785 Mpc in m
        Lsun = 3.846e26

        # Observed SED
        input   = np.loadtxt(refSED)
        obsWls  = input[:,0]
        obsFlux = input[:,1]
        obsErr  = input[:,2]
        chi2weight = input[:,3]

        obsFlux = obsFlux * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun
        obsErr  = obsErr * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun

        compareFlux = obsFlux[chi2weight>=0]
        compareErr  = obsErr[chi2weight>=0]

        # skirt SED
        input      = np.loadtxt(inpath+sed)
        modWls     = input[:,0]
        modFlux    = input[:,1]
        modDirect  = input[:,2]
        modStellarScatter = input[:,3]
        modDust    = input[:,4]
        modDustScatter = input[:,5]
        modTrans   = input[:,6]

        modFlux    = modFlux    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDirect  = modDirect  * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modStellarScatter = modStellarScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDust    = modDust    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modTrans   = modTrans   * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDustScatter = modDustScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun

        modBands = np.array([])
        for band in bands:
            filter = "Files/Filters/transmission_"+band+".dat"
            modBands = np.append(modBands, convolveFilter(modFlux,modWls,filter))

        self.chi_squared = np.sum(chi2weight[chi2weight>=0] * (compareFlux - modBands)**2/compareErr**2)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the flux differences
        self.write_differences()

        # Write the chi-squared value
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux-density differences for the current model ...")

        # Determine the path to the differences table
        path = filesystem.join(self.fit_res_path, self.simulation.name, "differences.dat")

        # Save the differences table
        tables.write(self.differences, path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding the chi squared value for the current model to the chi squared data file ...")

        # Open the chi squared file (in 'append' mode)
        resultfile = open(self.chi_squared_table_path, 'a')

        # Add a line to the chi squared file containing the simulation name and the chi squared value
        resultfile.write(self.simulation.name + " " + str(self.chi_squared) + "\n")

        # Close the output file
        resultfile.close()

# -----------------------------------------------------------------
