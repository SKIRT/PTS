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
from ..core.component import ModelingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem, tables

# -----------------------------------------------------------------

class ModelAnalyser(ModelingComponent):

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

        # The observed fluxes table
        self.fluxes = None

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

        # 2. Load the observed fluxes
        self.load_observed_fluxes()

        # 3. Calculate the differences
        self.calculate_differences()

        # 3. Calculate the chi squared for this model
        self.calculate_chi_squared()

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

    # -----------------------------------------------------------------

    def load_observed_fluxes(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the fluxes table
        fluxes_path = filesystem.join(self.phot_path, "fluxes.dat")

        # Load the fluxes table
        self.fluxes = tables.from_file(fluxes_path)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

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

        chi2 = np.sum(chi2weight[chi2weight>=0] * (compareFlux - modBands)**2/compareErr**2)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the model analyser ...")

        # Set the simulation to None
        self.simulation = None

# -----------------------------------------------------------------
