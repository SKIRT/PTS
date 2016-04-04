#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedFluxCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import tables
from ..tools.logging import log
from ..basics.filter import Filter
from ...modeling.core.sed import SED

# -----------------------------------------------------------------

class ObservedFluxCalculator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ObservedFluxCalculator, self).__init__()

        # -- Attributes --

        # The paths to the SED files produced by SKIRT
        self.sed_paths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The filters for which the fluxes should be calculated
        self.filters = None

        # The fluxes table
        self.fluxes = None

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None, filter_names=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param filter_names:
        :return:
        """

        # Obtain the paths to the SED files created by the simulation
        self.sed_paths = simulation.seddatpaths()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Create the filters
        self.create_filters()

        # Calculate the observed fluxes
        self.calculate()

        # Write the results
        if output_path is not None: self.write(output_path)

    # -----------------------------------------------------------------

    def create_filters(self):

        """
        This function ...
        :return:
        """

        # Initialize the list
        self.filters = []

        # Loop over the different filter names
        for filter_name in self.filter_names:

            # Create the filter
            filter = Filter.from_string(filter_name)

            # Add the filter to the list
            self.filters.append(filter)

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # Loop over the different SEDs
        for sed_path in self.sed_paths:

            # Load the simulated SED
            sed = SED.from_file(sed_path)

            # Get the wavelengths and flux densities
            wavelengths = sed.wavelengths("micron", asarray=True)
            fluxdensities = sed.fluxes("Jy", asarray=True)

            # Loop over the different filters
            for filter in self.filters:

                # Calculate the flux
                flux = filter.convolve(wavelengths, fluxdensities)



    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Writing the observed flux table ...")

        # Write the fluxes table
        tables.write(self.fluxes, output_path)

# -----------------------------------------------------------------
