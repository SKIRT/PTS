#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.stars.ionizing Contains the IonizingStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...core.tools import filesystem as fs
from ...core.basics.distribution import Distribution
from ...core.plot.distribution import DistributionPlotter
from ...magic.core.image import Image
from ...core.units.parsing import parse_unit as u
from ...magic.maps.ionizingstars.ionizing import IonizingStellarMapsMaker

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

class IonizingStellarMapMaker(MapsComponent):

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

        # Call the constructor of the base class
        super(IonizingStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The MIPS 24 micron image and error map IN SOLAR UNITS
        self.mips24 = None
        self.mips24_errors = None

        # The Halpha map IN SOLAR UNITS
        self.halpha = None

        # The NORMALIZED (to unity) DISK IMAGE
        self.disk = None

        # The maps of the corrected 24 micron emission
        self.corrected_24mu_maps = dict()

        # The distributions of the pixel values of the corrected 24mu maps
        self.corrected_24mu_distributions = dict()

        # The map of ionizing stars
        #self.map = None

        # NEW: THE MAPS OF IONIZING STARS
        self.maps = dict()

        # The path to the maps/ionizing/24mu directory
        self.maps_ionizing_24mu_path = None

        # The path to the maps/ionizing/solar directory
        self.maps_ionizing_solar_path = None

        # The region of the pixels used for plotting the distributions of pixel values
        self.distribution_region = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary frames
        self.load_frames()

        # Calculate the significance masks
        #self.calculate_significance()

        # 3. Make the map
        self.make_map()

        # ...
        self.create_distribution_region()
        self.make_distributions()

        # 4. Normalize the map
        #self.normalize_map()

        # Make the cutoff mask
        self.make_cutoff_mask()

        # 5. Cut-off map
        #self.cutoff_map()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(IonizingStellarMapMaker, self).setup()

        # Set paths
        #self.maps_ionizing_24mu_path = fs.create_directory_in(self.maps_ionizing_path, "24mu")
        #self.maps_ionizing_solar_path = fs.create_directory_in(self.maps_ionizing_path, "solar")

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the MIPS 24 micron image and convert to solar units
        self.load_mips()

        # Load the H alpha image and convert to solar units
        self.load_halpha()

        # Load the disk image and normalize to unity
        self.load_disk()

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.mips24_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("MIPS 24mu", self.config.mips24_significance), "MIPS_24mu")
        if self.config.halpha_significance > 0: self.significance.add_mask(self.get_halpha_significance_mask(self.config.halpha_significance), "Halpha")

    # -----------------------------------------------------------------



    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image
        self.halpha = self.masked_halpha_frame

        # Convert from erg/s to Lsun
        self.halpha.convert_to("Lsun")

    # -----------------------------------------------------------------




    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps that are converted to solar units
        self.write_solar()

        # Write the maps
        self.write_24mu_maps()

        # Write distribution region
        self.write_distribution_region()

        # Write histograms of corrected 24 micron pixels
        self.write_24mu_histograms()

        # Write the terms in the equation to make the map, normalized
        #self.write_terms()

        # Write the ionizing stars map
        #self.write_map()

        # Write the significance masks
        self.write_significance_masks()

    # -----------------------------------------------------------------

    def write_solar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps that have been converted to solar units ...")

        # Save
        path = fs.join(self.maps_ionizing_solar_path, "MIPS 24mu.fits")
        self.mips24.saveto(path)

        # Save H alpha in solar units
        path = fs.join(self.maps_ionizing_solar_path, "Halpha.fits")
        self.halpha.saveto(path)

    # -----------------------------------------------------------------

    def write_24mu_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the corrected 24 micron images ...")

        # Loop over the corrected 24 micron maps
        for factor in self.corrected_24mu_maps:

            # Determine path
            path = fs.join(self.maps_ionizing_24mu_path, str(factor) + ".fits")

            # Write
            self.corrected_24mu_maps[factor].saveto(path)

    # -----------------------------------------------------------------

    def write_distribution_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution region ...")

        path = fs.join(self.maps_ionizing_24mu_path, "histogram.reg")
        self.distribution_region.saveto(path)

    # -----------------------------------------------------------------

    def write_24mu_histograms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the histograms of the corrected 24 micron pixels in the specified region ...")

        # Create a distribution plotter
        plotter = DistributionPlotter()

        # Loop over the distributions
        for factor in self.corrected_24mu_distributions:

            # Determine path
            path = fs.join(self.maps_ionizing_24mu_path, str(factor) + " histogram.pdf")

            # Plot the distribution as a histogram
            plotter.add_distribution(self.corrected_24mu_distributions[factor], "Correction factor of " + str(factor))
            plotter.run(path)

            # Clear the distribution plotter
            plotter.clear()

    # -----------------------------------------------------------------

    def write_terms(self):

        """
        This function ...
        :return:
        """

        # BOTH TERMS ARE JUST SUMMMED WITH EQUAL WEIGHT (AFTER 0.031 FACTOR), SO IT MAKES SENSE TO COMPARE THE HISTOGRAMS

        # Calculate the terms

        # H alpha
        term_halpha = self.halpha

        # 24 micron
        term_24mu = 0.031 * self.corrected_24mu_maps[self.config.best_factor]

        # Determine paths
        halpha_path = fs.join(self.maps_ionizing_path, "halpha_term.fits")
        mips_24_path = fs.join(self.maps_ionizing_path, "MIPS24_term.fits")

        # Save the terms
        term_halpha.saveto(halpha_path)
        term_24mu.saveto(mips_24_path)

        # Distributions
        halpha_distribution = Distribution.from_values(term_halpha.data.flatten())
        mips24_distribution = Distribution.from_values(term_24mu.data.flatten())

        plotter = DistributionPlotter()

        plotter.add_distribution(halpha_distribution, "H alpha contribution")
        plotter.add_distribution(mips24_distribution, "corrected MIPS 24 contribution")

        # Run the plotter
        path = fs.join(self.maps_ionizing_path, "histograms.pdf")
        plotter.run(path)

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Write the ionizing stars map ...")

        # Write
        self.map.saveto(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.saveto(self.ionizing_stellar_significance_path)

# -----------------------------------------------------------------
