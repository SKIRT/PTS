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
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from ....core.tools import filesystem as fs
from ....core.basics.distribution import Distribution
from ....core.plot.distribution import DistributionPlotter

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * Unit("W")

# -----------------------------------------------------------------

class IonizingStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(IonizingStellarMapMaker, self).__init__(config)

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

        # The map of ionizing stars
        self.map = None

        # The path to the maps/ionizing/24mu directory
        self.maps_ionizing_24mu_path = None

        # The path to the maps/ionizing/solar directory
        self.maps_ionizing_solar_path = None

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

        # 3. Make the map
        self.make_map()

        # 4. Normalize the map
        self.normalize_map()

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

        # Set
        self.maps_ionizing_24mu_path = fs.create_directory_in(self.maps_ionizing_path, "24mu")

        self.maps_ionizing_solar_path = fs.create_directory_in(self.maps_ionizing_path, "solar")

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

    def load_mips(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the MIPS 24 micron image and converting to solar units ...")

        # Get MIPS 24 micron frame and error map
        self.mips24 = self.dataset.get_frame("MIPS 24mu")  # in original MJy/sr units
        self.mips24_errors = self.dataset.get_errors("MIPS 24mu")  # in original MJy/sr units

        ## CONVERT TO LSUN

        # Get the galaxy distance
        distance = self.galaxy_properties.distance

        # Get pixelscale and wavelength
        pixelscale = self.mips24.average_pixelscale
        wavelength = self.mips24.filter.pivotwavelength() * Unit("micron")

        # Conversion from MJy / sr to Jy / sr
        conversion_factor = 1e6

        # Conversion from Jy / sr to Jy / pix(2)
        conversion_factor *= (pixelscale ** 2).to("sr/pix2").value

        # Conversion from Jy / pix to W / (m2 * Hz) (per pixel)
        conversion_factor *= 1e-26

        # Conversion from W / (m2 * Hz) (per pixel) to W / (m2 * m) (per pixel)
        conversion_factor *= (speed_of_light / wavelength ** 2).to("Hz/m").value

        # Conversion from W / (m2 * m) (per pixel) [SPECTRAL FLUX] to W / m [SPECTRAL LUMINOSITY]
        conversion_factor *= (4. * np.pi * distance ** 2).to("m2").value

        # Conversion from W / m [SPECTRAL LUMINOSITY] to W [LUMINOSITY]
        conversion_factor *= wavelength.to("m").value

        # Conversion from W to Lsun
        conversion_factor *= 1. / solar_luminosity.to("W").value

        ## DO THE CONVERSION

        self.mips24 *= conversion_factor
        self.mips24_errors *= conversion_factor

    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image
        self.halpha = self.halpha_frame

        # Convert from erg/s to Lsun
        self.halpha.convert_to("Lsun")

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Get disk frame
        self.disk = self.disk_frame

        # Normalize the disk image
        self.disk.normalize()
        self.disk.unit = None

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # H-ALPHA HAS BEEN CONVERTED TO LSUN (ABOVE)

        # Inform the user
        log.info("Making the map of ionizing stars ...")

        # Loop over the different colour options
        for factor in (self.config.factor_range.linear(self.config.factor_nvalues, as_list=True) + [self.config.best_factor]):

            # Calculate the corrected 24 micron image
            corrected = self.make_corrected_24mu_map(factor)

            # Add the attenuation map to the dictionary
            self.corrected_24mu_maps[factor] = corrected

        # Young ionizing stars = Ha + 0.031 x MIPS24_corrected

        # Calculate ionizing stars map and ratio
        ionizing = self.halpha + 0.031 * self.corrected_24mu_maps[self.config.best_factor]

        # ionizing_ratio = self.ha / (0.031*mips_young_stars)

        # MASK NEGATIVE AND LOW SIGNAL-TO-NOISE PIXELS

        # Set pixels to zero with low signal-to-noise in the H Alpha image
        # ionizing[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0
        # ionizing_ratio[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0

        # Set pixels to zero with low signal-to-noise in the 24 micron image
        # ionizing[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0
        # ionizing_ratio[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        # ionizing[ionizing < 0.0] = 0.0
        # ionizing_ratio[ionizing < 0.0] = 0.0

        # New
        #ionizing[self.cutoff_masks["Halpha"]] = 0.0

        # Set the ionizing stars map
        self.map = ionizing

    # -----------------------------------------------------------------

    def make_corrected_24mu_map(self, factor):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the 24 micron emission map with a factor of " + str(factor) + " ...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        ## MIPS HAS BEEN CONVERTED TO LSUN (ABOVE)

        # typisch 20% en 35% respectievelijk
        # 48% voor MIPS 24 komt van Lu et al. 2014

        # Total contribution in solar units
        total_contribution = factor * self.mips24.sum()

        # Subtract the disk contribution to the 24 micron image
        new_mips = self.mips24 - total_contribution * self.disk # disk image is normalized

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_mips[new_mips < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        # new_mips[self.mips < self.config.ionizing_stars.mips_young_stars.mips_snr_level*self.mips_errors] = 0.0

        # Return the new 24 micron frame
        return new_mips

    # -----------------------------------------------------------------

    def normalize_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the map of ionizing stars ...")

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

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

        # Write the terms in the equation to make the map, normalized
        self.write_terms()

        # Write the ionizing stars map
        self.write_map()

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
        self.mips24.save(path)

        # Save H alpha in solar units
        path = fs.join(self.maps_ionizing_solar_path, "Halpha.fits")
        self.halpha.save(path)

    # -----------------------------------------------------------------

    def write_24mu_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Write the corrected 24 micron images ...")

        # Loop over the corrected 24 micron maps
        for factor in self.corrected_24mu_maps:

            # Determine path
            path = fs.join(self.maps_ionizing_24mu_path, str(factor) + ".fits")

            # Write
            self.corrected_24mu_maps[factor].save(path)

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
        term_halpha.save(halpha_path)
        term_24mu.save(mips_24_path)

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
        self.map.save(self.ionizing_stellar_map_path)

# -----------------------------------------------------------------
