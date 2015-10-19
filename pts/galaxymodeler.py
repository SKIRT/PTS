#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.galaxymodeler Model a galaxy by fitting using Astromagic and SKIRT
#
# An instance of the GalaxyModeler class in this module is responsible for taking reduced astronomical image data
# of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
# dust, star formation and old stars. Then, it uses these maps as input for SKIRT radiative transfer simulations
# and fits the output to the observed galaxy SED.

# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import numpy as np
import inspect

# Import astronomical modules
from astropy import units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic.tools import configuration
from astromagic import Image
import astromagic.utilities as iu
from astromagic.magic.galaxyextraction import GalaxyExtractor
from astromagic.magic.starextraction import StarExtractor
from astromagic.magic.skyextraction import SkyExtractor

# Import PTS modules
from pts.imagepreparation import ImagePreparation

# *****************************************************************

class GalaxyModeler(object):

    """
    An instance of the GalaxyModeler class in this module is responsible for taking reduced astronomical image data
    of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
    dust, star formation and old stars. Then, it uses these maps as input for SKIRT radiative transfer simulations
    and fits the output to the observed galaxy SED.
    """

    # *****************************************************************

    def __init__(self, path, filter_name=None, config=None):

        """
        The constructor ...
        :param directory:
        :param filter_name:
        :param plot:
        :param save:
        :param config:
        :return:
        """

        ### LOAD CONFIGURATION

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "galaxymodeler.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ### TEMPORARY

        self.config.decompose = False
        self.config.make_maps = False
        self.config.fit_sed = False

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ### SET PATHS

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = os.path.basename(path)

        # Get the full path to the 'data', 'prep' and 'in' directories
        self.data_path = os.path.join(path, self.config.data_dir)
        self.prep_path = os.path.join(path, self.config.prep_dir)
        self.in_path = os.path.join(path, self.config.in_dir)
        self.config_path = os.path.join(path, self.config.config_dir)
        self.extra_path = os.path.join(path, self.config.extra_dir)

        # Create the preparation and input directories if they were not yet present
        try: os.mkdir(self.prep_path)
        except OSError: pass
        try: os.mkdir(self.in_path)
        except OSError: pass

    # *****************************************************************

    def run(self):

        """
        This function runs the image preparation procedure
        :return:
        """

        # 1. Prepare
        if self.config.prepare: self.prepare_images()

        # 2. Fit bulge and disk
        if self.config.decompose: self.fit_bulge_and_disk()

        # 4. Make maps
        if self.config.make_maps: self.make_maps()

        # 5. Run SKIRT simulations, fit the SED
        if self.config.fit_sed: self.fit_sed()

    # *****************************************************************

    def prepare_images(self):

        """
        This function ...
        :return:
        """

        # Find the input FITS files
        self.find_input_files()

        # Loop over all filters for which we have an image
        for filter_name, path in self.image_paths.items():

            ### IF THE FINAL.FITS FILE EXISTS, SKIP THIS IMAGE

            final_path = os.path.join(self.prep_path, filter_name, "final.fits")
            if os.path.isfile(final_path): continue

            ### CONFIGURATION FOR THE PREPARATION

            # Determine the path to the default configuration file
            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
            default_config = os.path.join(directory, "config", "imagepreparation.cfg")

            # Look for a user configuration file
            config_path = os.path.join(self.config_path, filter_name + ".cfg")
            config = config_path if os.path.isfile(config_path) else None

            # Open the default configuration if no configuration file is specified, otherwise adjust the default
            # settings according to the user defined configuration file
            if config is None: config = configuration.open(default_config)
            else: config = configuration.open(config, default_config)

            # Set saving parameters for galaxy extractor
            config.galaxy_extraction.save_region = True
            config.galaxy_extraction.save_masked_frame = True
            config.galaxy_extraction.save_result = True
            config.galaxy_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "galaxy.reg")
            config.galaxy_extraction.saving.region_annotation = "name"
            config.galaxy_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_galaxies.fits")
            config.galaxy_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedgalaxies.fits")

            # Set saving parameters for star extractor
            config.star_extraction.save_region = True
            config.star_extraction.save_masked_frame = True
            config.star_extraction.save_result = True
            config.star_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "stars.reg")
            config.star_extraction.saving.region_annotation = "flux"
            config.star_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_stars.fits")
            config.star_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedstars.fits")

            ### OPENING THE IMAGE

            # Open the image
            image = Image(path)

            # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
            if image.frames.errors is None:

                error_path = os.path.join(self.data_path, filter_name + "_error.fits")
                if os.path.isfile(error_path): image.load_frames(error_path, 0, config.errors, "the error map")

            if image.frames.errors is None: log.warning("No error data found for " + filter_name)

            ### SELECTING THE APPROPRIATE FRAMES

            # Select the primary and errors frame
            image.deselect_all()
            image.frames[config.primary].select()
            if config.errors in image.frames: image.frames[config.errors].select()

            ### SETTING FLAGS

            # Set the extract_stars flag
            if image.wavelength < 10.0 * u.micron: config.extract_stars = True
            else: config.extract_stars = False

            ### SETTING THE UNITS

            # Set the unit
            unit = u.Unit(config.unit)
            image.set_unit(unit)

            ### SETTING THE FWHM

            # Set the FWHM of the PSF
            fwhm = config.fwhm * u.Unit(config.fwhm_unit) if config.fwhm is not None else None
            image.frames[config.primary].set_fwhm(fwhm)

            ### REMOVING NANS AND BAD REGIONS

            # Replace nans by zeros
            image.frames[config.primary].replace_nans(0.0)

            # Replace pixels in the 'extra' region with zeros
            extra_path = os.path.join(self.extra_path, image.name + '.reg')
            if os.path.isfile(extra_path):

                image.import_region(extra_path, "extra")
                image.regions.extra.select()
                image.create_mask()
                image.regions.extra.deselect()
                image.masks.extra.select()
                image.apply_masks(0.0)

            ### SETTING THE REFERENCE IMAGE

            # Set the path to the reference image for the rebinning
            config.rebinning.rebin_to = os.path.join(self.data_path, self.config.reference_image)

            ### PERFORMING THE PREPARATION

            # Create the preparation object
            preparation = ImagePreparation(config)

            # Run the preparation
            preparation.run(image)

            ### SAVING

            # Save the result
            image.save(final_path)

    # *****************************************************************

    def prepare_images_old(self):

        """
        This function prepares the images
        """

        config = self.config.preparation

        # Loop over all filters for which we have an image
        for filter_name, path in self.image_paths.items():

            # Path where the prepared images are being saved to
            filter_prep_path = os.path.join(self.prep_path, filter_name)

            # Open the image
            image = iu.open(path)

            # Set the unit
            unit = u.Unit(units[filter_name])
            image.primary.set_unit(unit)

            # Set the fwhm of the image, if it is not None
            if fwhmax[filter_name] is not None: image.primary.set_fwhm(fwhm[filter_name])

            # Mask NaNs, edges and extra user-defined regions
            extra_path = os.path.join(self.data_path, 'extra', image.name + '.reg')
            extra_path = extra_path if os.path.isfile(extra_path) else None

            # Mask nans and extra regions
            image.frames.primary.replace_nans(0.0)
            if extra_path is not None:

                image.import_region(extra_path, "extra")
                image.regions.extra.select()
                image.create_mask()
                image.regions.extra.deselect()
                image.masks.extra.select()
                image.apply_masks(0.0)

            # Extract galaxies from the image
            galaxyex = GalaxyExtractor(config.galaxyextraction)
            galaxyex.run(image.frames.primary)

            # Extract the stars from the image
            if image.wavelength < 10.0 * u.micron:

                starex = StarExtractor(config.starextraction)
                starex.run(image.frames.primary, galaxyex)

            else: starex = None

            # Subtract the sky
            assert (sky_subtracted[filter_name] == image.sky_subtracted)

            # Extract the sky from the image
            if not image.sky_subtracted:

                skyex = SkyExtractor()
                skyex.run(image.frames.primary, galaxyex, starex)

            # Determine whether galactic extinction should be taken into account
            # If the wavelength is smaller than 1 micron, such extinction is expected
            extinction = image.wavelength < 1.0 * u.micron
            if extinction: image.frames.primary *= 10**(0.4*attenuations[filter_name])

            exit()

            # Convert the units to MJy / sr
            unit = u.Unit("MJy/sr")
            image.primary.convert_to(unit)

            # Convolution
            if self.config.convolve and filter_name != self.config.convolve_to:

                reference_path = "Kernel_HiRes_" + aniano_names[filter_name] + "_to_" + aniano_names[self.config.convolve_to] + ".fits"
                reference_image = Image(reference_path)
                reference_frame = reference_image.frames.primary

                # Convolve the primary and errors frame (if present)
                image.frames.primary.convolve(reference_frame)
                if image.frames.errors is not None: image.frames.errors.convolve(reference_frame)

            # Rebinning
            if self.config.rebin and filter_name != self.config.rebin_to:

                reference_path = self.image_paths[self.config.rebin_to]
                reference_image = Image(reference_path)
                reference_frame = reference_image.frames.primary

                # Rebin the primary and errors frame (if present)
                image.frames.primary.rebin(reference_frame)
                if image.frames.errors is not None: image.frames.errors.rebin(reference_frame)

            # Set the uncertainties
            if filter_name != "2MASSH": iu.set_uncertainty(image, self.prep_path, "noise.reg")

            # If requested, save the result
            if self.config.save: iu.save(image, filter_prep_path, 'convolved_rebinned.fits')

            # Crop the interesting part of the image
            image.frames.primary.crop(350, 725, 300, 825)
            if image.frames.errors is not None: image.frames.errors.crop(350, 725, 300, 825)

            # If requested, save the result
            iu.save(image, filter_prep_path, 'final.fits')

    # *****************************************************************

    def fit_bulge_and_disk(self):

        """
        This function ...
        :return:
        """

        # TODO: do the bulge/disk fitting here

        # Set path of bulge and disk images
        bulge_path = os.path.join(self.prep_path, "Bulge", "M81_bulge_i59_total.fits")
        disk_path = os.path.join(self.prep_path, "Disk", "M81_disk_i59_total.fits")

        # Create list
        paths = {"Bulge": bulge_path, "Disk": disk_path}

        # For bulge and disk ...
        for name, path in paths.items():

            # Open the image
            image = iu.open(path)

            # Set the header of the image
            image.header["EQUINOX"] = 2000.0
            image.header["NAXIS"] = 2
            image.header["NAXIS1"] = 1000
            image.header["NAXIS2"] = 1000
            image.header["CRPIX1"] = 500.5
            image.header["CRPIX2"] = 500.5
            image.header["CRVAL1"] = 148.8883333
            image.header["CRVAL2"] = +69.06527778
            image.header["CD1_1"] = -4.77942772e-4
            image.header["CD1_2"] = 0.0
            image.header["CD2_1"] = 0.0
            image.header["CD2_2"] = 4.77942772e-4
            image.header["CROTA1"] = 0.0
            image.header["CROTA2"] = 0.0
            image.header["CTYPE1"] = 'RA---TAN'
            image.header["CTYPE2"] = 'DEC--TAN'

            # Convolve the bulge image to the PACS 160 resolution
            iu.convolve(image, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")

            # Rebin the convolved image to the frame of the PACS 160 image (just as we did with the other images)
            iu.rebin(image, self.data_path, 'PACS160.fits')
            image.crop(350, 725, 300, 825)

            # Save the convolved, rebinned and cropped bulge or disk image
            iu.save(image, os.path.join(self.prep_path, name), 'final.fits')

    # *****************************************************************

    def make_maps(self):

        """
        This function makes the maps of dust and stars ...
        """

        # Make the specific star formation rate map
        self.make_ssfr_map()

        # Make the FUV attenuation map
        self.make_attenuation_map()

        # Make the old stars map
        self.make_oldstars_map()

        # Make FUV and MIPS 24 emission maps
        self.make_fuv_and_24_maps()

        # Make ionizing stars map
        self.make_ionizing_stars_map()

    # *****************************************************************

    def make_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the prepared FUV and H images
        fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        h_path = os.path.join(self.prep_path, "2MASSH", "final.fits")

        # Open the FUV image and import the H image
        fuv_image = iu.open(fuv_path)
        fuv_image.import_datacube(h_path, "h")

        # Young non-ionizing stars (specific star formation rate) = GALEXFUV - H
        fuv_h_data = -2.5*(np.log10(fuv_image.frames.primary) - np.log10(fuv_image.frames.h))

        # Add the sSFR map to the FUV image
        fuv_image.add_frame(fuv_h_data, "ssfr")

        # Deselect all masks, regions and frames
        fuv_image.deselect_all()

        # Select the sSFR map
        fuv_image.frames.ssfr.select()

        # Mask nans in the sSFR map
        fuv_image.mask_nans()
        fuv_image.masks.nans.select()
        fuv_image.apply_masks(0.0)

        # Deselect all masks, regions and frames except for the sSFR frame
        iu.reset_selection(fuv_image)
        fuv_image.frames.primary.deselect()
        fuv_image.frames.ssfr.select()

        # Set the sSFR to zero in low signal-to-noise pixels
        fuv_image.frames.ssfr[(fuv_image.frames.h < 0.0) + (fuv_image.frames.primary < 10.0*fuv_image.frames.errors)] = 0.0

        # Save the sSFR map as FITS file
        iu.save(fuv_image, self.in_path, "ssfr.fits")

    # *****************************************************************

    def make_tir_map(self):

        """
        This function ...
        :return:
        """

        D_M81 = 3.6

        factor24 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/24e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor70 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/70e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor160 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/160e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        mips24_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        pacs70_path = os.path.join(self.prep_path, "PACS70", "final.fits")
        pacs160_path = os.path.join(self.prep_path, "PACS160", "final.fits")

        # Open the images
        mips24_image = iu.open(mips24_path)
        pacs70_image = iu.open(pacs70_path)
        pacs160_image = iu.open(pacs160_path)

        # Convert from MJy/sr to L_sun
        mips24_image.multiply(10.0**factor24)
        pacs70_image.multiply(10.0**factor70)
        pacs160_image.frames.errors.select()
        pacs160_image.multiply(10.0**factor160)

        # Galametz+2013 formula for Lsun units
        tir_data = (2.133*mips24_image.frames.primary) + \
                   (0.681*pacs70_image.frames.primary) + \
                   (1.125*pacs160_image.frames.primary)

        # Convert TIR from Lsun to W/m2
        factor = np.log10(3.846e26) - np.log10(4*np.pi) - (2*np.log10(D_M81*3.08567758e22))
        tir_data *= 10.0**factor

        # Return the TIR map
        return tir_data

    # *****************************************************************

    def make_attenuation_map(self):

        """
        This function ...
        :return:
        """

        # Dust = FUV attenuation = ratio of TIR and FUV luminosity

        fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")

        fuv_image = iu.open(fuv_path)

        # Convert the FUV map
        factor = - 20 + np.log10(3e8) - np.log10(0.153e-6) + (2*np.log10(2.85/206264.806247))
        fuv_converted_data = fuv_image.frames.primary * 10**factor

        tir_data = self.make_tir_map()

        tir_fuv_ratio_data = np.log10(tir_data/fuv_converted_data)

        x = tir_fuv_ratio_data
        x2 = np.power(tir_fuv_ratio_data, 2.0)
        x3 = np.power(tir_fuv_ratio_data, 3.0)
        x4 = np.power(tir_fuv_ratio_data, 4.0)

        #a_fuv_buat = (-0.0333*x3) + (0.3522*x2) + (1.1960*x) + 0.4967

        # Create an empty image
        a_fuv_cortese = np.zeros_like(tir_data)

        # Open the sSFR image
        ssfr_path = os.path.join(self.in_path, "ssfr.fits")
        ssfr_image = iu.open(ssfr_path)
        ssfr_data = ssfr_image.frames.primary

        #print (ssfr_data < 4.0).shape
        #print a_fuv_cortese.shape
        #print a_fuv_cortese[ssfr_data < 4.0]

        cnd = ssfr_data < 4.0
        a_fuv_cortese[cnd] = 0.50994 + (0.88311*x[cnd]) + (0.53315*x2[cnd]) + (0.04004*x3[cnd]) - (0.04883*x4[cnd])
        
        cnd = (ssfr_data >= 4.0) * (ssfr_data < 4.2)
        a_fuv_cortese[cnd] = 0.49867 + (0.86377*x[cnd]) + (0.51952*x2[cnd]) + (0.04038*x3[cnd]) - (0.04624*x4[cnd])
        
        cnd = (ssfr_data >= 4.2) * (ssfr_data < 4.6)
        a_fuv_cortese[cnd] = 0.49167 + (0.85201*x[cnd]) + (0.51152*x2[cnd]) + (0.04060*x3[cnd]) - (0.04475*x4[cnd])
        
        cnd = (ssfr_data >= 4.6) * (ssfr_data < 5.0)
        a_fuv_cortese[cnd] = 0.48223 + (0.83642*x[cnd]) + (0.50127*x2[cnd]) + (0.04092*x3[cnd]) - (0.04288*x4[cnd])
        
        cnd = (ssfr_data >= 5.0) * (ssfr_data < 5.4)
        a_fuv_cortese[cnd] = 0.46909 + (0.81520*x[cnd]) + (0.48787*x2[cnd]) + (0.04138*x3[cnd]) - (0.04050*x4[cnd])
        
        cnd = (ssfr_data >= 5.4) * (ssfr_data < 5.8)
        a_fuv_cortese[cnd] = 0.45013 + (0.78536*x[cnd]) + (0.47009*x2[cnd]) + (0.04210*x3[cnd]) - (0.03745*x4[cnd])
        
        cnd = (ssfr_data >= 5.8) * (ssfr_data < 6.3)
        a_fuv_cortese[cnd] = 0.42168 + (0.74191*x[cnd]) + (0.44624*x2[cnd]) + (0.04332*x3[cnd]) - (0.03362*x4[cnd])
        
        cnd = (ssfr_data >= 6.3) * (ssfr_data < 6.6)
        a_fuv_cortese[cnd] = 0.40210 + (0.71272*x[cnd]) + (0.43139*x2[cnd]) + (0.04426*x3[cnd]) - (0.03140*x4[cnd])
        
        cnd = (ssfr_data >= 6.6) * (ssfr_data < 6.9)
        a_fuv_cortese[cnd] = 0.37760 + (0.67674*x[cnd]) + (0.41420*x2[cnd]) + (0.04555*x3[cnd]) - (0.02900*x4[cnd])
        
        cnd = (ssfr_data >= 6.9) * (ssfr_data < 7.2)
        a_fuv_cortese[cnd] = 0.34695 + (0.63224*x[cnd]) + (0.39438*x2[cnd]) + (0.04739*x3[cnd]) - (0.02650*x4[cnd])
        
        cnd = (ssfr_data >= 7.2) * (ssfr_data < 7.5)
        a_fuv_cortese[cnd] = 0.30899 + (0.57732*x[cnd]) + (0.37157*x2[cnd]) + (0.05000*x3[cnd]) - (0.02399*x4[cnd])
        
        cnd = (ssfr_data >= 7.5) * (ssfr_data < 7.8)
        a_fuv_cortese[cnd] = 0.26302 + (0.51013*x[cnd]) + (0.34522*x2[cnd]) + (0.05377*x3[cnd]) - (0.02164*x4[cnd])
        
        cnd = (ssfr_data >= 7.8) * (ssfr_data < 8.1)
        a_fuv_cortese[cnd] = 0.20982 + (0.42980*x[cnd]) + (0.31431*x2[cnd]) + (0.05909*x3[cnd]) - (0.01957*x4[cnd])
        
        cnd = (ssfr_data >= 8.1) * (ssfr_data < 8.4)
        a_fuv_cortese[cnd] = 0.15293 + (0.33799*x[cnd]) + (0.27713*x2[cnd]) + (0.06638*x3[cnd]) - (0.01792*x4[cnd])
        
        cnd = (ssfr_data >= 8.4) * (ssfr_data < 8.8)
        a_fuv_cortese[cnd] = 0.09944 + (0.24160*x[cnd]) + (0.23161*x2[cnd]) + (0.07580*x3[cnd]) - (0.01671*x4[cnd])
        
        cnd = (ssfr_data >= 8.8) * (ssfr_data < 9.2)
        a_fuv_cortese[cnd] = 0.05822 + (0.15524*x[cnd]) + (0.17801*x2[cnd]) + (0.08664*x3[cnd]) - (0.01593*x4[cnd])
        
        cnd = (ssfr_data >= 9.2) * (ssfr_data < 9.6)
        a_fuv_cortese[cnd] = 0.03404 + (0.09645*x[cnd]) + (0.12452*x2[cnd]) + (0.09679*x3[cnd]) - (0.01548*x4[cnd])
        
        cnd = (ssfr_data >= 9.6) * (ssfr_data < 10.0)
        a_fuv_cortese[cnd] = 0.02355 + (0.06934*x[cnd]) + (0.08725*x2[cnd]) + (0.10339*x3[cnd]) - (0.01526*x4[cnd])
        
        cnd = (ssfr_data >= 9.6) * (ssfr_data < 10.0)
        a_fuv_cortese[cnd] = 0.02025 + (0.06107*x[cnd]) + (0.07212*x2[cnd]) + (0.10588*x3[cnd]) - (0.01517*x4[cnd])

        # Open the PACS160 image
        pacs160_path = os.path.join(self.prep_path, "PACS160", "final.fits")
        pacs160_image = iu.open(pacs160_path)

        # Set pixels to zero in some cases
        a_fuv_cortese[(ssfr_data < 0) + (ssfr_data > 10.5) + (fuv_image.frames.primary <= 0) + (pacs160_image.frames.primary < 5.0*pacs160_image.frames.errors)] = 0.0

        # Make sure all pixels are larger or equal to zero
        a_fuv_cortese[a_fuv_cortese < 0.0] = 0.0

        attenuation_image = iu.new("attenuation")
        attenuation_image.add_frame(a_fuv_cortese, "primary")
        attenuation_image.frames.primary.select()

        # Save FUV attenuation map as FITS file
        iu.save(attenuation_image, self.in_path, "fuv_attenuation.fits")

    # *****************************************************************

    def make_oldstars_map(self):

        """
        This function ...
        :return:
        """

        # Old stars = IRAC3.6 - bulge
        # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

        irac_path = os.path.join(self.prep_path, "IRACI1", "final.fits")
        irac_image = iu.open(irac_path)

        #bulge_path = os.path.join(self.prep_path, "Bulge", "M81_bulge_i59_total.fits")
        bulge_path = os.path.join(self.prep_path, "Bulge", "final.fits")
        irac_image.import_datacube(bulge_path, "bulge")

        irac_image.frames.bulge.select()

        totalbulge = irac_image.frames.bulge.sum
        flux3_6 = 60541.038
        factor = 0.54 * flux3_6 / totalbulge

        # Multiply the bulge frame with the factor
        #irac_image.multiply(factor)

        # Do the subtraction
        #irac_image.subtract()

        oldstars_data = irac_image.frames.primary - factor*irac_image.frames.bulge

        oldstars_data[(irac_image.frames.primary < 10.0*irac_image.frames.errors) + (oldstars_data < 0.0)] = 0.0

        # Add old stars frame
        irac_image.add_frame(oldstars_data, "oldstars")

        irac_image.deselect_all()

        irac_image.frames.oldstars.select()

        # Save old stars map as FITS file
        iu.save(irac_image, self.in_path, "old_stars.fits")

    # *****************************************************************

    def make_fuv_and_24_maps(self):

        """
        This function ...
        :return:
        """

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        #disk_path = os.path.join(self.prep_path, "Disk", "M81_disk_i59_total.fits")
        disk_path = os.path.join(self.prep_path, "Disk", "final.fits")

        disk_image = iu.open(disk_path)

        totaldisk = disk_image.frames.primary.sum

        fluxFUV = 855.503
        flux24 = 27790.448

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factorFUV = 0.2*fluxFUV/totaldisk
        factor24 = 0.48*flux24/totaldisk

        fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        fuv_image = iu.open(fuv_path)

        mips_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        mips_image = iu.open(mips_path)

        new_fuv_data = fuv_image.frames.primary - factorFUV*disk_image.frames.primary
        new_mips_data = mips_image.frames.primary - factor24*disk_image.frames.primary

        # Set zero where negative
        new_fuv_data[new_fuv_data < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_fuv_data[fuv_image.frames.primary < 10.0*fuv_image.frames.errors] = 0.0

        # Set zero where negative
        new_mips_data[new_mips_data < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_mips_data[mips_image.frames.primary < 10.0*mips_image.frames.errors] = 0.0

        new_fuv_image = iu.new("fuv")
        new_mips_image = iu.new("mips")

        new_fuv_image.add_frame(new_fuv_data, "primary")
        new_mips_image.add_frame(new_mips_data, "primary")
        new_mips_image.add_frame(mips_image.frames.errors, "errors")

        new_fuv_image.frames.primary.select()
        new_mips_image.frames.primary.select()
        new_mips_image.frames.errors.select()

        # Save the new images
        iu.save(new_fuv_image, self.in_path, "fuv.fits")
        iu.save(new_mips_image, self.in_path, "mips.fits")

    # *****************************************************************

    def make_ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        #Young ionizing stars = Ha + 0.031 x MIPS24

        D_M81 = 3.6

        # Convert to Lsun
        factor_ha = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/0.657894736e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        ha_path = os.path.join(self.prep_path, "Ha", "final.fits")
        ha_image = iu.open(ha_path)

        ha_image.multiply(10**factor_ha)

        new_mips_path = os.path.join(self.in_path, "mips.fits")
        new_mips_image = iu.open(new_mips_path)

        ionizing_stars_data = ha_image.frames.primary + 0.031*new_mips_image.frames.primary
        ionizing_stars_ratio = ha_image.frames.primary / (0.031*new_mips_image.frames.primary)

        low_snr = (ha_image.frames.primary < 10.0*ha_image.frames.errors) + (new_mips_image.frames.primary < 10.0*new_mips_image.frames.errors)

        ionizing_stars_data[low_snr] = 0.0
        ionizing_stars_ratio[low_snr] = 0.0

        ionizing_stars_data[ionizing_stars_data < 0.0] = 0.0
        ionizing_stars_ratio[ionizing_stars_data < 0.0] = 0.0

        ionizing_stars_image = iu.new("ionizingstars")

        ionizing_stars_image.add_frame(ionizing_stars_data, "primary")
        ionizing_stars_image.add_frame(ionizing_stars_ratio, "ratio")

        ionizing_stars_image.frames.primary.select()
        ionizing_stars_image.frames.ratio.select()

        # Save the new image
        iu.save(ionizing_stars_image, self.in_path, "ionizingstars.fits")

    # *****************************************************************
    
    def fit_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # *****************************************************************

    def find_input_files(self):

        """
        This function ...
        :return:
        """

        # Get a list of files in the data directory
        files = [f for f in os.listdir(self.data_path) if os.path.isfile(os.path.join(self.data_path,f))]

        # Create a dictionary holding the path of each valid FITS file with a key that represents the filter
        self.image_paths = dict()

        # Loop over all files in the data directory
        for filename in files:

            # Ignore non-FITS files or hidden files
            if not filename.endswith(".fits") or filename.startswith("."): continue

            # Ignore error maps
            if "error" in filename: continue

            # Get the name of the file without the extension
            base_filename = os.path.splitext(filename)[0]

            # If a filtername was specified, only add the file that corresponds to this filter
            if filter_name is not None:

                if filter_name.lower() == base_filename.lower():

                    self.image_paths[filter_name] = os.path.join(self.data_path, filename)
                    break

            # If no filtername was specified, add each FITS file found in the data directory to the dictionary
            else: self.image_paths[base_filename] = os.path.join(self.data_path, filename)

            # If intermediate results should be saved, create a seperate directory for each filter
            if self.config.save:

                # Create the directory if it was not yet present
                try: os.mkdir(os.path.join(self.prep_path, base_filename))
                except OSError: pass

# *****************************************************************



