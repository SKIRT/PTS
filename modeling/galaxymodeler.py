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
from .imagepreparation import ImagePreparation
from .mapmaker import MapMaker
from .sedfitter import SEDFitter
from .galaxydecomposer import GalaxyDecomposer

# *****************************************************************

class GalaxyModeler(object):

    """
    An instance of the GalaxyModeler class in this module is responsible for taking reduced astronomical image data
    of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
    dust, star formation and old stars. Then, it uses these maps as input for SKIRT radiative transfer simulations
    and fits the output to the observed galaxy SED.
    """

    # *****************************************************************

    def __init__(self, path, config=None):

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
        self.ignore_path = os.path.join(path, self.config.ignore_dir)

    # *****************************************************************

    def run(self):

        """
        This function runs the image preparation procedure
        :return:
        """

        # 1. Prepare
        if self.config.prepare: self.prepare_images()

        # 2. Fit bulge and disk
        if self.config.decompose: self.decompose()

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
            #config.galaxy_extraction.save_masked_frame = True
            #config.galaxy_extraction.save_result = True
            config.galaxy_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "galaxies.reg")
            config.galaxy_extraction.saving.region_annotation = "name"
            #config.galaxy_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_galaxies.fits")
            #config.galaxy_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedgalaxies.fits")

            # Set saving parameters for star extractor
            config.star_extraction.save_region = True
            #config.star_extraction.save_masked_frame = True
            config.star_extraction.save_result = True
            config.star_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "stars.reg")
            config.star_extraction.saving.region_annotation = "flux"
            #config.star_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_stars.fits")
            config.star_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedstars.fits")

            # Set saving parameters for sky extractor
            config.sky_extraction.save_masked_frame = True
            config.sky_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "sky_mask.fits")

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

            ### LOOKING FOR REGIONS TO IGNORE FOR THE EXTRACTION ALGORITHMS

            # If a region file exists with the name of this image, add its path to the configuration for the preparation
            ignore_path = os.path.join(self.ignore_path, image.name + ".reg")
            if os.path.isfile(ignore_path): config.star_extraction.ignore_region = ignore_path

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

    def decompose(self):

        """
        This function ...
        :return:
        """

        # Create a GalaxyDecomposer object
        decomposer = GalaxyDecomposer()

        # Run the decomposition
        decomposer.run()



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

        # Create a MapMaker object
        maker = MapMaker()

        # Run the map maker
        maker.run()

    # *****************************************************************
    
    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Create a SEDFitter object
        fitter = SEDFitter()

        # Run the SED fitting
        fitter.run()

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
            if self.config.preparation.filter_name is not None:

                if self.config.preparation.filter_name.lower() == base_filename.lower():

                    self.image_paths[self.config.preparation.filter_name] = os.path.join(self.data_path, filename)
                    break

            # If no filtername was specified, add each FITS file found in the data directory to the dictionary
            else: self.image_paths[base_filename] = os.path.join(self.data_path, filename)

            # If intermediate results should be saved, create a seperate directory for each filter
            if self.config.save:

                # Create the directory if it was not yet present
                try: os.mkdir(os.path.join(self.prep_path, base_filename))
                except OSError: pass

# *****************************************************************



