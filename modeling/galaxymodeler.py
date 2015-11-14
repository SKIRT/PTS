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

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os.path
import inspect
import shutil

# Import astronomical modules
from astropy import units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic.tools import configuration
from astromagic import Image
from astromagic.core.frames import Frame

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
        self.manual_path = os.path.join(path, self.config.manual_dir)

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

    def clear_output(self):

        """
        This function ...
        :return:
        """

        shutil.rmtree('/home/me/test')

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

            # Set log level for the different children of the ImagePreparation object, if cascading is enabled
            if self.config.logging.cascade:

                # Galaxy extractor
                config.galaxy_extraction.logging.level = self.config.logging.level
                config.galaxy_extraction.logging.cascade = self.config.logging.cascade
                config.galaxy_extraction.logging.path = self.config.logging.path

                # Star extractor
                config.star_extraction.logging.level = self.config.logging.level
                config.star_extraction.logging.cascade = self.config.logging.cascade
                config.star_extraction.logging.path = self.config.logging.path

                # Sky extractor
                config.sky_extraction.logging.level = self.config.logging.level
                config.sky_extraction.logging.cascade = self.config.logging.cascade
                config.sky_extraction.logging.path = self.config.logging.path

            # Set saving parameters for galaxy extractor
            config.galaxy_extraction.save_table = True
            config.galaxy_extraction.save_region = True
            #config.galaxy_extraction.save_masked_frame = True
            #config.galaxy_extraction.save_result = True
            config.galaxy_extraction.saving.table_path = os.path.join(self.prep_path, filter_name, "galaxies.txt")
            config.galaxy_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "galaxies.reg")
            config.galaxy_extraction.saving.region_annotation = "name"
            #config.galaxy_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_galaxies.fits")
            #config.galaxy_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedgalaxies.fits")

            # Set saving parameters for star extractor
            config.star_extraction.save_table = True
            config.star_extraction.save_region = True
            #config.star_extraction.save_masked_frame = True
            config.star_extraction.save_result = True
            config.star_extraction.saving.table_path = os.path.join(self.prep_path, filter_name, "stars.txt")
            config.star_extraction.saving.region_path = os.path.join(self.prep_path, filter_name, "stars.reg")
            config.star_extraction.saving.region_annotation = "flux"
            #config.star_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "masked_stars.fits")
            config.star_extraction.saving.result_path = os.path.join(self.prep_path, filter_name, "extractedstars.fits")

            # Set saving parameters for sky extractor
            config.sky_extraction.save_masked_frame = True
            config.sky_extraction.save_clipped_masked_frame = True
            config.sky_extraction.save_histogram = False
            config.sky_extraction.saving.masked_frame_path = os.path.join(self.prep_path, filter_name, "sky_mask.fits")
            config.sky_extraction.saving.clipped_masked_frame_path = os.path.join(self.prep_path, filter_name, "clipped_sky_mask.fits")
            config.sky_extraction.saving.histogram_path = os.path.join(self.prep_path, filter_name, "sky_histogram.pdf")

            # Temporary: do not include calibration errors
            config.uncertainties.add_calibration_error = False

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

            # Set the correct_for_extinction flag
            if image.wavelength < 1.0 * u.micron: config.correct_for_extinction = True
            else: config.correct_for_extinction = False

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

                # Set the path of the extra region in the SkyExtractor's configuration
                config.sky_extraction.extra_region = extra_path

            ### SETTING THE REFERENCE IMAGE

            # Set the path to the reference image for the rebinning
            config.rebinning.rebin_to = os.path.join(self.data_path, self.config.reference_image + ".fits")

            # Set the 'rebin' and 'convolve' flags
            if image.name == self.config.reference_image:

                log.info("This is the reference image, will not be rebinned or convolved")
                config.rebin = False
                config.convolve = False

            ### SETTING THE NOISE REGION

            # Set the path to the noise region
            noise_path = os.path.join(self.prep_path, "noise.reg")
            config.uncertainties.noise_path = noise_path

            ### LOOKING FOR REGIONS TO IGNORE FOR THE EXTRACTION ALGORITHMS

            # If a region file exists with the name of this image, add its path to the configuration for the preparation
            ignore_path = os.path.join(self.ignore_path, image.name + ".reg")
            if os.path.isfile(ignore_path): config.star_extraction.ignore_region = ignore_path

            ### LOOKING FOR MANUAL STAR REGIONS

            # If a region file exists with the name of this image, add its path to the configuration for the preparation
            manual_path = os.path.join(self.manual_path, image.name + ".reg")
            if os.path.isfile(manual_path): config.star_extraction.manual_region = manual_path

            ### PERFORMING THE PREPARATION

            # Create the preparation object
            preparation = ImagePreparation(config)

            # Run the preparation
            preparation.run(image)

            ### SAVING

            # Save the result
            image.save(final_path)

    # *****************************************************************

    def decompose(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "decomposer.cfg")

        # Create config for the GalaxyDecomposer object
        config = configuration.open(default_config)

        # Set log level for the decomposer, if cascading is enabled
        if self.config.logging.cascade:

            config.logging.level = self.config.logging.level
            config.logging.cascade = self.config.logging.cascade
            config.logging.path = self.config.logging.path

        # Create a GalaxyDecomposer object
        decomposer = GalaxyDecomposer(config)

        # Run the decomposition
        decomposer.run()

        # TODO: do the bulge/disk fitting here

        # Set path of bulge and disk images
        bulge_path = os.path.join(self.prep_path, "Bulge", "bulge.fits")
        disk_path = os.path.join(self.prep_path, "Disk", "disk.fits")

        # Create list
        paths = {"Bulge": bulge_path, "Disk": disk_path}

        # For bulge and disk ...
        for name, path in paths.items():

            # Open the frame
            frame = Frame.from_file(path)

            # Convolve the frame to the PACS 160 resolution
            kernels_dir = os.path.expanduser("~/Kernels")
            kernel_path = os.path.join(kernels_dir, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")
            kernel = Frame.from_file(kernel_path)
            frame.convolve(kernel)

            # Rebin the convolved frame to the PACS 160 frame (just as we did with the other images)
            reference_path = os.path.join(self.data_path, "PACS160.fits")
            reference = Frame.from_file(reference_path)
            frame.rebin(reference)

            # Finally, crop the image
            frame.frames.primary.crop(350, 725, 300, 825)

            # Save the convolved, rebinned and cropped bulge or disk frame
            component_path = os.path.join(self.prep_path, name, 'final.fits')
            frame.save(component_path)

    # *****************************************************************

    def make_maps(self):

        """
        This function makes the maps of dust and stars ...
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "mapmaker.cfg")

        # Create config for the MapMaker object
        config = configuration.open(default_config)

        ### CONFIGURATION FOR THE PREPARATION

        # Set log level for the map maker, if cascading is enabled
        if self.config.logging.cascade:

            config.logging.level = self.config.logging.level
            config.logging.cascade = self.config.logging.cascade
            config.logging.path = self.config.logging.path

        # Set the names of the primary and error frames
        config.primary = self.config.primary
        config.errors = self.config.errors

        # Open the prepared reference image
        config.cutoff.reference_path = os.path.join(self.prep_path, self.config.reference_image, "final.fits")

        # Set the path to the low signal-to-noise cutoff mask file
        config.saving.cutoff_mask_path = os.path.join(self.prep_path, self.config.reference_image, "cutoff_mask.fits")
        config.saving.cutoff_mask_segments_path = os.path.join(self.prep_path, self.config.reference_image, "cutoff_mask_segments.fits")
        config.saving.cutoff_mask_holes_path = os.path.join(self.prep_path, self.config.reference_image, "cutoff_mask_holes.fits")
        config.saving.cutoff_mask_with_holes_path = os.path.join(self.prep_path, self.config.reference_image, "cutoff_mask_with_holes.fits")

        # Set the paths to the processed images
        config.h_path = os.path.join(self.prep_path, "2MASSH", "final.fits")
        config.fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        config.ha_path = os.path.join(self.prep_path, "Ha", "final.fits")
        config.irac_path = os.path.join(self.prep_path, "IRACI1", "final.fits")
        config.mips_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        config.pacsblue_path = os.path.join(self.prep_path, "PACS70", "final.fits")
        config.pacsred_path = os.path.join(self.prep_path, "PACS160", "final.fits")
        config.disk_path = os.path.join(self.prep_path, "Disk", "final.fits")
        config.bulge_path = os.path.join(self.prep_path, "Bulge", "final.fits")

        # Set the paths to the cutoff maps
        config.saving.h_cutoff_path = os.path.join(self.prep_path, "2MASSH", "cutoff.fits")
        config.saving.fuv_cutoff_path = os.path.join(self.prep_path, "GALEXFUV", "cutoff.fits")
        config.saving.ha_cutoff_path = os.path.join(self.prep_path, "Ha", "cutoff.fits")
        config.saving.irac_cutoff_path = os.path.join(self.prep_path, "IRACI1", "cutoff.fits")
        config.saving.mips_cutoff_path = os.path.join(self.prep_path, "MIPS24", "cutoff.fits")
        config.saving.pacsblue_cutoff_path = os.path.join(self.prep_path, "PACS70", "cutoff.fits")
        config.saving.pacsred_cutoff_path = os.path.join(self.prep_path, "PACS160", "cutoff.fits")
        config.saving.disk_cutoff_path = os.path.join(self.prep_path, "Disk", "cutoff.fits")
        config.saving.bulge_cutoff_path = os.path.join(self.prep_path, "Bulge", "cutoff.fits")

        # Set the paths to the maps converted to solar luminosities
        config.conversion.ha_output_path = os.path.join(self.in_path, "solar", "ha.fits")
        config.conversion.ha_errors_output_path = os.path.join(self.in_path, "solar", "ha_errors.fits")
        config.conversion.mips_output_path = os.path.join(self.in_path, "solar", "mips.fits")
        config.conversion.mips_errors_output_path = os.path.join(self.in_path, "solar", "mips_errors.fits")
        config.conversion.pacsblue_output_path = os.path.join(self.in_path, "solar", "pacsblue.fits")
        config.conversion.pacsred_output_path = os.path.join(self.in_path, "solar", "pacsred.fits")

        # Set the paths to the output maps
        config.dust.output_path = os.path.join(self.in_path, "dust.fits")
        config.dust.ssfr.output_path = os.path.join(self.in_path, "ssfr.fits")  # Temporary ...
        config.dust.ssfr.color_output_path = os.path.join(self.in_path, "fuv_h_color.fits") # Temporary ...
        config.dust.ssfr.with_nans_output_path = os.path.join(self.in_path, "ssfr_withnans.fits")  # Temporary ...
        config.dust.tir_to_fuv_output_path = os.path.join(self.in_path, "tir_to_fuv.fits") # Temporary ...
        config.old_stars.output_path = os.path.join(self.in_path, "old_stars.fits")
        config.ionizing_stars.output_path = os.path.join(self.in_path, "ionizing_stars.fits")
        config.non_ionizing_stars.output_path = os.path.join(self.in_path, "non_ionizing_stars.fits")

        ### PERFORMING THE MAP MAKING

        # Create a MapMaker object
        maker = MapMaker(config)

        # Run the map maker
        maker.run()

    # *****************************************************************
    
    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Create config for SEDFitter ...
        config = None

        # Set log level for the SED fitter, if cascading is enabled
        if self.config.logging.cascade:

            config.logging.level = self.config.logging.level
            config.logging.cascade = self.config.logging.cascade
            config.logging.path = self.config.logging.path

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



