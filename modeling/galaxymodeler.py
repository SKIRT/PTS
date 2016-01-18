#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.galaxymodeler Model a galaxy by fitting using Astromagic and SKIRT
#
# An instance of the GalaxyModeler class in this module is responsible for taking reduced astronomical image data
# of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
# dust, star formation and old stars. Then, it uses these maps as input for SKIRT radiative transfer simulations
# and fits the output to the observed galaxy SED.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import astronomical modules
from astropy import units as u

# Import the relevant AstroMagic classes and modules
from ..magic.core import Image, Frame

# Import the relevant PTS classes and modules
from .core import ImagePreparation, MapMaker, SEDFitter, PhotoMeter
from ..core.basics.configurable import Configurable
from ..core.tools import filesystem, tables, time

# -----------------------------------------------------------------

class GalaxyModeler(Configurable):

    """
    An instance of the GalaxyModeler class in this module is responsible for taking reduced astronomical image data
    of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
    dust, star formation and old stars. Then, it uses these maps as input for SKIRT radiative transfer simulations
    and fits the output to the observed galaxy SED.
    """

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param directory:
        :param filter_name:
        :param plot:
        :param save:
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(GalaxyModeler, self).__init__(config, "modeling")

        # -- Temporary --

        self.config.decompose = False
        self.config.make_maps = False
        self.config.fit_sed = False

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Modeler instance
        modeler = cls(arguments.config)

        # Logging
        if arguments.debug:

            modeler.config.logging.level = "DEBUG"
            modeler.config.logging.cascade = True

        modeler.config.write_steps = arguments.steps

        if arguments.report: modeler.config.logging.path = time.unique_name("log") + ".txt"

        # Set the input and output path
        #modeler.config.input_path = arguments.input_path
        #modeler.config.output_path = arguments.output_path

        # Check if a specific stage is defined
        if arguments.stage is not None:

            modeler.config.prepare = False
            modeler.config.decompose = False
            modeler.config.make_maps = False
            modeler.config.fit_sed = False

            if arguments.stage == "preparation": modeler.config.prepare = True
            elif arguments.stage == "photometry": modeler.config.do_photometry = True
            elif arguments.stage == "mapmaking": modeler.config.make_maps = True
            elif arguments.stage == "fitting": modeler.config.fit_sed = True
            else: raise ValueError("Unkown stage (choose 'preparation', 'decomposition', 'mapmaking' or 'fitting')")

        # Return the new instance
        return modeler

    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function runs the image preparation procedure
        :return:
        """

        # 1. Call the setup function
        self.setup(path)

        # 2. Prepare the images
        if self.config.prepare: self.prepare_images()

        # 3. Perform the photometry
        if self.config.do_photometry: self.do_photometry()

        # 3. Fit bulge and disk
        if self.config.decompose: self.decompose()

        # 4. Make maps
        if self.config.make_maps: self.make_maps()

        # 5. Run SKIRT simulations, fit the SED
        if self.config.fit_sed: self.fit_sed()

    # -----------------------------------------------------------------

    def setup(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # -- Children --

        # Create the preparation object
        self.add_child("image_preparer", ImagePreparation, self.config.preparation)
        self.add_child("photometer", PhotoMeter)
        self.add_child("map_maker", MapMaker)
        self.add_child("sed_fitter", SEDFitter)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(GalaxyModeler, self).setup()

        # -- Attributes --

        self.image_preparer.config.write_steps = self.config.write_steps

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

    # -----------------------------------------------------------------

    def prepare_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Preparing the images ...")

        # Load image information
        info_path = os.path.join(self.data_path, "info.dat")
        info_table = tables.from_file(info_path)

        # Loop over all files found in the data directory
        for file_path in filesystem.files_in_path(self.data_path, extension="fits", not_contains="error"):

            # Get the file name
            file_name = os.path.basename(file_path)
            image_name = os.path.splitext(file_name)[0]

            # Inform the user
            self.log.debug("Preparing the " + image_name + " image")

            image_output_path = os.path.join(self.prep_path, image_name)

            # Get the corresponding index in the information table
            info_index = tables.find_index(info_table, image_name)
            # TEMP: skip image if not defined in table !!
            if info_index is None: continue

            # Check whether this image already has a prepared image
            final_path = os.path.join(image_output_path, "result.fits")
            if filesystem.is_file(final_path): continue

            # Open the image
            image = Image(file_path)
            self.load_error_frame_and_select(image, image_name)

            # -- Units --

            # Set the unit
            unit = u.Unit(info_table["Unit"][info_index])
            image.set_unit(unit)

            # -- The FWHM of the PSF --

            # Set the FWHM of the PSF
            fwhm = info_table["FWHM"][info_index] * u.Unit(info_table["Unit of FWHM"][info_index]) if not info_table["FWHM"].mask[info_index] else None
            image.frames["primary"].set_fwhm(fwhm)

            # -- Extra region mask --

            extra_path = os.path.join(self.data_path, "extra", image_name + ".reg")
            # Check if a region is present
            # Check whether the extra region file is present
            if filesystem.is_file(extra_path): self.image_preparer.config.extra_region_path = extra_path
            else: self.image_preparer.config.extra_region_path = None

            # -- The reference image (for convolution and rebinning)

            # Set the path to the reference image for the rebinning
            reference_path = os.path.join(self.data_path, self.config.reference_image + ".fits")

            if file_path == reference_path:

                self.image_preparer.config.rebin = False
                self.image_preparer.config.convolve = False

            else:

                self.image_preparer.config.rebin = True
                self.image_preparer.config.convolve = True

                self.image_preparer.config.rebin_to = reference_path
                self.image_preparer.config.convolve_to = self.config.reference_image.replace(' ', '')

            # -- Noise --

            # Set the path to the noise region
            noise_path = os.path.join(self.data_path, "noise.reg")
            self.image_preparer.config.uncertainties.noise_path = noise_path

            # -- Output --

            # Create the output directory if it does not exist for this image
            filesystem.create_directory(image_output_path)
            self.image_preparer.config.output_path = image_output_path

            # -- Result --

            # Save the result
            self.image_preparer.config.writing.result_path = "result.fits"

            # Run the image preparation
            self.image_preparer.run(image)

    # -----------------------------------------------------------------

    def calculate_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Calculating the photometry for the different bands ...")

        # Loop over all directories found in the prep directory
        for directory_path in filesystem.directories_in_path(self.prep_path):

            name = os.path.basename(directory_path)

            # Debug info
            self.log.info("Calculating the photometry for " + name +  " ...")

            path = os.path.join(directory_path, "result.fits")

            if not filesystem.is_file(path): raise ValueError("No resulting image present for " + name)

            # Check whether a result FITS file is present
            frame = Frame.from_file(path)

            self.photometer.run(frame, None)

    # -----------------------------------------------------------------

    def decompose(self):

        """
        This function ...
        :return:
        """

        # Create a GalaxyDecomposer object
        #decomposer = GalaxyDecomposer()

        # Set log level for the decomposer, if cascading is enabled
        #if self.config.logging.cascade:

        #    decomposer.config.logging.level = self.config.logging.level
        #    decomposer.config.logging.cascade = self.config.logging.cascade
        #    decomposer.config.logging.path = self.config.logging.path

        # TODO: do the bulge/disk fitting here
        # Run the decomposition
        #decomposer.run()

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

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function makes the maps of dust and stars ...
        """

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

        # Run the map maker
        self.map_maker.run()

    # -----------------------------------------------------------------
    
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

        # Run the SED fitting
        self.sed_fitter.run()

    # -----------------------------------------------------------------

    def load_error_frame_and_select(self, image, image_name):

        """
        This function ...
        :param image:
        :param name:
        :return:
        """

        # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
        if image.frames.errors is None:

            error_path = os.path.join(self.data_path, image_name + " error.fits")
            if os.path.isfile(error_path): image.load_frames(error_path, 0, "errors", "the error map")

        # Still no errors frame
        if image.frames.errors is None: self.log.warning("No error data found for " + image_name + ".fits")

        image.deselect_all()
        image.frames["primary"].select()
        if "errors" in image.frames: image.frames["errors"].select()

    # -----------------------------------------------------------------

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

# -----------------------------------------------------------------
