#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.datapreparation Contains the DataPreparer class
#
# Info ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import astronomical modules
from astropy import units as u
from astroquery.irsa_dust import IrsaDust

# Import the relevant AstroMagic classes and modules
from ...magic import ImageImporter

# Import the relevant PTS classes and modules
from .component import ModelingComponent
from .imagepreparation import ImagePreparer
from ...core.tools import filesystem, tables, time

# -----------------------------------------------------------------

irsa_names = {"SDSS u": "SDSS u",
              "SDSS g": "SDSS g",
              "SDSS r": "SDSS r",
              "SDSS i": "SDSS i",
              "SDSS z": "SDSS z",
              "2MASS J": "2MASS J",
              "2MASS H": "2MASS H",
              "2MASS K": "2MASS Ks",
              "IRAC I1": "IRAC-1",
              "IRAC I2": "IRAC-2",
              "IRAC I3": "IRAC-3",
              "IRAC I4": "IRAC-4",
              "WISE W1": "WISE-1",
              "WISE W2": "WISE-2"}

# -----------------------------------------------------------------

aniano_names = {"GALEX FUV": "GALEX_FUV",
                "GALEX NUV": "GALEX_NUV",
                "SDSS u": "BiGauss_02.0",
                "SDSS g": "BiGauss_02.0",
                "SDSS r": "BiGauss_02.0",
                "Ha": "Gauss_03.0",
                "SDSS i": "BiGauss_02.0",
                "SDSS z": "BiGauss_02.0",
                "2MASS J": "Gauss_03.5",
                "2MASS H": "Gauss_03.0",
                "2MASS K": "Gauss_03.5",
                "WISE W1": "WISE_FRAME_3.4",
                "IRAC I1": "IRAC_3.6",
                "IRAC I2": "IRAC_4.5",
                "WISE W2": "WISE_FRAME_4.6",
                "IRAC I3": "IRAC_5.8",
                "IRAC I4": "IRAC_8.0",
                "WISE W3": "WISE_FRAME_11.6",
                "WISE W4": "WISE_FRAME_22.1",
                "MIPS 24": "MIPS_24",
                "PACS 70": "PACS_70",
                "PACS 160": "PACS_160",
                "SPIRE PSW": None,
                "SPIRE PLW": None}

# -----------------------------------------------------------------

calibration_errors = {"GALEX FUV": "0.05 mag",
                      "GALEX NUV": "0.03 mag",
                      "SDSS u": "2%",
                      "SDSS g": "2%",
                      "SDSS r": "2%",
                      "Ha": "5%",
                      "SDSS i": "2%",
                      "SDSS z": "2%",
                      "2MASS J": "0.03 mag",
                      "2MASS H": "0.03 mag",
                      "2MASS K": "0.03 mag",
                      "WISE W1": "2.4%",
                      "IRAC I1": "1.8%",
                      "IRAC I2": "1.9%",
                      "WISE W2": "2.8%",
                      "IRAC I3": "2.0%",
                      "IRAC I4": "2.1%",
                      "WISE W3": "4.5%",
                      "WISE W4": "5.7%",
                      "MIPS 24": "4%",
                      "PACS 70": "5%",
                      "PACS 160": "5%",
                      "SPIRE PSW": "4%",
                      "SPIRE PMW": "4%",
                      "SPIRE PLW": "4%"}

# -----------------------------------------------------------------

class DataPreparer(ModelingComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DataPreparer, self).__init__(config)

        # The list of images
        self.images = []

        # Information about the images
        self.attenuations = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Modeler instance
        preparer = cls(arguments.config)

        # Logging
        if arguments.debug:

            preparer.config.logging.level = "DEBUG"
            preparer.config.logging.cascade = True

        preparer.config.preparation.write_steps = arguments.steps

        # Set the reference image
        if arguments.reference is not None: preparer.reference_image = arguments.reference

        # Set the input and output path
        preparer.config.path = arguments.path
        preparer.config.input_path = os.path.join(arguments.path, "data")
        preparer.config.output_path = os.path.join(arguments.path, "prep")

        # A single image can be specified so the preparation is only run with that image
        preparer.config.single_image = arguments.image

        # Set logging path
        if arguments.report: preparer.config.logging.path = os.path.join(preparer.config.output_path, time.unique_name("log") + ".txt")

        # Return the new instance
        return preparer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function runs the data preparation ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the input data
        self.load_data()

        # 3. Get attenuations
        self.get_attenuations()

        # 2. Prepare the images
        self.prepare()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # -- Children --

        # Create the preparation object
        self.add_child("importer", ImageImporter, self.config.importation)
        self.add_child("image_preparer", ImagePreparer, self.config.preparation)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(DataPreparer, self).setup()

        # -- Fixed properties for the image preparer (valid for all target images)

        # Set the path to the reference image for the rebinning
        reference_path = os.path.join(self.data_path, self.config.reference_image + ".fits")

        # Set the path of the rebinning reference path and the kernel image
        self.image_preparer.config.rebinning.rebin_to = reference_path

        # Save the result
        self.image_preparer.config.write_result = True
        self.image_preparer.config.writing.result_path = "result.fits"

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Loading the images ...")

        # Load image information
        info_path = os.path.join(self.data_path, "info.dat")
        info_table = tables.from_file(info_path)

        # Loop over all files found in the data directory
        for image_path, image_name in filesystem.files_in_path(self.data_path, extension="fits", not_contains="error", names=True):

            # If only a single image must be prepared, check if this image matches the specified image name
            if self.config.single_image is not None and image_name != self.config.single_image: continue

            # Determine the output path for this image
            image_output_path = os.path.join(self.prep_path, image_name)

            # Check whether this image already has a prepared image
            final_path = os.path.join(image_output_path, "result.fits")
            if filesystem.is_file(final_path): continue

            # Get the corresponding index in the information table
            info_index = tables.find_index(info_table, image_name)
            if info_index is None:
                # TEMP: skip image if not defined in table !!
                self.log.warning("No information about " + image_name + ": skipping")
                continue

            # Get image properties such as the unit and the FWHM of the PSF
            unit = u.Unit(info_table["Unit"][info_index]) if not info_table["Unit"].mask[info_index] else None
            fwhm = info_table["FWHM"][info_index] * u.Unit(info_table["FWHM unit"][info_index]) if not info_table["FWHM"].mask[info_index] else None

            # Set the path to the region of bad pixels
            bad_region_path = os.path.join(self.data_path, "bad", image_name + ".reg")
            if not filesystem.is_file(bad_region_path): bad_region_path = None

            # Import the image
            importer = ImageImporter()
            importer.run(image_path, bad_region_path, unit=unit, fwhm=fwhm)

            # Add the image that has to be processed to the list
            self.images.append(importer.image)

            # Clear the image importer
            self.importer.clear()

    # -----------------------------------------------------------------

    def get_attenuations(self):

        """
        This function ...
        :return:
        """

        # Download the extinction table
        table = IrsaDust.get_extinction_table(self.galaxy_name)

        # Loop over all images
        for image in self.images:

            # Get the image name
            image_name = image.name

            # GALEX bands
            if "GALEX" in image.name:

                # Get the A(V) / E(B-V) ratio
                v_band_index = tables.find_index(table, "CTIO V")
                av_ebv_ratio = table["A_over_E_B_V_SandF"][v_band_index]

                # Get the attenuation of the V band A(V)
                attenuation_v = table["A_SandF"][v_band_index]

                # Determine the factor
                if "NUV" in image.name: factor = 8.0
                elif "FUV" in image.name: factor = 7.9
                else: raise ValueError("Unsure which GALEX band this is")

                # Calculate the attenuation
                attenuation = factor * attenuation_v / av_ebv_ratio

                # Set the attenuation
                self.attenuations[image_name] = attenuation

            # Fill in the Ha attenuation manually
            elif "Ha" in image.name: self.attenuations[image_name] = 0.174

            # Other bands for which attenuation is listed by IRSA
            elif image_name in irsa_names:

                irsa_name = irsa_names[image_name]

                # Find the index of the corresponding table entry
                index = tables.find_index(table, irsa_name)

                # Get the attenuation and add it to the dictionary
                attenuation = table["A_SandF"][index]
                self.attenuations[image_name] = attenuation

            # All other bands: set attenuation to zero
            else: self.attenuations[image_name] = 0.0

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Preparing the images ...")

        # Loop over the input images (and masks)
        for image in self.images:

            # Set the attenuation value
            self.image_preparer.config.attenuation = self.attenuations[image.name]

            # Check if there is a file stating the indices of stars that should be treated as exceptions during the extraction procedure
            exceptions_path = os.path.join(self.data_path, "exceptions", image.name + ".txt")
            if filesystem.is_file(exceptions_path): self.image_preparer.config.extraction.exceptions_path = exceptions_path
            else: self.image_preparer.config.extraction.exceptions_path = None

            # If this image is not the reference image, set the appropriate options for rebinning and convolution
            # or this image does not need to be convolved (e.g. SPIRE images)
            if image.name == self.config.reference_image or image.name not in aniano_names:

                self.image_preparer.config.rebin = False
                self.image_preparer.config.convolve = False

            # Images that do need to be convolved
            else:

                # Set the path to the convolution kernel
                this_aniano_name = aniano_names[image.name]
                reference_aniano_name = aniano_names[self.config.reference_image]
                kernel_path = os.path.join(self.kernels_path, "Kernel_HiRes_" + this_aniano_name + "_to_" + reference_aniano_name + ".fits")
                self.image_preparer.config.convolution.kernel_path = kernel_path

                self.image_preparer.config.rebin = True
                self.image_preparer.config.convolve = True

            # Set the path to the noise region
            noise_path = os.path.join(self.data_path, "noise.reg")
            self.image_preparer.config.uncertainties.noise_path = noise_path
            self.image_preparer.config.uncertainties.calibration_error = calibration_errors[image.name]

            # Determine the output path for this image
            image_output_path = os.path.join(self.prep_path, image.name)

            # Create the output directory if it does not exist for this image
            filesystem.create_directory(image_output_path)
            self.image_preparer.config.output_path = image_output_path

            # Run the image preparation
            self.image_preparer.run(image)

            # Clear the image preparer
            self.image_preparer.clear()

# -----------------------------------------------------------------
