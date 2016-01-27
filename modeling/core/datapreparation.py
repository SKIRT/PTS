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

# Import the relevant AstroMagic classes and modules
from ...magic.core import Image

# Import the relevant PTS classes and modules
from .component import ModelingComponent
from .imagepreparation import ImagePreparer
from ...core.tools import filesystem, tables, time

# -----------------------------------------------------------------

class DataPreparer(ModelingComponent):

    """
    This class ...
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
        super(DataPreparer, self).__init__(config)

        # The list of images
        self.images = []

        # Information about the images
        self.attenuations = dict()
        self.aniano_names = dict()

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

        # Set the input and output path
        modeler.config.input_path = arguments.input_path
        modeler.config.output_path = arguments.output_path

        # A single image can be specified so the preparation is only run with that image
        modeler.config.single_image = arguments.image

        # Set logging path
        if arguments.report: modeler.config.logging.path = os.path.join(modeler.config.output_path, time.unique_name("log") + ".txt")

        # Return the new instance
        return modeler

    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function runs the data preparation ...
        :return:
        """

        # 1. Call the setup function
        self.setup(path)

        # 2. Load the input data
        self.load_data()

        # 2. Prepare the images
        self.prepare()

    # -----------------------------------------------------------------

    def setup(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # -- Children --

        # Create the preparation object
        self.add_child("preparer", ImagePreparer)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(DataPreparer, self).setup(path)

        # -- Attributes --

        # Create the prep path if it does not exist yet
        filesystem.create_directory(self.prep_path)

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
        for file_path in filesystem.files_in_path(self.data_path, extension="fits", not_contains="error"):

            # Get the file name
            file_name = os.path.basename(file_path)
            image_name = os.path.splitext(file_name)[0]

            # If only a single image must be prepared, check if this image matches the specified image name
            if self.config.single_image is not None and image_name != self.config.single_image: continue

            # Determine the output path for this image
            image_output_path = os.path.join(self.prep_path, image_name)

            # Get the corresponding index in the information table
            info_index = tables.find_index(info_table, image_name)
            if info_index is None: continue  # TEMP: skip image if not defined in table !!

            # Check whether this image already has a prepared image
            final_path = os.path.join(image_output_path, "result.fits")
            if filesystem.is_file(final_path): continue

            # Inform the user
            self.log.debug("Importing the " + image_name + " image")

            # Open the image
            image = Image(file_path)
            self.load_error_frame_and_select(image, image_name)

            # Set image properties such as the unit and the FWHM of the PSF
            unit = u.Unit(info_table["Unit"][info_index])
            fwhm = info_table["FWHM"][info_index] * u.Unit(info_table["Unit of FWHM"][info_index]) if not info_table["FWHM"].mask[info_index] else None
            image.set_unit(unit)
            image.set_fwhm(fwhm)

            # Add the image that has to be processed to the list
            self.images.append(image)

            # Set the attenuation
            self.attenuations[image.name] = info_table["Attenuation"][info_index]

            # Set the name for the Aniano kernel
            self.aniano_names[image.name] = info_table["Aniano name"][info_index]

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Preparing the images ...")

        # Loop over the input images
        for image in self.images:

            # Set the attenuation value
            self.preparer.config.attenuation = self.attenuations[image.name]

            # Set extra mask
            extra_path = os.path.join(self.data_path, "extra", image.name + ".reg")
            # Check if a region is present
            # Check whether the extra region file is present
            if filesystem.is_file(extra_path): self.image_preparer.config.extra_region_path = extra_path
            else: self.image_preparer.config.extra_region_path = None

            # If this image is not the reference image, set the appropriate options for rebinning and convolution
            if image.name != self.config.reference:

                # Set the path to the convolution kernel
                this_aniano_name = self.aniano_names[image.name]
                reference_aniano_name = self.aniano_names[self.config.reference]
                kernel_path = os.path.join(self.kernels_path, "Kernel_HiRes_" + this_aniano_name + "_to_" + reference_aniano_name + ".fits")
                self.image_preparer.config.convolution.kernel_path = kernel_path

                self.image_preparer.config.rebin = True
                self.image_preparer.config.convolve = True

            else:

                self.image_preparer.config.rebin = True
                self.image_preparer.config.convolve = True

            # Set the path to the noise region
            noise_path = os.path.join(self.data_path, "noise.reg")
            self.image_preparer.config.uncertainties.noise_path = noise_path

            # Determine the output path for this image
            image_output_path = os.path.join(self.prep_path, image.name)

            # Create the output directory if it does not exist for this image
            filesystem.create_directory(image_output_path)
            self.image_preparer.config.output_path = image_output_path

            # Run the image preparation
            self.preparer.run(image)

            # Clear the image preparer
            self.preparer.clear()

# -----------------------------------------------------------------
