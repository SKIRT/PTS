#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.single Contains the SingleImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class SingleImagePreparer(Configurable):

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
        super(SingleImagePreparer, self).__init__(config)

        # -- Attributes --

        self.attenuation = None
        self.calibration = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Get the galactic extinction
        self.get_extinction()

        # 3. Get calibration error
        self.get_calibration()

        # 4. Get convolution kernel
        self.get_kernel()

        # 5. Get rebinning WCS
        self.get_wcs()

        # 6. Prepare the image
        self.prepare()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SingleImagePreparer, self).setup()

    # -----------------------------------------------------------------

    def get_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # Get the galactic extinction for this image
        self.attenuation = GalacticExtinction(center_coordinate).extinction_for_filter(image.filter)

    # -----------------------------------------------------------------

    def get_calibration(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # Get the calibration error
        self.calibration = CalibrationError.from_filter(image.filter)

    # -----------------------------------------------------------------

    def get_kernel(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking up the necessary kernel file ...")

        # Get the filter to which to convolve to
        convolve_to_filter = Filter.from_string(arguments.convolve_to)

        # Create an AnianoKernels instance
        kernels = AnianoKernels()

        # Get the path to the appropriate convolution kernel
        kernel_path = kernels.get_kernel_path(image.filter, convolve_to_filter, fwhm=fwhm)

        # Set the kernel path
        #arguments.kernel = kernel_path
        self.kernel = kernel_path

    # -----------------------------------------------------------------

    def get_wcs(self):

        """
        This function ...
        :return:
        """

        # Determine the absolute path to the reference image
        self.rebin_to = fs.absolute(arguments.rebin_to)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Create an ImagePreparer instance
        preparer = ImagePreparer.from_arguments(arguments)

        # Run the image preparation
        preparer.run(image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments, visualisation_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
