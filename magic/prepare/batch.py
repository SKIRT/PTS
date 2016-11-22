#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.batch Contains the BatchImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from multiprocessing import Pool

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.filter import Filter
from ...core.basics.configurable import Configurable
from ..misc.calibration import CalibrationError
from ..misc.extinction import GalacticExtinction
from ..core.frame import Frame
from ..core.dataset import DataSet
from ..region.list import SkyRegionList

# -----------------------------------------------------------------

class BatchImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(BatchImagePreparer, self).__init__(config)

        # The frames
        self.frames = dict()

        # The errormaps
        self.errormaps = dict()

        # The attenuations
        self.attenuations = dict()

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

    # -----------------------------------------------------------------

    def add_frame(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.frames:

            wcs = self.frames[name].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.frames: boxes_region.append(self.frames[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    @property
    def center_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.bounding_box.center

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the galactic attenuation
        self.get_extinction()

        # Getting calibration error
        #CalibrationError.from_filter(image.filter)

        # 3. Prepare the images
        self.prepare()

        # 4. Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchImagePreparer, self).setup()

        # Load the frames
        self.load_frames()

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Create new dataset
        if self.config.dataset.endswith(".fits"):

            # Load the frame
            frame = Frame.from_file(self.config.dataset)

            # Determine the name for this image
            name = str(frame.filter)

            # Add the frame
            self.add_frame(name, frame)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the frames
            self.frames = dataset.get_frames()

        # Invalid value for 'dataset'
        else: raise ValueError("Parameter 'dataset' must be filename of a dataset file (.dat) or a FITS file (.fits)")

    # -----------------------------------------------------------------

    def get_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galactic extinction ...")

        # Create the galactic extinction calculator
        extinction = GalacticExtinction(self.center_coordinate)

        # Loop over the frames
        for label in self.frames:

            # Get the filter
            fltr = self.frames[label].filter

            # Get the extinction
            self.attenuations[label] = extinction.extinction_for_filter(fltr)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the images ...")

        for label in self.frames:

            # Check if the intermediate results have already been produced for this image and saved to the
            # corresponding preparation subdirectory
            extracted_path = fs.join(output_path, "extracted.fits")
            corrected_path = fs.join(output_path, "corrected_for_extinction.fits")
            converted_path = fs.join(output_path, "converted_unit.fits")
            convolved_path = fs.join(output_path, "convolved.fits")
            rebinned_path = fs.join(output_path, "rebinned.fits")
            subtracted_path = fs.join(output_path, "sky_subtracted.fits")

            self.extracted_paths[label] = extracted_path
            self.corrected_paths[label] = corrected_path
            self.converted_paths[label] = converted_path
            self.convolved_paths[label] = convolved_path
            self.rebinned_paths[label] = rebinned_path
            self.subtracted_paths[label] = subtracted_path

        # 1. Extract stars and galaxies from the image
        self.extract_sources()

        # 2. If requested, correct for galactic extinction
        self.correct_for_extinction()

        # 3. If requested, convert the unit
        self.convert_units()

        # 4. If requested, convolve
        self.convolve()

        # 5. Rebin
        self.rebin()

        # 6.
        self.subtract_sky()

        # 7.
        self.calculate_calibration_uncertainties()

        # 8.
        self.calculate_errormaps()

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.frames:

            # Check if the source-extracted image is present
            if fs.is_file(self.extracted_paths[label]): continue

            # Execute
            result = self.pool.apply_async(_extract_sources, args=(,))  # All simple types (strings) ?

            # Add the result
            results.append(result)

            # Get and set the dust masses

        #self.dust_masses = [result.get() for result in results]

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.frames:

            corrected_path = fs.join(output_path, "corrected_for_extinction.fits")

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        for label in self.frames:

            converted_path = fs.join(output_path, "converted_unit.fits")

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        for label in self.frames:

            convolved_path = fs.join(output_path, "convolved.fits")

            if label == self.config.convolution_reference: continue

            if fs.is_file(convolved_path): continue

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # 5. If requested, rebin
        # rebin

        for label in self.frames:

            rebinned_path = fs.join(output_path, "rebinned.fits")

            if label == self.config.rebinning_reference: continue

            if fs.is_file(rebinned_path): continue

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # 6. If requested, subtract the sky
        # subtract_sky

        for label in self.frames:

            subtracted_path = fs.join(output_path, "sky_subtracted.fits")

    # -----------------------------------------------------------------

    def calculate_calibration_uncertainties(self):

        """
        This function ...
        :return:
        """


        for label in self.frames:

            calibration_error = CalibrationError.from_filter(self.frames[label].filter)

        # 7. Calculate the calibration uncertainties
        # calculate_calibration_uncertainties

    # -----------------------------------------------------------------

    def calculate_errormaps(self):

        """
        This function ...
        :return:
        """

        # 8. If requested, set the uncertainties
        # set_uncertainties

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------

def _extract_sources():

    """
    This function ...
    :return:
    """

# -----------------------------------------------------------------
