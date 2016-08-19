#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.analyser Contains the MosaicAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from .component import DataComponent
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class MosaicAnalyser(DataComponent):
    
    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(MosaicAnalyser, self).__init__(config)

        # The task which has created the mosaics
        self.task = None

        # The results
        self.mosaics = dict()
        self.poisson_frames = dict()
        self.relative_poisson_frames = dict()

        # The references
        self.references = dict()

        # Origin of the mosaics
        self.origin = None

        # ...
        self.check_paths = dict()

    # -----------------------------------------------------------------

    @classmethod
    def for_task(cls, task):

        """
        This function ...
        :param task:
        :return:
        """

        # Create the instance
        analyser = cls()

        # Set the modeling path as the working path for this class
        analyser.config.path = task.analysis_info["modeling_path"]

        # Set the task
        analyser.task = task

        # Return the instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the results from the mosaicing
        self.load_results()

        # 3. Load the image as obtained from other source
        self.load_references()

        # Rebin if necessary
        self.rebin()

        # 4. Make directory for each band that has been processed
        self.make_directories()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def load_results(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the result of the mosaicing procedure ...")

        # Loop over the FITS files found in the output directory
        for path, name in fs.files_in_path(self.task.local_output_path, extension="fits", contains=self.ngc_id_nospaces, returns=["path", "name"]):

            # Split
            splitted = name.split("_")

            assert splitted[0] == self.ngc_id_nospaces

            if self.origin is None: self.origin = splitted[1]
            else: assert self.origin == splitted[1]
            band = splitted[2]

            # Determine band id
            band_id = self.origin + "_" + band

            # Check type
            if splitted[-1] == "errors": self.poisson_frames[band_id] = Frame.from_file(path)
            elif splitted[-1] == "relerrors": self.relative_poisson_frames[band_id] = Frame.from_file(path)
            elif splitted[-1] == "swarp": continue
            else: self.mosaics[band_id] = Frame.from_file(path) # the mosaic

    # -----------------------------------------------------------------

    def load_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the reference mosaics ...")

        # Loop over the bands
        for band_id in self.mosaics:

            # Determine the path to the image downloaded from the DustPedia archive
            path = fs.join(self.data_images_paths[self.origin], self.ngc_id_nospaces + "_" + band_id + ".fits")

            # Load the reference image
            self.references[band_id] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning to the reference WCS if necessary ...")

        # Loop over the bands
        for band_id in self.mosaics:

            reference_wcs = self.references[band_id].wcs
            mosaic_wcs = self.mosaics[band_id].wcs

            if reference_wcs == mosaic_wcs: continue

            # Rebin mosaic
            self.mosaics[band_id] /= self.mosaics[band_id].pixelarea.to("sr").value
            self.mosaics[band_id].rebin(reference_wcs)
            self.mosaics[band_id] *= self.mosaics[band_id].pixelarea.to("sr").value

            # Rebin poisson frame
            self.poisson_frames[band_id] /= self.poisson_frames[band_id].pixelarea.to("sr").value
            self.poisson_frames[band_id].rebin(reference_wcs)
            self.poisson_frames[band_id] *= self.poisson_frames[band_id].pixelarea.to("sr").value

            # Rebin relative poisson frame
            self.relative_poisson_frames[band_id].rebin(reference_wcs)

    # -----------------------------------------------------------------

    def make_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making directories for checking the mosaicing ...")

        # Loop over the bands
        for band_id in self.mosaics:

            # Set the path
            check_band_path = fs.join(self.data_images_paths[self.origin], "check_" + band_id)
            self.check_paths[band_id] = check_band_path

            # Create the directory
            fs.create_directory(check_band_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write poisson error maps
        self.write_poisson()

        # For checking :

        # Write mosaics
        self.write_mosaics()

        # Write relative poisson maps
        self.write_relative_poisson()

        # Write differences with reference
        self.write_difference()

        # Write relative differences with reference
        self.write_relative_difference()

    # -----------------------------------------------------------------

    def write_mosaics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mosaics created by PTS ...")

        # Loop over the mosaic frames
        for band_id in self.mosaics:

            # Determine path
            path = fs.join(self.check_paths[band_id], self.ngc_id_nospaces + "_" + band_id + ".fits")

            # Save the mosaic frame
            self.mosaics[band_id].save(path)

    # -----------------------------------------------------------------

    def write_poisson(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the poisson error maps ...")

        # Loop over the frames
        for band_id in self.poisson_frames:

            # Determine path
            path = fs.join(self.data_images_paths[self.origin], self.ngc_id_nospaces + "_" + band_id + "_poisson.fits")

            # For GALEX
            if "GALEX" in band_id:

                # Construct the poisson frame in Jy/pix from the relative poisson frame
                poisson_frame_jy = self.relative_poisson_frames[band_id] * self.references[band_id]
                poisson_frame_jy.unit = "Jy/pix"

                # Save the poisson error frame
                poisson_frame_jy.save(path)

            # For SDSS
            else:

                # Save the poisson error frame
                self.poisson_frames[band_id].save(path)

    # -----------------------------------------------------------------

    def write_relative_poisson(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the relative poisson error maps ...")

        # Loop over the frames
        for band_id in self.relative_poisson_frames:

            # Determine path
            path = fs.join(self.check_paths[band_id], self.ngc_id_nospaces + "_" + band_id + "_relpoisson.fits")

            # Save the relative poisson error frame
            self.relative_poisson_frames[band_id].save(path)

    # -----------------------------------------------------------------

    def write_difference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the difference with the reference (DustPedia) mosaics ...")

        # Loop over the bands
        for band_id in self.references:

            # Determine path
            path = fs.join(self.check_paths[band_id], self.ngc_id_nospaces + "_" + band_id + "_difference.fits")

            # Calculate the difference
            difference = self.references[band_id] - self.mosaics[band_id]

            # Write
            difference.save(path)

    # -----------------------------------------------------------------

    def write_relative_difference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the relative difference with the reference (DustPedia) mosaics ...")

        # Loop over the bands
        for band_id in self.references:

            # Determine path
            path = fs.join(self.check_paths[band_id], self.ngc_id_nospaces + "_" + band_id + "_reldifference.fits")

            # Calculate the difference
            reldifference = self.references[band_id] / self.mosaics[band_id]

            # Write
            reldifference.save(path)

# -----------------------------------------------------------------
