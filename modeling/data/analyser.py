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

        # The references
        self.references = dict()

        # Origin of the mosaics
        self.origin = None

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

        # 4. Verify the mosaicing result with the reference
        self.verify()

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
            elif splitted[-1] == "relerrors": continue
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

        for band_id in self.mosaics:

            # Determine the path to the image downloaded from the DustPedia archive
            path = fs.join(self.data_images_paths[self.origin], self.ngc_id_nospaces + "_" + band_id + ".fits")

            # Load the reference image
            self.references[band_id] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def verify(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Verifying the mosaiced images and poisson frame ...")

        # Convert relative poisson frames into absolute poisson errors in Jansky
        for band_id in self.poisson_frames:

            # Calculate absolute poisson frame
            self.poisson_frames[band_id] *= self.references[band_id]

            # Set the unit of the poisson frame
            self.poisson_frames[band_id].unit = self.references[band_id].unit

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write mosaics
        self.write_mosaics()

        # Write poisson error maps
        self.write_poisson()

    # -----------------------------------------------------------------

    def write_mosaics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mosaics created by PTS ...")

        mosaic_pts_path = fs.join(self.data_images_paths[self.origin], "pts")
        fs.create_directory(mosaic_pts_path)

        # Loop over the mosaic frames
        for band_id in self.mosaics:

            # Determine path
            path = fs.join(mosaic_pts_path, self.ngc_id_nospaces + "_" + band_id + ".fits")

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

            # Save the poisson error frame
            self.poisson_frames[band_id].save(path)

# -----------------------------------------------------------------
