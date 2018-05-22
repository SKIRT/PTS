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
from ...core.basics.log import log
from .component import DataComponent
from ...magic.core.frame import Frame, get_filter
from ...magic.core.image import Image
from ...core.launch.pts import load_task
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

class MosaicAnalyser(DataComponent):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MosaicAnalyser, self).__init__(*args, **kwargs)

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the results from the mosaicing
        if not self.has_results: self.load_results()

        # Check results
        self.check_results()

        # 3. Load the image as obtained from other source
        self.load_references()

        # 4. Rebin if necessary
        self.rebin()

        # 5. Make directory for each band that has been processed
        self.make_directories()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(MosaicAnalyser, self).setup(**kwargs)

        # Check for image that is specified
        if "image" in kwargs: self.setup_from_image(kwargs.pop("image"))

        # Check for image path that is specified
        elif self.config.image_path is not None:

            # Load the image
            image = Image.from_file(self.config.image_path)
            self.setup_from_image(image, self.config.band_id)

        # From directory with multiple images
        if self.config.images_path is not None: self.setup_from_images_path()

        # From seperate mosaic, errors and relerrors files
        if "mosaic" in kwargs and "errors" in kwargs and "relerrors" in kwargs: self.setup_from_mosaics(kwargs.pop("mosaic"), kwargs.pop("errors"), kwargs.pop("relerrors"))
        elif self.config.out_path is not None: self.setup_from_out(self.config.out_path)

        # Check whether task is specified
        if "task" in kwargs: task = kwargs.pop("task")
        elif self.config.host_id is not None:

            if self.config.task_id is None: raise ValueError("Task ID is not specified")

            # Load the task
            task = load_task(self.config.host_id, self.config.task_id)

        elif self.config.task_id is not None: raise ValueError("Task ID is specified but host ID is not")
        else: task = None

        # If the task is not None
        if task is not None: self.task = task

    # -----------------------------------------------------------------

    def setup_from_image(self, image, band_id=None):

        """
        This function ...
        :param image:
        :param band_id: 
        :return: 
        """

        # Get band id
        if image.filter is not None: band_id = image.filter_name.replace(" ", "_")
        if band_id is None: raise ValueError("Band ID must be specified")

        # Get frames
        mosaic_frame = image.frames["primary"]
        mosaic_errors = image.frames["errors"]

        self.mosaics[band_id] = mosaic_frame
        self.poisson_frames[band_id] = mosaic_errors

        self.origin = band_id.split("_")[0]

    # -----------------------------------------------------------------

    def setup_from_images_path(self):

        """
        This function ...
        :return: 
        """

        # Find files
        paths = fs.files_in_path(self.config.images_path, extension="fits")

        origin = None

        # Loop over the paths
        for path in paths:

            name = fs.strip_extension(fs.name(path))

            fltr = get_filter(path)
            if fltr is None:
                if "mosaic_jansky" in name: fltr = parse_filter(name.split("mosaic_jansky_")[1])
                else: raise NotImplementedError("Don't know how to proceed")

            band_id = str(fltr).replace(" ", "_")
            origin_i = band_id.split("_")[0]

            if origin is None: origin = origin_i
            elif origin_i != origin: raise IOError("Found files of different origin (instrument)")

            # Load image
            image = Image.from_file(path)

            # Get frames
            mosaic_frame = image.frames["primary"]
            mosaic_errors = image.frames["errors"]

            # Set maps
            self.mosaics[band_id] = mosaic_frame
            self.poisson_frames[band_id] = mosaic_errors
            if "relerrors" in image.frames: self.relative_poisson_frames[band_id] = image.frames["relerrors"]

        # Set origin
        self.origin = origin

    # -----------------------------------------------------------------

    def setup_from_mosaics(self, mosaic, errors, relerrors=None):

        """
        This funciton ...
        :param mosaic: 
        :param errors: 
        :param relerrors: 
        :return: 
        """

        # Check filter
        if mosaic.filter is None: raise ValueError("Mosaic filter is not defined")
        if errors.filter is not None and mosaic.filter != errors.filter: raise ValueError("Filter of mosaic is not the same as filter of errors")
        if relerrors is not None and relerrors.filter is not None and mosaic.filter != relerrors.filter: raise ValueError("Filter of mosaic is not the same as filter of relerrors")

        # Determine the band ID
        band_id = mosaic.filter_name.replace(" ", "_")

        # Set
        self.mosaics[band_id] = mosaic
        self.poisson_frames[band_id] = errors
        if relerrors is not None: self.relative_poisson_frames[band_id] = relerrors

        # Set origin
        self.origin = band_id.split("_")[0]

    # -----------------------------------------------------------------

    def setup_from_out(self, out_path):

        """
        This function ...
        :param out_path: 
        :return: 
        """

        # Find files
        paths = fs.files_in_path(out_path, extension="fits")

        # Determine origin
        origin = None
        for path in paths:
            name = fs.strip_extension(fs.name(path))
            origin_i = name.split("_")[1]
            if origin is None: origin = origin_i
            elif origin != origin_i: raise IOError("Found files of different origin (instrument)")

        # Set the origin
        self.origin = origin

        # Loop over the files
        for path in paths:

            name = fs.strip_extension(fs.name(path))
            band_id = origin + "_" + name.split("_")[2]

            maptype = name.split("_")[-1]

            # Check type
            if maptype == "errors": self.poisson_frames[band_id] = Frame.from_file(path)
            elif maptype == "relerrors": self.relative_poisson_frames[band_id] = Frame.from_file(path)
            elif maptype == "swarp": continue
            else: self.mosaics[band_id] = Frame.from_file(path)  # the mosaic

    # -----------------------------------------------------------------

    @property
    def has_results(self):

        """
        This function ...
        :return: 
        """

        return len(self.mosaics) > 0

    # -----------------------------------------------------------------

    def load_results(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the result of the mosaicing procedure ...")

        # Loop over the FITS files found in the output directory
        for path, name in fs.files_in_path(self.task.local_output_path, extension="fits", contains=self.ngc_name_nospaces, returns=["path", "name"]):

            # Split
            splitted = name.split("_")

            assert splitted[0] == self.ngc_name_nospaces

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

    def check_results(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the results ...")

        # Loop over the bands
        for band in self.mosaics:

            # Check whether errors is present
            if band not in self.poisson_frames: raise ValueError("Poisson frame for " + band + " is not found")

            # Check whether relative errors is present
            if band not in self.relative_poisson_frames:

                mosaic_errors = self.poisson_frames[band]
                mosaic_frame = self.mosaics[band]

                # Create the relative poisson frame
                # Calculate the relative error map
                relerrors = mosaic_errors / mosaic_frame
                relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
                relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

                # Add the relative errors frame
                self.relative_poisson_frames[band] = relerrors

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
            path = fs.join(self.data_images_paths[self.origin], self.ngc_name_nospaces + "_" + band_id + ".fits")

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
            path = fs.join(self.check_paths[band_id], self.ngc_name_nospaces + "_" + band_id + ".fits")

            # Save the mosaic frame
            self.mosaics[band_id].saveto(path)

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
            path = fs.join(self.data_images_paths[self.origin], self.ngc_name_nospaces + "_" + band_id + "_poisson.fits")

            # For GALEX
            if "GALEX" in band_id:

                # Construct the poisson frame in Jy/pix from the relative poisson frame
                poisson_frame_jy = self.relative_poisson_frames[band_id] * self.references[band_id]
                poisson_frame_jy.unit = "Jy/pix"

                # Save the poisson error frame
                poisson_frame_jy.saveto(path)

            # For SDSS
            else:

                # Save the poisson error frame
                self.poisson_frames[band_id].saveto(path)

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
            path = fs.join(self.check_paths[band_id], self.ngc_name_nospaces + "_" + band_id + "_relpoisson.fits")

            # Save the relative poisson error frame
            self.relative_poisson_frames[band_id].saveto(path)

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
            path = fs.join(self.check_paths[band_id], self.ngc_name_nospaces + "_" + band_id + "_difference.fits")

            # Calculate the difference
            difference = self.references[band_id] - self.mosaics[band_id]

            # Write
            difference.saveto(path)

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
            path = fs.join(self.check_paths[band_id], self.ngc_name_nospaces + "_" + band_id + "_reldifference.fits")

            # Calculate the difference
            reldifference = self.references[band_id] / self.mosaics[band_id]

            # Write
            reldifference.saveto(path)

# -----------------------------------------------------------------
