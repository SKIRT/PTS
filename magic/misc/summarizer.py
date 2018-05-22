#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.summerizer Contains the ImageSummarizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.basics.configurable import Configurable
from ...core.tools import formatting as fmt
from ..core.image import Image

# -----------------------------------------------------------------

class ImageSummarizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImageSummarizer, self).__init__(*args, **kwargs)

        self.image_paths = []

        self.pixelscales = dict()
        self.fwhms = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Find the simulations
        self.find()

        # List the simulations
        if self.config.list: self.list()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding images ...")

        search_path = self.config.path

        # Search for ski files
        for path, name in fs.files_in_path(search_path, extension="fits", recursive=self.config.recursive, returns=["path", "name"], exact_name="initialized"):

            # Determine the prefix
            #prefix = name

            #dirpath = fs.directory_of(path)

            # Add the path
            #ski_files[prefix].append(dirpath)

            self.image_paths.append(path)

            image = Image.from_file(path)

            pixelscale = image.average_pixelscale
            fwhm = image.fwhm

            self.pixelscales[path] = pixelscale
            self.fwhms[path] = fwhm

    # -----------------------------------------------------------------

    @property
    def image_paths_ordered_by_pixelscale(self):

        """
        This function ...
        :return:
        """

        return sorted(self.image_paths, key=lambda path: self.pixelscales[path])

    # -----------------------------------------------------------------

    @property
    def image_paths_ordered_by_fwhm(self):

        """
        This function ...
        :return:
        """

        return sorted(self.image_paths, key=lambda path: self.fwhms[path])

    # -----------------------------------------------------------------

    def list(self):

        """
        This function ...
        :return:
        """

        print(fmt.green + " images" + fmt.reset + ":")

        print("")

        counter = 1

        if self.config.pixelscale: paths = self.image_paths_ordered_by_pixelscale
        elif self.config.fwhm: paths = self.image_paths_ordered_by_fwhm
        else: paths = self.image_paths

        for path in paths:

            relpath = path.split(self.config.path)[1][1:]

            print(fmt.underlined + relpath + fmt.reset)
            print("")

            image = Image.from_file(path)

            nframes = image.nframes
            nregions = image.nregions
            nmasks = image.nmasks
            nsegments = image.nsegments

            nx = image.xsize
            ny = image.ysize

            pixelscale = image.average_pixelscale
            fwhm = image.fwhm

            fltr = image.filter

            print(" * " + fmt.bold + "Number of frames: " + fmt.reset + str(nframes))
            print(" * " + fmt.bold + "Number of regions: " + fmt.reset + str(nregions))
            print(" * " + fmt.bold + "Number of masks: " + fmt.reset + str(nmasks))
            print(" * " + fmt.bold + "Number of segmentation maps: " + fmt.reset + str(nsegments))
            print(" * " + fmt.bold + "Number of x pixels: " + fmt.reset + str(nx))
            print(" * " + fmt.bold + "Number of y pixels: " + fmt.reset + str(ny))
            print(" * " + fmt.bold + "Pixelscale: " + fmt.reset + str(pixelscale))
            print(" * " + fmt.bold + "FWHM: " + fmt.reset + str(fwhm))
            print(" * " + fmt.bold + "Filter: " + fmt.reset + str(fltr))
            print(" * " + fmt.bold + "Sky subtracted: " + fmt.reset + str(image.sky_subtracted))
            print(" * " + fmt.bold + "Sources extracted: " + fmt.reset + str(image.source_extracted))
            print(" * " + fmt.bold + "Extinction corrected: " + fmt.reset + str(image.extinction_corrected))
            print("")

            counter += 1

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
