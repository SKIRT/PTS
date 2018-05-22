#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.core.test.implementation import TestImplementation
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.modeling.tests.base import m81_data_path
from pts.core.filter.filter import parse_filter
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.core.frame import Frame
from pts.magic.core.list import FrameList
from pts.core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

description = "Test homogenization of images with different PSF and coordinate system"

# -----------------------------------------------------------------

filter_names = ["GALEX FUV", "IRAC I3", "MIPS 160", "Planck_550"]

# -----------------------------------------------------------------

# Determine the path to the headers directory
headers_path = fs.join(m81_data_path, "headers")

# -----------------------------------------------------------------

def get_coordinate_system(fltr):

    """
    This function ...
    :param fltr: 
    :return: 
    """

    for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):
        header_fltr = parse_filter(name)
        if header_fltr == fltr: return CoordinateSystem.from_header_file(path)
    return None

# -----------------------------------------------------------------

class HomogenizeTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HomogenizeTest, self).__init__(*args, **kwargs)

        # The galaxy name
        #self.galaxy_name = None

        # Frames
        self.frames = None

        # DustPedia database
        self.database = None

        # Combinations
        self.first = None
        self.second = None

        # The distance
        self.distance = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the frames
        self.get_frames()

        # 3. Combine
        self.combine()

        # 4. Test the result
        self.test()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(HomogenizeTest, self).setup(**kwargs)

        # Set galaxy name
        #sample = DustPediaSample()
        #self.galaxy_name = sample.get_name(self.config.galaxy)

        # Connect to the database
        #self.database = DustPediaDatabase()
        #username, password = get_account()
        #self.database.login(username, password)

        # Set distance
        self.distance = 10. * u("Mpc")

        #print(self.database.get_image_urls(self.galaxy_name))

    # -----------------------------------------------------------------

    def get_frames(self):

        """
        This function ...
        :return: 
        """

        # Inform theb user
        log.info("Getting the frames ...")

        # Get the frames
        #self.frames = self.database.get_framelist_for_filters(self.galaxy_name, filter_names)

        self.frames = FrameList()

        # Loop over the filter names
        for filter_name in filter_names:

            # Parse filter
            fltr = parse_filter(filter_name)

            # Get wcs
            wcs = get_coordinate_system(fltr)

            # Generate frame
            #frame = Frame.random(wcs.shape, wcs=wcs, filter=fltr)
            frame = Frame.ones(wcs.shape, wcs=wcs, filter=fltr)
            frame.unit = "Jy"

            # Add the frame
            self.frames.append(frame)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    @property
    def highest_resolution_filters(self):

        """
        This function ...
        :return: 
        """

        return sorted(self.filters, key=lambda fltr: self.frames[fltr].pixelscale.average)[0:2]

    # -----------------------------------------------------------------

    @property
    def lowest_resolution_filters(self):

        """
        This function ...
        :return: 
        """

        return sorted(self.filters, key=lambda fltr: self.frames[fltr].pixelscale.average)[2:]

    # -----------------------------------------------------------------

    def combine(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Combining ...")

        # Combine the first pair
        self.combine_first()

        # Combine the second pair
        self.combine_second()

        # Combine last
        self.combine_last()

    # -----------------------------------------------------------------

    def combine_first(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Doing first combination ...")

        fltr_a, fltr_b = self.highest_resolution_filters
        frames = self.frames[fltr_a, fltr_b]

        frames.convolve_to_highest_fwhm()
        frames.rebin_to_highest_pixelscale()
        frames.convert_to_same_unit(unit="W/micron/m2") # here conversion related to frequency/wavelength is necessary
        frames.convert_to_same_unit(unit="W/m2", density=True) # here conversion from spectral flux density to flux is necessary

        # Set the frames
        #self.first = frames

        # Now do something, like adding them
        self.first = np.random.random() * frames[fltr_a].get_log10() + np.random.random() * frames[fltr_b].get_log10()

        # Set the unit
        self.first.unit = "W/m2"

    # -----------------------------------------------------------------

    def combine_second(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Doing second combination ...")

        fltr_a, fltr_b = self.lowest_resolution_filters
        frames = self.frames[fltr_a, fltr_b]

        frames.convolve_to_highest_fwhm()
        frames.rebin_to_highest_pixelscale()
        frames.convert_to_same_unit(unit="MJy/sr") # here conversion to surface brightness, related to the pixelscale is necessary

        # Set the frames
        #self.second = frames

        frames.convert_to_same_unit(unit="Lsun", density=True, distance=self.distance) # convert to solar luminosity

        # Now do something
        self.second = np.random.random() * frames[fltr_a] + np.random.random() * frames[fltr_b]

        # Set unit
        self.second.unit = "Lsun"

    # -----------------------------------------------------------------

    def combine_last(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Doing last combination")

    # -----------------------------------------------------------------

    def test(self):

        """
        This fucntion ...
        :return: 
        """

        # Inform the user
        log.info("Testing ...")

# -----------------------------------------------------------------