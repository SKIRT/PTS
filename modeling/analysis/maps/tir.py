#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.tir Contains the TIRMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....magic.maps.tir.single import SingleBandTIRMapMaker
from ....magic.maps.tir.multi import MultiBandTIRMapMaker
from ....magic.core.list import FrameList
from ....core.tools.utils import lazyproperty
from ...maps.tir import singleband_filter_names, multiband_filter_names
from ....core.filter.filter import parse_filter
from .component import MapsAnalysisComponent

# -----------------------------------------------------------------

class TIRMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TIRMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the maps
        self.make_maps()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TIRMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def available_filters_singleband(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the colours
        for filter_name in singleband_filter_names:

            # Parse fltr
            fltr = parse_filter(filter_name)

            # If no image is avilalbe for this filters, skip
            if not self.has_frame_for_filter(fltr): continue

            # otherwise, add to the list of filters
            filters.append(fltr)

        # Return the available filters
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def available_filters_multiband(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the colours
        for filter_name in multiband_filter_names:

            # Parse filter
            fltr = parse_filter(filter_name)

            # If no image is avilalbe for this filters, skip
            if not self.has_frame_for_filter(fltr): continue

            # otherwise, add to the list of filters
            filters.append(fltr)

        # Return the available filters
        return filters

    # -----------------------------------------------------------------

    def load_frames_singleband(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        frames = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_singleband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frame.distance = self.galaxy_distance
            frames.append(frame, fltr)

        # Return the frames
        return frames

    # -----------------------------------------------------------------

    def load_frames_multiband(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        # Frames and error maps
        frames = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_multiband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frame.distance = self.galaxy_distance
            frames.append(frame, fltr)

        # Return the frames
        return frames

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR map ...")

        # 1. Make maps based on a single band
        self.make_tir_maps_single()

        # 2. Make maps based on multiple bands
        self.make_tir_maps_multi()

    # -----------------------------------------------------------------

    def make_tir_maps_single(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on a single band ...")

        # Set the method name
        method_name = "single"

        # Create
        maker = SingleBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames = self.load_frames_singleband()
        maker.run(frames=frames, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_tir_maps_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on multiple bands ...")

        # Set method name
        method_name = "multi"

        # Create
        maker = MultiBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames = self.load_frames_multiband()
        maker.run(frames=frames, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.tir_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

# -----------------------------------------------------------------
