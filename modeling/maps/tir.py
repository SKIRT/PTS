#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.tir Contains the TIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapMakingComponent
from ...core.basics.log import log
from ...magic.maps.tir.single import SingleBandTIRMapMaker
from ...magic.maps.tir.multi import MultiBandTIRMapMaker
from ...core.filter.filter import parse_filter
from ...magic.core.list import FrameList
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

methods = ["single", "multi"]

# -----------------------------------------------------------------

singleband_filter_names = ["IRAC I4", "MIPS 24mu", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250"]
multiband_filter_names = ["MIPS 24mu", "Pacs blue", "Pacs green", "Pacs red", "SPIRE 250"]

# -----------------------------------------------------------------

class TIRMapMaker(MapMakingComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(TIRMapMaker, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_tir_path

    # -----------------------------------------------------------------

    @property
    def make_single(self):

        """
        This function ...
        :return:
        """

        return "single" in self.config.methods

    # -----------------------------------------------------------------

    @property
    def make_multi(self):

        """
        This function ...
        :return:
        """

        return "multi" in self.config.methods

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make maps based on a single band
        if self.make_single: self.make_maps_single()

        # 3. Make maps based on multiple bands
        if self.make_multi: self.make_maps_multi()

        # 4. Writing
        self.write()

        # 5. Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TIRMapMaker, self).setup(**kwargs)

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

    def load_data_singleband(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        frames = FrameList()
        errormaps = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_singleband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frames.append(frame, fltr)

            # Load the error map
            if not self.has_errormap_for_filter(fltr): continue
            errors = self.get_errormap_for_filter(fltr)
            errormaps.append(errors, fltr)

        # Return the frames and errors maps
        return frames, errormaps

    # -----------------------------------------------------------------

    def load_data_multiband(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the data ...")

        # Frames and error maps
        frames = FrameList()
        errormaps = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_multiband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frames.append(frame, fltr)

            # Load the error map, if present
            if not self.has_errormap_for_filter(fltr): continue
            errors = self.get_errormap_for_filter(fltr)
            errormaps.append(errors, fltr)

        # Return the frames and error maps
        return frames, errormaps

    # -----------------------------------------------------------------

    def make_maps_single(self):

        """
        This function ...
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
        frames, errors = self.load_data_singleband()
        maker.run(frames=frames, errors=errors, maps=current, method_name=method_name, region_of_interest=self.truncation_ellipse)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_maps_multi(self):

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
        frames, errors = self.load_data_multiband()
        maker.run(frames=frames, errors=errors, maps=current, method_name=method_name, region_of_interest=self.truncation_ellipse)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def write(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the origins
        self.write_origins()

        # Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_tir()

# -----------------------------------------------------------------
