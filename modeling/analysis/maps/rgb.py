#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.rgb Contains the RGBMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.visualization import make_lupton_rgb

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from .component import MapsAnalysisComponent
from ....core.filter.filter import parse_filter
from ....core.tools.strings import find_delimiter
from ....core.tools.utils import lazyproperty
from ....magic.core.image import Image

# -----------------------------------------------------------------

uv_rgb_strings = ["FUV-NUV-u"]
optical_rgb_strings = ["u-g-i", "u-g-r", "u-g-z"]
fir_rgb_strings = ["Mips24-Pacs70-Pacs100", "Mips24-Pacs70-Pacs160", "Pacs70-Pacs100-Pacs160", "SPIRE250-SPIRE350-SPIRE500"]

# -----------------------------------------------------------------

class RGBMapsAnalyser(MapsAnalysisComponent):

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
        super(RGBMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # Make the maps
        self.make_maps()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RGBMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB map ...")

        # UV
        self.make_uv_maps()

        # OPTICAL
        self.make_optical_maps()

        # FIR
        self.make_fir_maps()

    # -----------------------------------------------------------------

    @lazyproperty
    def available_uv_rgb_names(self):

        """
        This function ...
        :return:
        """

        names = []

        # Loop over the strings
        for string in uv_rgb_strings:

            # Get the filters
            for fltr in get_filters(string, "-"):

                # If either one of the images is not available, we can not calculate the RGB
                if not self.simulated_dataset.has_frame_for_filter(fltr): break

            # Break not encountered
            else: names.append(string)

        # Return names
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def available_optical_rgb_names(self):

        """
        This function ...
        :return:
        """

        names = []

        # Loop over the strings
        for string in optical_rgb_strings:

            # Get the filters
            for fltr in get_filters(string, "-"):

                # If either one of the images is not available, we cannot calculate the RGB
                if not self.simulated_dataset.has_frame_for_filter(fltr): break

            # Break not encountered
            else: names.append(string)

        # Return names
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def available_fir_rgb_names(self):

        """
        This function ...
        :return:
        """

        names = []

        # Loop over the strings
        for string in fir_rgb_strings:

            # Get the filters
            for fltr in get_filters(string, "-"):

                # If either one of the images is not available, we cannot calculate the RGB
                if not self.simulated_dataset.has_frame_for_filter(fltr): break

            # Break not encountered
            else: names.append(string)

        # Return names
        return names

    # -----------------------------------------------------------------

    def make_uv_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB maps based on UV images ...")

        # Set method name
        method_name = "UV"

        # The RGB maps
        maps = dict()
        origins = dict()

        # Loop over the names
        for name in self.available_uv_rgb_names:

            # Get the filters
            fltr_a, fltr_b, fltr_c = get_filters(name, "-")

            # Load and uniformize to same resolution and unit of Jansky
            a, b, c = self.get_frames_for_filters(fltr_a, fltr_b, fltr_c, uniformize=True, unit="Jy")

            # Make the map
            rgb = make_lupton_rgb(a, b, c, Q=self.config.uv_softening, stretch=self.config.uv_stretch)
            rgb_image = Image.from_3d_array(rgb, name=name, wcs=a.wcs, names=["r", "g", "b"])

            # Add the map
            maps[name] = rgb_image
            origins[name] = [fltr_a, fltr_b, fltr_c]

        # Set the maps
        self.maps[method_name] = maps

        # Set the origins
        self.origins[method_name] = origins

        # Set the methods
        #self.methods[method_name] = methods

    # -----------------------------------------------------------------

    def make_optical_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB maps based on optical images ...")

        # Set method name
        method_name = "Optical"

        # The RGB maps
        maps = dict()
        origins = dict()

        # Loop over the names
        for name in self.available_optical_rgb_names:

            # Get the filters
            fltr_a, fltr_b, fltr_c = get_filters(name, "-")

            # Load and uniformize to same resolution and unit of Jansky
            a, b, c = self.get_frames_for_filters(fltr_a, fltr_b, fltr_c, uniformize=True, unit="Jy")

            # Make the map
            rgb = make_lupton_rgb(a, b, c, Q=self.config.optical_softening, stretch=self.config.optical_stretch)
            rgb_image = Image.from_3d_array(rgb, name=name, wcs=a.wcs, names=["r", "g", "b"])

            # Add the map
            maps[name] = rgb_image
            origins[name] = [fltr_a, fltr_b, fltr_c]

        # Set the maps
        self.maps[method_name] = maps

        # Set the origins
        self.origins[method_name] = origins

        # Set the methods
        # self.methods[method_name] = methods

    # -----------------------------------------------------------------

    def make_fir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB maps based on FIR images ...")

        # Set method name
        method_name = "TIR"

        # The RGB maps
        maps = dict()
        origins = dict()

        # Loop over the names
        for name in self.available_fir_rgb_names:

            # Get the filters
            fltr_a, fltr_b, fltr_c = get_filters(name, "-")

            # Load and uniformize to same resolution and unit of Jansky
            a, b, c = self.get_frames_for_filters(fltr_a, fltr_b, fltr_c, uniformize=True, unit="Jy")

            # Make the map
            rgb = make_lupton_rgb(a, b, c, Q=self.config.fir_softening, stretch=self.config.optical_stretch)
            rgb_image = Image.from_3d_array(rgb, name=name, wcs=a.wcs, names=["r", "g", "b"])

            # Add the map
            maps[name] = rgb_image
            origins[name] = [fltr_a, fltr_b, fltr_c]

        # Set the maps
        self.maps[method_name] = maps

        # Set the origins
        self.origins[method_name] = origins

        # Set the methods
        # self.methods[method_name] = methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):
        return self.rgb_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the RGB maps
        self.write_maps()

        # Write origins
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

        # Plot the maps
        self.plot_maps()

# -----------------------------------------------------------------

def get_filters(string, delimiter="auto"):

    """
    This function ...
    :param string:
    :param delimiter:
    :return:
    """

    # Find the delimiter
    if delimiter == "auto": delimiter = find_delimiter(string)

    # Split and create
    return map(parse_filter, string.split(delimiter))

# -----------------------------------------------------------------
