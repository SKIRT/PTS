#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.old Contains the OldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...magic.maps.oldstars.disk import DiskOldStellarMapMaker
from ...magic.maps.oldstars.bulge import BulgeOldStellarMapMaker
from ...magic.maps.oldstars.total import TotalOldStellarMapMaker
from ...magic.core.list import FrameList

# -----------------------------------------------------------------

methods = ["disk", "bulge", "total"]

# -----------------------------------------------------------------

class OldStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(OldStellarMapMaker, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_old_path

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make disk map
        self.make_disk_map()

        # 3. Make total map
        self.make_total_map()

        # 4. Make bulge map
        self.make_bulge_map()

        # 5. Writing
        self.write()

        # 6. Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(OldStellarMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_disk_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the old stellar disk map ...")

        # Set the method name
        method_name = "disk"

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = DiskOldStellarMapMaker()


        if "IRAC I1" in self.frame_list:

            # Get the I1 frame
            i1 = self.get_frame_for_filter(self.i1_filter)
            frames = FrameList(i1)

            # Get the bulge frame
            bulge = self.bulge_frame
            bulges = FrameList(i1=bulge)

        else:

            # Get the I2 frame
            i2 = self.get_frame_for_filter(self.i2_filter)
            frames = FrameList(i2)

            # Get the bulge frame
            bulge = self.bulge_frame
            bulges = FrameList(i2=bulge)

        # Run
        maker.run(frames=frames, bulges=bulges, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_total_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the total map of the old stars ...")

        # Set the method name
        method_name = "total"

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = TotalOldStellarMapMaker()

        if "IRAC I1" in self.frame_list:

            # Get the I1 frame
            i1 = self.get_frame_for_filter(self.i1_filter)
            frames = FrameList(i1)

        else:

            # Get the I2 frame
            i2 = self.get_frame_for_filter(self.i2_filter)
            frames = FrameList(i2)

        # Run
        maker.run(frames=frames, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_bulge_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of the old stellar bulge ...")

        # Set the method name
        method_name = "bulge"

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = BulgeOldStellarMapMaker()

        if "IRAC I1" in self.frame_list:

            # Get the bulge frame
            bulge = self.bulge_frame
            bulge.wcs = self.get_frame_for_filter(self.i1_filter).wcs
            # bulge.filter = self.i1_filter
            bulges = FrameList(i1=bulge)

        else:

            # Get the bulge frame
            bulge = self.bulge_frame
            bulge.wcs = self.get_frame_for_filter(self.i2_filter).wcs
            # bulge.filter = self.i2_filter
            bulges = FrameList(i2=bulge)

        # Run
        maker.run(bulges=bulges, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    @property
    def scale(self):

        """
        Thisf unction ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps(scale=self.scale)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the NaNs masks
        self.plot_nans()

# -----------------------------------------------------------------
