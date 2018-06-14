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
from .component import MapMakingComponent
from ...magic.maps.oldstars.disk import DiskOldStellarMapMaker
from ...magic.maps.oldstars.bulge import BulgeOldStellarMapMaker
from ...magic.maps.oldstars.total import TotalOldStellarMapMaker
from ...magic.core.list import FrameList

# -----------------------------------------------------------------

methods = ["disk", "bulge", "total"]

# -----------------------------------------------------------------

class OldStellarMapMaker(MapMakingComponent):

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

    @property
    def disk(self):

        """
        This function ...
        :return:
        """

        return "disk" in self.config.methods

    # -----------------------------------------------------------------

    @property
    def total(self):

        """
        Thisf unction ...
        :return:
        """

        return "total" in self.config.methods

    # -----------------------------------------------------------------

    @property
    def bulge(self):

        """
        This function ...
        :return:
        """

        return "bulge" in self.config.methods

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 2. Make disk map
        if self.disk: self.make_disk_map()

        # 3. Make total map
        if self.total: self.make_total_map()

        # 4. Make bulge map
        if self.bulge: self.make_bulge_map()

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

        # Clear?
        if self.config.clear: self.clear_current_all_method(method_name)

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = DiskOldStellarMapMaker()

        # Get the frame
        nuv = self.get_frame_for_filter(self.decomposition_filter)
        frames = FrameList(nuv)

        # Get the bulge frame
        bulge = self.bulge_frame
        bulges = FrameList()
        bulges.append(bulge, self.decomposition_filter)

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

        # Clear?
        if self.config.clear: self.clear_current_all_method(method_name)

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = TotalOldStellarMapMaker()

        # Get the I1 frame
        nuv = self.get_frame_for_filter(self.decomposition_filter)
        frames = FrameList(nuv)

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

        # Clear?
        if self.config.clear: self.clear_current_all_method(method_name)

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = BulgeOldStellarMapMaker()

        # Get the bulge frame
        bulge = self.bulge_frame
        bulge.wcs = self.get_frame_for_filter(self.decomposition_filter).wcs
        bulges = FrameList()
        bulges.append(bulge, self.decomposition_filter)

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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_old()

# -----------------------------------------------------------------
