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
from ...core.tools.logging import log
from .component import MapsComponent
from ...magic.maps.oldstars.disk import DiskOldStellarMapMaker
from ...magic.maps.oldstars.bulge import BulgeOldStellarMapMaker
from ...magic.maps.oldstars.total import TotalOldStellarMapMaker
from ...magic.core.list import FrameList

# -----------------------------------------------------------------

class OldStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(OldStellarMapMaker, self).__init__(config, interactive)

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

        # Create the maker
        maker = DiskOldStellarMapMaker()

        # Get the I1 frame
        i1 = self.get_frame_for_filter(self.i1_filter)
        frames = FrameList(i1)

        # Get the bulge frame
        bulge = self.bulge_frame
        bulges = FrameList(i1=bulge)

        # Run
        maker.run(frames=frames, bulges=bulges)

        # Set the maps
        self.maps["disk"] = maker.maps

        # Set the origins
        self.origins["disk"] = maker.origins

    # -----------------------------------------------------------------

    def make_total_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the total map of the old stars ...")

        # Create the maker
        maker = TotalOldStellarMapMaker()

        # Get the I1 frame
        i1 = self.get_frame_for_filter(self.i1_filter)
        frames = FrameList(i1)

        # Run
        maker.run(frames=frames)

        # Set the maps
        self.maps["total"] = maker.maps

        # Set the origins
        self.origins["total"] = maker.origins

    # -----------------------------------------------------------------

    def make_bulge_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of the old stellar bulge ...")

        # Create the maker
        maker = BulgeOldStellarMapMaker()

        # Get the bulge frame
        bulge = self.bulge_frame
        bulges = FrameList(i1=bulge)

        # Run
        maker.run(bulges=bulges)

        # Set the maps
        self.maps["bulge"] = maker.maps

        # Set the origins
        self.origins["bulge"] = maker.origins

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

# -----------------------------------------------------------------
