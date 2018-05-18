#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.energy.cell Contains the CellEnergyAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log

# -----------------------------------------------------------------

class CellEnergyAnalyser(AnalysisComponent):
    
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
        super(CellEnergyAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Get the dust absorptions
        self.get_absorptions()

        # Get the stellar emissions
        self.get_emissions()

        # Get the maps
        self.get_maps()

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
        super(CellEnergyAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def get_absorptions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_emissions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the stellar emission luminosities ...")

        # Old bulge
        self.get_bulge_emissions()

        # Old disk
        self.get_disk_emissions()

        # Young
        self.get_young_emissions()

        # Ionizing
        self.get_ionizing_emissions()

        # Total
        self.get_total_emissions()

    # -----------------------------------------------------------------

    def get_bulge_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_disk_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_young_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_ionizing_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_total_emissions(self):

        """
        Thisn function ...
        :return:
        """

        # Combine bulge, disk, young and ionizing

    # -----------------------------------------------------------------

    def get_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the maps of the energy imbalance ...")

        # Get the map
        if self.needs_map: self.get_map()

        # Get the map in the midplane
        if self.needs_map_midplane: self.get_map_midplane()

    # -----------------------------------------------------------------

    def get_map(self):

        """
        Thisn function ...
        :return:
        """

        # Load
        if self.has_map: self.load_map()

        # Create
        else: self.create_map()

    # -----------------------------------------------------------------

    def load_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of the energy imbalance ...")

        # Load
        self.map = Frame.from_file(self.map_path)

    # -----------------------------------------------------------------

    def create_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the energy imbalance ...")

    # -----------------------------------------------------------------

    def get_map_midplane(self):

        """
        Thisn function ...
        :return:
        """

        # Load
        if self.has_map_midplane: self.load_map_midplane()

        # Create
        else: self.create_map_midplane()

    # -----------------------------------------------------------------

    def load_map_midplane(self):

        """
        Thisn function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_map_midplane(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        if self.do_write_absorptions: self.write_absorptions()

        if self.do_write_emissions: self.write_emissions()

        if self.do_write_map: self.write_map()

        if self.do_write_map_midplane: self.write_map_midplane()

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return: 
        """

    # -----------------------------------------------------------------

    def write_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_map_midplane(self):

        """
        Thisn function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
