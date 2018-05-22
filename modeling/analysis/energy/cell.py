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

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Load the run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

    # -----------------------------------------------------------------

    def get_absorptions(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_absorptions: self.load_absorptions()

        # Create
        else: self.create_absorptions()

    # -----------------------------------------------------------------

    def load_absorptions(self):

        """
        This function ...
        :return:
        """



    # -----------------------------------------------------------------

    def create_absorptions(self):

        """
        This function ...
        :return:
        """

        total_absorptions = self.total_contribution_absorption_luminosities
        total_unit = self.total_contribution_absorption_unit
        total_conversion = total_unit.conversion_factor("W")
        total_absorptions_watt = total_absorptions * total_conversion

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

        if self.has_bulge_emissions: self.load_bulge_emissions()

        else: self.create_bulge_emissions()

    # -----------------------------------------------------------------

    def load_bulge_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_bulge_emissions(self):

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

        if self.has_disk_emissions: self.load_disk_emissions()

        else: self.create_disk_emissions()

    # -----------------------------------------------------------------

    def load_disk_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_disk_emissions(self):

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

        # Load
        if self.has_young_emissions: self.load_young_emissions()

        # Create
        else: self.create_young_emissions()

    # -----------------------------------------------------------------

    def load_young_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_young_emissions(self):

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

        # Load
        if self.has_ionizing_emissions: self.load_ionizing_emissions()

        # Create
        else: self.create_ionizing_emissions()

    # -----------------------------------------------------------------

    def load_ionizing_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_ionizing_emissions(self):

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
        if self.has_total_emissions: self.load_total_emissions()

        # Create
        else: self.create_total_emissions()

    # -----------------------------------------------------------------

    def load_total_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_total_emissions(self):

        """
        Thisn function ...
        :return:
        """

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

        # Inform the user
        log.info("")

    # -----------------------------------------------------------------

    def create_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the energy imbalance in the midplane ...")

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_disk_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_young_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_ionizing_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_total_emissions(self):

        """
        Thisn function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_map_midplane(self):

        """
        Thisn function ...
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

        # Absorptions
        if self.do_write_absorptions: self.write_absorptions()

        # Emissions
        self.write_emissions()

        # Map
        if self.do_write_map: self.write_map()

        # Map of midplane
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

        # Bulge
        if self.do_write_bulge_emissions: self.write_bulge_emissions()

        # Disk
        if self.do_write_disk_emissions: self.write_disk_emissions()

        # Young
        if self.do_write_young_emissions: self.write_young_emissions()

        # Ionizing
        if self.do_write_ionizing_emissions: self.write_ionizing_emissions()

        # Total
        if self.do_write_total_emissions: self.write_total_emissions()

    # -----------------------------------------------------------------

    def write_bulge_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_disk_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_young_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_ionizing_emissions(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_total_emissions(self):

        """
        Thisn function ...
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
