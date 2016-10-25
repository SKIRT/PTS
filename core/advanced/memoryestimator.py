#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.advanced MemoryEstimator Contains the MemoryEstimator classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.skifile import SkiFile
from ..tools import introspection
from ..tools import filesystem as fs
from ..advanced.dustgridtool import DustGridTool

# -----------------------------------------------------------------

# A Gigabyte is 1,073,741,824 (2^30) bytes
bytes_per_gigabyte = 1073741824.

# -----------------------------------------------------------------

class MemoryEstimator(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(MemoryEstimator, self).__init__(config)

        # The ski file
        self.ski = None

        # The path to a temporary directory
        self.temp_path = None

        # Properties of the simulation
        self.ncells = None
        self.nwavelengths = None
        self.npixels = None
        self.ncomponents = None
        self.nitems = None
        self.npopulations = None
        self.dust_emission = None
        self.self_absorption = None
        self.transient_heating = None

        # The calculated serial and parallel part of the memory
        self.serial_memory = None
        self.parallel_memory = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Estimate
        self.estimate()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MemoryEstimator, self).setup()

        # Load the ski file
        self.ski = self.config.ski if isinstance(self.config.ski, SkiFile) else SkiFile(self.config.ski)

        # Path to temporary directory
        if not fs.is_directory(introspection.pts_temp_dir): fs.create_directory(introspection.pts_temp_dir)
        self.temp_path = fs.create_directory_in(introspection.pts_temp_dir, "memory_estimator")

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # Get the number of dust cells
        self.get_ncells()

        # Get the number of wavelengths
        self.get_nwavelengths()

        # Get the number of instrument pixels
        self.get_npixels()

        # Get other properties
        self.get_other()

        # Estimate depending on the type of simulation
        if self.ski.oligochromatic(): self.estimate_oligo()
        elif self.ski.panchromatic(): self.estimate_pan()
        else: raise ValueError("Invalid ski file")

    # -----------------------------------------------------------------

    def get_ncells(self):

        """
        This function ...
        :return:
        """

        # Get the number of dust cells
        if self.ski.treegrid():
            if self.config.ncells is not None:
                self.ncells = self.config.ncells
            elif self.config.probe: self.estimate_ncells()
            else: raise ValueError("The number of dust cells could not be determined (probing is disabled)")
        else: self.ncells = self.ski.ncells()

    # -----------------------------------------------------------------

    def get_nwavelengths(self):

        """
        This function ...
        :return:
        """

        # Get the number of wavelengths
        if self.ski.wavelengthsfile():
            if self.config.nwavelengths is not None:
                self.nwavelengths = self.config.nwavelengths
            elif self.config.input is not None: self.nwavelengths = self.ski.nwavelengthsfile(self.config.input)
            else: raise ValueError("Wavelength file is used but input directory (or paths) not specified, nwavelengths also not passed in configuration")
        else: self.nwavelengths = self.ski.nwavelengths()

    # -----------------------------------------------------------------

    def get_npixels(self):

        """
        This function ...
        :return:
        """

        self.npixels = 0
        for name, instrument_type, npixels in self.ski.npixels(self.nwavelengths):
            self.npixels += npixels

    # -----------------------------------------------------------------

    def get_other(self):

        """
        This function ...
        :return:
        """

        # Get the number of dust components
        self.ncomponents = self.ski.ncomponents()

        # Get the number of items in the dust library
        self.nitems = self.ski.nlibitems()

        # Get the number of dust populations (all dust mixes combined)
        self.npopulations = self.ski.npopulations()

        # Get other simulation properties
        self.dust_emission = self.ski.dustemission()
        self.self_absorption = self.ski.dustselfabsorption()
        self.transient_heating = self.ski.transientheating()

    # -----------------------------------------------------------------

    def estimate_ncells(self):

        """
        This function ...
        :return:
        """

        # Create the dust grid tool
        tool = DustGridTool()

        # Get the dust grid statistics
        statistics = tool.get_statistics(self.ski, self.temp_path, self.config.input, "test")

        # Get the number of dust cells
        self.ncells = statistics.ncells

    # -----------------------------------------------------------------

    def estimate_oligo(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def estimate_pan(self):

        """
        This function ...
        :return:
        """

        # Estimate the parallel part
        self.estimate_parallel()

        # Estimate the serial part
        self.estimate_serial()

    # -----------------------------------------------------------------

    def estimate_parallel(self):

        """
        This function ...
        :return:
        """

        # Size of the parallel tables (in number of values)
        #table_size = self.nwavelengths * self.ncells

        tables_bytes = 0

        if self.dust_emission:
            tables_bytes += self.nwavelengths * self.ncells
            if self.self_absorption:
                tables_bytes += self.nwavelengths * self.ncells

        if self.dust_emission:
            if self.ncomponents > 1: tables_bytes += self.ncells + self.nwavelengths * self.ncells
            else: tables_bytes += self.ncells + self.nwavelengths * self.nitems

        # Size of the instruments
        instruments_size = self.npixels

        # Memory requirements for the parallel tables (in GB)
        #tables_memory = 8 * 3 * table_size / bytes_per_gigabyte

        tables_memory = 8 * tables_bytes / bytes_per_gigabyte

        # Memory requirements for the instruments (in GB)
        instruments_memory = 8 * instruments_size / bytes_per_gigabyte

        # Determine the parallel memory requirement
        self.parallel_memory = tables_memory + instruments_memory

        # TODO: determine the serial memory requirement

    # -----------------------------------------------------------------

    def estimate_serial(self):

        """
        This function ...
        :return:
        """

        # Overhead
        Ndoubles = 50e6 + (self.nwavelengths + self.ncells + self.ncomponents + self.npopulations) * 10

        # Dust system # serial part
        Ndoubles += (self.ncomponents + 1) * self.ncells

        # DUST SYSTEM # parallel part
        #if dustEmission:
        #    Ndoubles += Nlambda * Ncells
        #    if selfAbsorption:
        #        Ndoubles += Nlambda * Ncells

        # Dust grid tree
        if self.ski.treegrid():
            Ndoubles += self.ncells * 1.2 * 50

        # Dust mixes
        Ndoubles += self.nwavelengths * self.npopulations * 3
        Ndoubles += (self.nwavelengths + self.npopulations) * 10

        # Dust library PARALLEL
        #if dustEmission:
        #    Ndoubles += Ncells + Nlambda * Nitems

        # Transient heating
        if self.dust_emission and self.transient_heating:

            NT = 2250.
            Ndoubles += self.ncomponents * (self.nwavelengths + 1) * NT
            Ndoubles += self.npopulations * 5. / 8. * NT * NT

        Nbytes = Ndoubles * 8

        self.serial_memory = Nbytes / bytes_per_gigabyte

# -----------------------------------------------------------------
