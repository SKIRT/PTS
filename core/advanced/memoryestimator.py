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
from ..advanced.dustgridtool import get_statistics
from ..tools import formatting as fmt
from ..basics.log import log
from ..simulation.memory import MemoryRequirement
from ..units.parsing import parse_unit as u
from ..tools import time

# -----------------------------------------------------------------

# A Gigabyte is 1,073,741,824 (2^30) bytes
bytes_per_gigabyte = 1073741824.

# -----------------------------------------------------------------

def estimate_memory(ski_path, input_path=None, ncells=None):

    """
    This function ...
    :param ski_path:
    :param input_path:
    :param ncells:
    :return:
    """

    # The memory estimator
    estimator = MemoryEstimator()

    # Configure the memory estimator
    estimator.config.ski = ski_path
    estimator.config.input = input_path
    estimator.config.show = False
    estimator.config.ncells = ncells

    # Estimate the memory
    estimator.run()

    # Get the memory requirement
    return estimator.memory

# -----------------------------------------------------------------

class MemoryEstimator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(MemoryEstimator, self).__init__(*args, **kwargs)

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Estimate
        self.estimate()

        # Show
        if self.config.show: self.show()

        # Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MemoryEstimator, self).setup(**kwargs)

        # Load the ski file
        self.ski = self.config.ski if isinstance(self.config.ski, SkiFile) else SkiFile(self.config.ski)

        # Path to temporary directory
        if not fs.is_directory(introspection.pts_temp_dir): fs.create_directory(introspection.pts_temp_dir)
        self.temp_path = introspection.create_temp_dir(time.unique_name("memory_estimator"))

    # -----------------------------------------------------------------

    @property
    def memory(self):

        """
        This function ...
        :return:
        """

        return MemoryRequirement(self.serial_memory, self.parallel_memory)

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the memory requirements ...")

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

        # Inform the user
        log.info("Determining the number of dust cells ...")

        # Get the number of dust cells
        if self.ski.treegrid():

            if self.config.ncells is not None: self.ncells = self.config.ncells
            elif self.config.probe: self.estimate_ncells()
            else: raise ValueError("The number of dust cells could not be determined (probing is disabled)")

        else: self.ncells = self.ski.ncells()

        # Debugging
        log.debug("The number of dust cells is " + str(self.ncells))

    # -----------------------------------------------------------------

    def get_nwavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the number of wavelength points ...")

        # Get the number of wavelengths
        if self.ski.wavelengthsfile():

            if self.config.nwavelengths is not None: self.nwavelengths = self.config.nwavelengths
            elif self.config.input is not None: self.nwavelengths = self.ski.nwavelengthsfile(self.config.input)
            else: raise ValueError("Wavelength file is used but input directory (or paths) not specified, nwavelengths also not passed in configuration")

        else: self.nwavelengths = self.ski.nwavelengths()

        # Debugging
        log.debug("The number of wavelengths is " + str(self.nwavelengths))

    # -----------------------------------------------------------------

    def get_npixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the total number of instrument pixels ...")

        self.npixels = 0
        for name, instrument_type, npixels in self.ski.npixels(self.nwavelengths):
            self.npixels += npixels

        # Debugging
        log.debug("The number of instrument pixels is " + str(self.npixels))

    # -----------------------------------------------------------------

    def get_other(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining other simulation properties ...")

        # Get the number of dust components
        self.ncomponents = self.ski.ncomponents()

        # Get the number of items in the dust library
        self.nitems = self.ski.nlibitems(ncells=self.ncells)

        # Get the number of dust populations (all dust mixes combined)
        self.npopulations = self.ski.npopulations()

        # Get other simulation properties
        self.dust_emission = self.ski.dustemission()
        self.self_absorption = self.ski.dustselfabsorption()
        self.transient_heating = self.ski.transientheating()

        # Debugging
        log.debug("The number of dust components is " + str(self.ncomponents))
        log.debug("The number of library items is " + str(self.nitems))
        log.debug("The number of dust populations is " + str(self.npopulations))
        log.debug("Dust emission is enabled" if self.dust_emission else "Dust emission is disabled")
        log.debug("Dust self-absorption is enabled" if self.self_absorption else "Dust self-absorption is disabled")
        log.debug("Transient heating is enabled" if self.transient_heating else "Transient heating is disabled")

    # -----------------------------------------------------------------

    def estimate_ncells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the number of dust cells ...")

        # Debugging
        log.debug("Running a test simulation in the temporary directory '" + self.temp_path + "' ...")

        # Get the dust grid statistics
        statistics = get_statistics(self.ski, self.temp_path, self.config.input, "test")

        # Get the number of dust cells
        self.ncells = statistics.ncells

    # -----------------------------------------------------------------

    def estimate_oligo(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the total memory usage for an oligochromatic simulation ...")

    # -----------------------------------------------------------------

    def estimate_pan(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the total memory usage for a panchromatic simulation ...")

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

        # Inform the user
        log.info("Estimating the parallel part of the memory requirement ...")

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
        self.parallel_memory = (tables_memory + instruments_memory) * u("GB")

        # TODO: determine the serial memory requirement

    # -----------------------------------------------------------------

    def estimate_serial(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the serial part of the memory requirement ...")

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

        # Set serial part of the memory requirement
        self.serial_memory = (Nbytes / bytes_per_gigabyte) * u("GB")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print(fmt.green + self.ski.prefix + fmt.reset + ":")

        treegrid = self.ski.treegrid()
        npackages = self.ski.packages()

        print("")
        print(" npackages:", npackages)
        print(" nwavelengths:", self.nwavelengths)
        print(" treegrid:", treegrid)
        print(" ncells:", self.ncells)
        print(" library items:", self.nitems)
        print(" npopulations:", self.npopulations)
        print(" selfabsorption:", self.self_absorption)
        print(" transient heating:", self.transient_heating)

        print("")

        print(" - serial part:", self.serial_memory)
        print(" - parallel", self.parallel_memory)

        print("")

        for nproc in self.config.nprocesses:

            if nproc == 1: print(fmt.underlined + str(nproc) + " process" + fmt.reset + ":")
            else: print(fmt.underlined + str(nproc) + " processes" + fmt.reset + ":")
            print("")

            print(" - serial part:", self.serial_memory)
            print(" - parallel part:", self.parallel_memory / float(nproc))
            print(" - memory per process:", self.serial_memory + self.parallel_memory / float(nproc))
            print(" - total memory (all processes):", self.serial_memory * float(nproc) + self.parallel_memory)
            #print("")
            print(" - total memory (no data parallelization):", (self.serial_memory + self.parallel_memory)*float(nproc))

            print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
