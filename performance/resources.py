#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package performance.resources Estimate the resources needed for a particular simulation

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import shutil

# Import the relevant PTS modules
from ..core.skirtexec import SkirtExec
from ..core.parameters import SkirtParameters

# -----------------------------------------------------------------

class ResourceEstimator(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Configuration

        # ...

        ## Attributes

        # Create the SKIRT execution context
        self.skirt = SkirtExec("~/Development/SKIRT/release/SKIRTmain/skirt")

        # Set the log file instance to None initially
        self.log_file = None

        # Set the parameters to None initially
        self.parameters = None

        # Set the path to the temporary directory to None initially
        self.temp_path = None

    # -----------------------------------------------------------------

    def run(self, ski_path, processes=1, threads=1):

        """
        This function ...
        :return:
        """

        # Make the temporary directory
        self.make_temp(ski_path)

        # Set the parameters
        self.set_parameters(ski_path, processes, threads)

        # Run the simulation
        self.simulate()

        # Remove the output
        self.clear()

    # -----------------------------------------------------------------

    def make_temp(self, ski_path):

        """
        This function ...
        :return:
        """

        # Determine the path to the temporary directory
        base_path = os.path.dirname(ski_path) if "/" in ski_path else os.getcwd()
        self.temp_path = os.path.join(base_path, "temp")

        # Create the temporary directory if necessary
        if not os.path.exists(self.temp_path): os.makedirs(self.temp_path)

    # -----------------------------------------------------------------

    def set_parameters(self, ski_path, processes, threads):

        """
        This function ...
        :return:
        """

        # Create a SkirtParameters object
        self.parameters = SkirtParameters()

        # Adjust the parameters
        self.parameters.ski_pattern = ski_path
        self.parameters.parallel.processes = processes
        self.parameters.parallel.threads = threads
        self.parameters.output_path = self.temp_path
        self.parameters.emulate = True
        self.parameters.logging.verbose = False
        self.parameters.logging.memory = False
        self.parameters.logging.allocation = False
        self.parameters.single = True

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Run SKIRT in emulation mode
        simulation = self.skirt.run(self.parameters, silent=True)

        # Get the simulation's log file
        self.log_file = simulation.logfile()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Remove the temporary directory if it exists
        if os.path.exists(self.temp_path): shutil.rmtree(self.temp_path)

    # -----------------------------------------------------------------

    @property
    def memory(self):

        """
        This function ...
        :return:
        """

        return self.log_file.peak_memory

    # -----------------------------------------------------------------

    @property
    def walltime(self):

        """
        This function ...
        :return:
        """

        return None

# -----------------------------------------------------------------

## Return the estimated memory usage for a given ski file (in gigabytes)
def estimate_memory(skifile, inputpath=None):

    # Get the number of wavelengths (force all calculations to floating point, avoiding integer overrun)
    Nlambda = None
    try:

        Nlambda = float(skifile.nwavelengths())

    except ValueError:

        # Get the name of the wavelength data file
        filename = skifile.wavelengthsfile()

        # Open the wavelengths file
        filepath = os.path.join(inputpath, filename)
        with open(filepath) as file:

            first = file.readlines()[0]

        Nlambda = float(first.split(" ")[0])

    # Get the number of dust cells
    try:
        Ncells = skifile.ncells()
    except ValueError:

        raise ValueError("Disabled for now")
        #Ncells = float(raw_input("\033[91m" + "\033[1m" + "?  How many dust cells do you expect for this simulation? ... \033[0m"))

    # Get the number of dust components
    Ncomps = skifile.ncomponents()

    # Get the number of items in the dust library
    Nitems = skifile.nlibitems()

    # Get the number of dust populations (all dust mixes combined)
    Npops = skifile.npopulations()

    # Get other simulation properties
    dustEmission = skifile.dustemission()
    selfAbsorption = skifile.dustselfabsorption()
    transientHeating = skifile.transientheating()

    # Overhead
    Ndoubles = 50e6 + (Nlambda + Ncells + Ncomps + Npops) * 10

    # Instruments
    for instrument in skifile.npixels(Nlambda):

        Ndoubles += instrument[2]

    # Dust system
    Ndoubles += (Ncomps+1) * Ncells
    if dustEmission:
        Ndoubles += Nlambda * Ncells
        if selfAbsorption:
            Ndoubles += Nlambda * Ncells

    # Dust grid tree
    if True:

        Ndoubles += Ncells * 1.2 * 50

    # Dust mixes
    Ndoubles += Nlambda * Npops * 3
    Ndoubles += (Nlambda + Npops) * 10

    # Dust library
    if dustEmission:
        Ndoubles += Ncells + Nlambda*Nitems

    # Transient heating
    if dustEmission and transientHeating:

        NT = 2250.
        Ndoubles += Ncomps * (Nlambda+1) * NT
        Ndoubles += Npops * 5./8.*NT*NT

    # 8 bytes in a double
    Nbytes = Ndoubles * 8

    # Return the number of gigabytes
    return Nbytes/1e9

# -----------------------------------------------------------------