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

# Import the relevant PTS classes and modules
from ..basics import Configurable
from ..simulation import SkiFile, SkirtExec, SkirtArguments
from ..extract import TimeLineExtractor

# -----------------------------------------------------------------

class ResourceEstimator(Configurable):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        #super(ResourceEstimator, self).__init__(config)

        ## Attributes

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Set the simulation instance to None initially
        self.simulation = None

        # Set the ski file instance to None initially
        self.ski_file = None

        # Set the log file instance to None initially
        self.log_file = None

        # Set the arguments to None initially
        self.arguments = None

        # Set the path to the temporary directory to None initially
        self.temp_path = None

    # -----------------------------------------------------------------

    def run(self, ski_path, processes=1, threads=1):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(ski_path, processes, threads)

        # 2. Make the temporary directory
        self.make_temp(ski_path)

        # 3. Set the parameters for the simulation
        self.set_parameters()

        # 4. Run the simulation
        self.simulate()

        # 5. Remove the output
        self.clear()

    # -----------------------------------------------------------------

    def setup(self, ski_path, processes, threads):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        #super(ResourceEstimator, self).setup()

        # Adjust settings
        self.ski_file = SkiFile(ski_path)
        self.processes = processes
        self.threads = threads

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

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Create a SkirtArguments object
        self.arguments = SkirtArguments()

        # Adjust the parameters
        self.arguments.ski_pattern = self.ski_file.path
        self.arguments.parallel.processes = self.processes
        self.arguments.parallel.threads = self.threads
        self.arguments.output_path = self.temp_path
        self.arguments.emulate = True
        self.arguments.logging.verbose = False
        self.arguments.logging.memory = False
        self.arguments.logging.allocation = False
        self.arguments.single = True

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Run SKIRT in emulation mode
        self.simulation = self.skirt.run(self.arguments, silent=True)

        # Get the simulation's ski file and log file
        self.log_file = self.simulation.log_file
        self.ski_file = self.simulation.ski_file

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

        # Assume each process consumes the same amount of memory
        return self.log_file.peak_memory * self.processes

    # -----------------------------------------------------------------

    @property
    def walltime(self):

        """
        This function ...
        :return:
        """

        # Create a TimeLineExtractor instance
        extractor = TimeLineExtractor()
        extractor.run(self.simulation)

        # Determine the number of photon packages defined in the ski file
        packages = self.ski_file.packages

        # Calculate the rate of photon packages launched per second
        rate = self.log_file.packages("stellar") / extractor.stellar

        # Calculate the rate of photon pakcages launched per second without recording absorption
        rate_noabsorption = self.log_file.packages("dustem") / extractor.dustem

        # Estimate the walltime for the stellar emission phase
        stellar_walltime = self.ski_file.packages / rate

        # Estimate the walltime for the dust self-absorption phase
        factors = [1./10., 1./3., 1.]
        cycles = [3, 3, 3]  # this is a guess !!
        total_pp = (factors[0] * cycles[0] + factors[1] * cycles[1] + factors[2] * cycles[2]) * packages
        selfabs_walltime = total_pp * rate

        # Estimate the walltime for the dust emission phase
        dustem_walltime = self.ski_file.packages * self.ski_file.emission_boost / rate_noabsorption

        # Return an estimation for the walltime
        return extractor.duration_without(["stellar", "dust"]) + stellar_walltime + selfabs_walltime + dustem_walltime

# -----------------------------------------------------------------