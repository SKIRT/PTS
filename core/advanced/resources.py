#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.resources Estimate the resources needed for a particular simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..simulation.skifile import SkiFile
from ..simulation.execute import SkirtExec
from ..simulation.arguments import SkirtArguments
from ..extract.timeline import TimeLineExtractor
from ..tools import time
from ..basics.log import log
from ..tools import filesystem as fs

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

        # -- Attributes --

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Set the simulation instance to None initially
        self.simulation = None

        # Set the ski file instance to None initially
        self.ski_file = None

        # The number of processes and threads
        self.processes = None
        self.threads = None

        # Set the log file instance to None initially
        self.log_file = None

        # Set the arguments to None initially
        self.arguments = None

        # Set the path to the temporary directory to None initially
        self.temp_path = None

        # Set the timeline extractor to None initially
        self.extractor = None

    # -----------------------------------------------------------------

    def run(self, ski_path, processes=1, threads=1):

        """
        This function ...
        :param ski_path:
        :param processes:
        :param threads:
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

        # 5. Extract the timeline information
        self.extract()

        # 6. Remove the output
        self.clear()

    # -----------------------------------------------------------------

    def setup(self, ski_path, processes, threads):

        """
        This function ...
        :param ski_path:
        :param processes:
        :param threads:
        :return:
        """

        # Adjust settings
        self.ski_file = SkiFile(ski_path)
        self.processes = processes
        self.threads = threads

    # -----------------------------------------------------------------

    def make_temp(self, ski_path):

        """
        This function ...
        :param ski_path:
        :return:
        """

        # Determine the path to the temporary directory
        base_path = fs.directory_of(ski_path) if "/" in ski_path else fs.cwd()
        temp_name = time.unique_name("temp")
        self.temp_path = fs.join(base_path, temp_name)

        # Create the temporary directory if necessary
        if not fs.is_directory(self.temp_path): fs.create_directory(self.temp_path, recursive=True)

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

        # Inform the user
        log.info("Running the simulation in emulation mode to estimate the required resources ...")

        # Run SKIRT in emulation mode
        self.simulation = self.skirt.run(self.arguments, silent=True)

        # Get the simulation's ski file and log file
        self.log_file = self.simulation.log_file
        self.ski_file = self.simulation.ski_file

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Create a TimeLineExtractor instance
        self.extractor = TimeLineExtractor()
        self.extractor.run(simulation=self.simulation)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Remove the temporary directory if it exists
        if fs.is_directory(self.temp_path): fs.remove_directory(self.temp_path)

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

        walltime_emission, walltime_other = self.splitted_walltime

        # Return an estimation for the walltime
        return walltime_other + walltime_emission

    # -----------------------------------------------------------------

    @property
    def splitted_walltime(self):

        """
        This function ...
        :return:
        """

        # Determine the number of photon packages defined in the ski file
        packages = self.ski_file.packages()

        # Calculate the rate of photon packages launched per second
        rate = self.log_file.stellar_packages / self.extractor.stellar

        # Calculate the rate of photon packages launched per second without recording absorption
        rate_noabsorption = self.log_file.dust_packages / self.extractor.dustemission_photons

        # Estimate the walltime for the stellar emission phase
        stellar_walltime = packages / rate

        # Estimate the walltime for the dust self-absorption phase
        factors = [1./10., 1./3., 1.]
        cycles = [3, 3, 3]  # this is a guess !!
        total_pp = (factors[0] * cycles[0] + factors[1] * cycles[1] + factors[2] * cycles[2]) * packages
        selfabs_walltime = total_pp * rate

        # Estimate the walltime for the dust emission phase
        dustem_walltime = packages * self.ski_file.emission_boost / rate_noabsorption

        #print("stellar walltime: ", stellar_walltime)
        #print("selfabs walltime: ", selfabs_walltime)
        #print("dustem walltime: ", dustem_walltime)
        #print("without stellar and dust: ", self.extractor.duration_without(["stellar", "dust"]))

        return stellar_walltime + selfabs_walltime + dustem_walltime, self.extractor.duration_without(["stellar", "dust"])

    # -----------------------------------------------------------------

    def walltime_for(self, processes, threads):

        """
        This function ...
        :param processes:
        :param threads:
        :return:
        """

        walltime_emission, walltime_other = self.splitted_walltime

        processors = processes * threads

        # Return an estimation for the walltime
        return walltime_other + walltime_emission / processors

# -----------------------------------------------------------------
