#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.fitskirt Contains the FitSKIRTLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import math
import subprocess

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection, monitoring
from ...core.tools import parallelization as par
from ...core.simulation.parallelization import Parallelization

# -----------------------------------------------------------------

class FitSKIRTDefinition(object):

    """
    This function ...
    """

    def __init__(self, ski_path, fski_path, output_path, input_path=None, name=None):

        """
        The constructor ...
        :param ski_path:
        :param fski_path:
        :param output_path:
        :param input_path:
        :param name:
        """

        # Options for the ski file pattern
        self.ski_path = ski_path
        self.fski_path = fski_path

        # The input and output paths
        self.input_path = input_path
        self.output_path = output_path

        # A name for this simulation
        self.name = name

# -----------------------------------------------------------------

class FitSKIRT(object):

    """
    This function ...
    """

    def __init__(self, path=None, mpi_style="generic"):

        """
        This function ...
        :param path:
        :param mpi_style:
        """

        # Set the FitSKIRT path
        self.path = path if path is not None else ""

        if not self.path.endswith("skirt"): self.path = fs.join(self.path, "fitskirt")
        if self.path != "fitskirt": self.path = os.path.realpath(os.path.expanduser(self.path))

        if self.path == "fitskirt":
            if introspection.skirt_is_present(): self.path = introspection.fitskirt_path
            else: raise EnvironmentError("FitSKIRT is not installed or not in the PATH environment variable")

        # Indicate we are not running yet
        self._process = None

        # Set the MPI style
        self.mpi_style = mpi_style.lower()

    # -----------------------------------------------------------------

    def run(self, definition, parallelization, wait=True, silent=False):

        """
        This function ...
        :param definition:
        :param parallelization:
        :param wait:
        :param silent:
        :return:
        """

        # Create the arguments
        #arguments = FitSkirtArguments.from_definition(definition_or_arguments, logging_options, parallelization)

        # Check whether MPI is present on this system if multiple processe are requested
        if parallelization.nprocesses > 1 and not introspection.has_mpi():
            log.warning("No mpirun executable: not running")
            return []

        # Determine the MPI command
        if self.mpi_style == "lsf":
            scheduler = True
            mpi_command = "mpirun -lsf"
        elif self.mpi_style == "generic":
            scheduler = False
            mpi_command = "mpirun"
        else: raise ValueError("Invalid MPI style")

        # Get the command string
        #command = arguments.to_command(self.path, mpi_command, scheduler)

        # mkdir sample_output
        # fitskirt -t 1 tutorial.fski -o ./sample_output
        # mpirun -n 20 fitskirt -t 1 -o fit1 tutorial.fski

        # Construct command
        if parallelization.nprocesses > 1: parts = [mpi_command, "-n", str(parallelization.nprocesses), self.path]
        else: parts = [self.path]
        parts += [definition.fski_path]
        parts += ["-t", str(parallelization.nthreads), "-i", definition.input_path, "-o", definition.output_path]

        # Debugging
        command = " ".join(parts)
        log.debug("The command to launch FitSKIRT is: '" + command + "'")

        # Launch the FitSKIRT command
        if wait:
            self._process = None
            if silent: subprocess.call(parts, stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w'))
            else: subprocess.call(parts)
        else: self._process = subprocess.Popen(parts, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # Return the list of simulations so that their results can be followed up
        #simulations = arguments.simulations(simulation_names=simulation_names)

        # Check whether FitSKIRT has started
        returncode = self._process.poll() if self._process is not None else None
        if wait or returncode is not None: # when wait=True, or returncode is not None, FitSKIRT executable should have finished

            # Check presence of log file
            #if not fs.is_file(simulations.logfilepath()): raise RuntimeError("SKIRT executable has stopped but log file is not present")

            # Check presence of output files
            if fs.is_empty(definition.output_path): raise RuntimeError("FitSKIRT executable has stopped but no output present")

        # Return the list of simulations
        #return simulations

# -----------------------------------------------------------------

class FitSKIRTLauncher(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(FitSKIRTLauncher, self).__init__(config, interactive)

        # The FitSKIRT execution context
        self.fitskirt = FitSKIRT()

        # The definition
        self.definition = None

        # The parallelization scheme
        self.parallelization = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the simulation definition (if necessary)
        if self.definition is None: self.create_definition()

        # 2. Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()
        else: self.check_parallelization()

        # 3. Launch the simulation
        self.launch()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(FitSKIRTLauncher, self).setup(**kwargs)

        # Setup the remote execution context
        #if self.config.remote is not None:
        #    self.remote = SkirtRemote()
        #    self.remote.setup(self.config.remote, self.config.cluster)

        # Create output directory
        if self.config.create_output and not fs.is_directory(self.config.output): fs.create_directory(self.config.output)

    # -----------------------------------------------------------------

    def create_definition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the FitSKIRT definition ...")

        # Create the simulation definition
        self.definition = FitSKIRTDefinition(self.config.ski, self.config.fski, self.config.output, self.config.input)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        threads_per_core = par.nthreads_per_core()

        # Check whether MPI is available on this system
        if introspection.has_mpi():

            # Pure MPI
            processes = par.ncores() * threads_per_core
            threads = 1
            threads_per_core = 1

        # No MPI available
        else:

            processes = 1
            threads = min(int(math.ceil(monitoring.free_cpus())), 1)

        # Debugging
        log.debug("The number of thread per core is " + str(threads_per_core))
        log.debug("The number of processes is " + str(processes))

        # Set the parallelization scheme
        self.parallelization = Parallelization.from_processes_and_threads(processes, threads, threads_per_core=threads_per_core)

    # -----------------------------------------------------------------

    def check_parallelization(self):

        """
        This function checks whether the parallelization scheme that is asked by the user is possible given the
        number of cores and hyperthreads per core on the remote host.
        Returns:
        """

        # Inform the user
        log.info("Checking the parallelization scheme ...")

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching FitSKIRT ...")

        # Launch remotely or locally
        if self.config.remote is not None: self.launch_remote()
        else: self.launch_local()

    # -----------------------------------------------------------------

    def launch_remote(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def launch_local(self):

        """
        This function ...
        :return:
        """

        # INform the user
        log.info("Launching FitSKIRT locally ...")

        # Run FitSKIRT
        self.fitskirt.run(self.definition, silent=False, wait=True, parallelization=self.parallelization)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
