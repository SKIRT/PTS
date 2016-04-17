#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.parameterexploration Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem, time, tables
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.filter import Filter
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.launch.options import SchedulingOptions
from ...core.launch.parallelization import Parallelization

# -----------------------------------------------------------------

class ParameterExplorer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ParameterExplorer, self).__init__(config)

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski file
        self.ski = None

        # The parameter ranges
        self.young_luminosities = None
        self.ionizing_luminosities = None
        self.dust_masses = None

        # The table with the parameter values for each simulation
        self.table = None

        # A dictionary with the average runtime for the different remote hosts
        self.runtimes = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ParameterExplorer instance
        explorer = cls(arguments.config)

        # Set the modeling path
        explorer.config.path = arguments.path

        # Set options for the young stellar population
        if arguments.young_nvalues is not None: explorer.config.young_stars.nvalues = arguments.young_nvalues
        if arguments.young_range is not None:
            explorer.config.young_stars.rel_min = arguments.young_range[0]
            explorer.config.young_stars.rel_max = arguments.young_range[1]
        if arguments.young_log: explorer.config.young_stars.scale = "log"
        else: explorer.config.young_stars.scale = "linear"

        # Set options for the ionizing stellar population
        if arguments.ionizing_nvalues is not None: explorer.config.ionizing_stars.nvalues = arguments.ionizing_nvalues
        if arguments.ionizing_range is not None:
            explorer.config.ionizing_stars.rel_min = arguments.ionizing_range[0]
            explorer.config.ionizing_stars.rel_max = arguments.ionizing_range[1]
        if arguments.ionizing_log: explorer.config.ionizing_stars = "log"
        else: explorer.config.ionizing_stars.scale = "linear"

        # Set options for the dust component
        if arguments.dust_nvalues is not None: explorer.config.dust.nvalues = arguments.dust_nvalues
        if arguments.dust_range is not None:
            explorer.config.dust.rel_min = arguments.dust_range[0]
            explorer.config.dust.rel_max = arguments.dust_range[1]
        if arguments.dust_log: explorer.config.dust.scale = "log"
        else: explorer.config.dust.scale = "linear"

        # Return the new instance
        return explorer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the current parameter table
        self.load_table()

        # 3. Load the ski file
        self.load_ski()

        # 4. Set the ranges of the different fit parameters
        self.set_parameter_ranges()

        # Set the parallelization
        self.set_parallelization()

        # Estimate the runtimes
        self.estimate_runtimes()

        # 6. Launch the simulations for different parameter values
        self.simulate()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExplorer, self).setup()

        # Get the names of the filters for which we have photometry
        filter_names = []
        fluxes_table_path = filesystem.join(self.phot_path, "fluxes.dat")
        fluxes_table = tables.from_file(fluxes_table_path)
        # Loop over the entries in the fluxes table, get the filter
        for entry in fluxes_table:
            # Get the filter
            filter_id = entry["Instrument"] + "." + entry["Band"]
            filter_names.append(filter_id)

        # Set options for the BatchLauncher
        self.launcher.config.analysis.extraction.path = self.fit_res_path
        self.launcher.config.analysis.misc.path = self.fit_res_path # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.plotting.path = self.fit_plot_path # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.extraction.timeline = True # extract the simulation timeline
        self.launcher.config.analysis.plotting.seds = True  # Plot the output SEDs
        self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        self.launcher.config.analysis.misc.images = True  # Make observed images
        self.launcher.config.analysis.misc.observation_filters = filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.plotting.format = "png" # plot in PNG format so that an animation can be made from the fit SEDs
        self.launcher.config.shared_input = True   # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = True # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        #self.launcher.config.remotes = ["nancy"]   # temporary; only use Nancy

    # -----------------------------------------------------------------

    def load_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the parameter table ...")

        # Load the parameter table
        self.table = tables.from_file(self.parameter_table_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Open the ski file (created by InputInitializer)
        self.ski = SkiFile(self.fit_ski_path)

    # -----------------------------------------------------------------

    def set_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parameter ranges ...")

        # Get the current values in the ski file prepared by InputInitializer
        young_luminosity, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
        ionizing_luminosity, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
        dust_mass = self.ski.get_dust_component_mass(0)

        # Set the parameter ranges
        self.set_young_luminosity_range(young_luminosity)
        self.set_ionizing_luminosity_range(ionizing_luminosity)
        self.set_dust_mass_range(dust_mass)

    # -----------------------------------------------------------------

    def set_young_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the FUV luminosity of the young stars ...")

        # Set the range of the FUV luminosity of the young stellar population
        min = self.config.young_stars.rel_min * luminosity
        max = self.config.young_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.young_stars.scale == "linear":
            self.young_luminosities = np.linspace(min, max, num=self.config.young_stars.nvalues, endpoint=True)
        elif self.config.young_stars.scale == "logarithmic":
            self.young_luminosities = np.logspace(min, max, num=self.config.young_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the young stellar luminosity values")

    # -----------------------------------------------------------------

    def set_ionizing_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the FUV luminosity of the ionizing stars ...")

        # Determine the minimum and maximum luminosity
        min = self.config.ionizing_stars.rel_min * luminosity
        max = self.config.ionizing_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.ionizing_stars.scale == "linear":
            self.ionizing_luminosities = np.linspace(min, max, num=self.config.ionizing_stars.nvalues, endpoint=True)
        elif self.config.ionizing_stars.scale == "log":
            self.ionizing_luminosities = np.logspace(min, max, num=self.config.ionizing_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the ionizing stellar luminosity values")

    # -----------------------------------------------------------------

    def set_dust_mass_range(self, mass):

        """
        This function ...
        :param mass:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the dust mass ...")

        # Set the dust mass range
        min = self.config.dust.rel_min * mass
        max = self.config.dust.rel_max * mass

        # Create a linear or logarithmic range of dust masses
        if self.config.dust.scale == "linear":
            self.dust_masses = np.linspace(min, max, num=self.config.dust.nvalues, endpoint=True)
        elif self.config.dust.scale == "log":
            self.dust_masses = np.logspace(min, max, num=self.config.dust.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the dust mass values")

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the correponding system).
        :return:
        """

        # Loop over the IDs of the hosts used by the batch launcher that use a scheduling system
        for host in self.launcher.scheduler_hosts:

            # Get the number of cores per node for this host
            cores_per_node = host.clusters[host.cluster_name].cores

            # Determine the number of cores corresponding to 4 full nodes
            cores = cores_per_node * 4

            # Use 1 core for each process (assume there is enough memory)
            processes = cores

            # Determine the number of threads per core
            if host.use_hyperthreading: threads_per_core = host.clusters[host.cluster_name].threads_per_core
            else: threads_per_core = 1

            # Create a Parallelization instance
            parallelization = Parallelization()

            # Set the parallelization properties
            parallelization.cores = cores
            parallelization.threads_per_core = threads_per_core
            parallelization.processes = processes

            # Set the parallelization for this host
            self.launcher.set_parallelization_for_host(host.id, parallelization)

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtimes based on the results of previous runs ...")

        # Get the number of photon packages (per wavelength) for this batch of simulations
        current_packages = self.ski.packages()

        # Load the runtime table
        runtimes_table = tables.from_file(self.runtime_table_path)

        # Keep a list of all the runtimes recorded for a certain remote host
        runtimes_for_hosts = defaultdict(list)

        # Loop over the entries in the runtime table
        # "Simulation name", "Remote host", "Cores", "Hyperthreads per core", "Processes", "Packages", "Runtime"
        for i in range(len(runtimes_table)):

            host_id = runtimes_table["Host id"][i]
            #cluster_name = runtimes_table["Cluster name"][i]
            cores = runtimes_table["Cores"][i]
            threads_per_core = runtimes_table["Hyperthreads per core"][i]
            processes = runtimes_table["Processes"][i]
            packages = runtimes_table["Packages"][i]
            runtime = runtimes_table["Runtime"][i]

            parallelization = Parallelization()
            parallelization.cores = cores
            parallelization.threads_per_core = threads_per_core
            parallelization.processes = processes

            parallelization_for_host = self.launcher.parallelization_for_host(host_id)

            if parallelization == parallelization_for_host and packages == current_packages:
                runtimes_for_hosts[host_id].append(runtime)

        # For each remote host, determine the mean runtime
        for host_id in runtimes_for_hosts: self.runtimes[host_id] = np.mean(np.array(runtimes_for_hosts[host_id]))

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Create a dictionary for the scheduling options for the different remotes with a scheduling system
        sched_options = dict()

        # Create scheduling options for each remote host with a scheduling system
        for host in self.launcher.scheduler_hosts:

            # Initialize a scheduling options object
            scheduling_options = SchedulingOptions()

            # A simulation with the current number of photon packages and with the same parallelization has already
            # been executed (and analysed) on this remote host
            if host.id in self.runtimes: scheduling_options.walltime = self.runtimes[host.id] * 1.2

            # No matching simulation could be found on this remote host, estimate the walltime based on simulations
            # with other number of photon packages and/or other parallelization on this remote host, or else simulations
            # performed on other hosts
            else:

                # 1. Search for timeline files from simulations ran on this remote host

                # 2. Search for the average runtime with the current number of photon packages and parallelization
                #    on other remote hosts

                # 3. Search for timeline files from simulations ran on another remote host

                pass

            # Add the scheduling options to the dictionary
            sched_options[host.id] = scheduling_options

        # Create a FUV filter object
        fuv = Filter.from_string("FUV")

        # Loop over the different values of the young stellar luminosity
        for young_luminosity in self.young_luminosities:

            # Loop over the different values of the ionizing stellar luminosity
            for ionizing_luminosity in self.ionizing_luminosities:

                # Loop over the different values of the dust mass
                for dust_mass in self.dust_masses:

                    # Create a unique name for this combination of parameter values
                    simulation_name = time.unique_name()

                    # Change the parameter values in the ski file
                    self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv)
                    self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv)
                    self.ski.set_dust_component_mass(0, dust_mass)

                    # Determine the directory for this simulation
                    simulation_path = filesystem.join(self.fit_out_path, simulation_name)

                    # Create the simulation directory
                    filesystem.create_directory(simulation_path)

                    # Create an 'out' directory within the simulation directory
                    output_path = filesystem.join(simulation_path, "out")
                    filesystem.create_directory(output_path)

                    # Put the ski file with adjusted parameters into the simulation directory
                    ski_path = filesystem.join(simulation_path, self.galaxy_name + ".ski")
                    self.ski.saveto(ski_path)

                    # Create the SKIRT arguments object
                    arguments = create_arguments(ski_path, self.fit_in_path, output_path)

                    # Debugging
                    log.debug("Adding a simulation to the queue with:")
                    log.debug(" - ski path: " + arguments.ski_pattern)
                    log.debug(" - output path: " + arguments.output_path)

                    # Put the parameters in the queue and get the simulation object
                    self.launcher.add_to_queue(arguments, simulation_name)

                    # Set scheduling options (for the different remote hosts with a scheduling system)
                    for host_id in sched_options: self.launcher.set_scheduling_options(host_id, simulation_name, sched_options[host_id])

                    # Add an entry to the parameter table
                    self.table.add_row([simulation_name, young_luminosity, ionizing_luminosity, dust_mass])

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Set the path to the reference SED (for plotting the simulated SED against the reference points)
            simulation.analysis.plotting.reference_sed = filesystem.join(self.phot_path, "fluxes.dat")

            # Save the simulation object
            simulation.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Create and write a table with the parameter values for each simulation
        self.write_parameter_table()

    # -----------------------------------------------------------------

    def write_parameter_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter table ...")

        # Set the units of the parameter table
        self.table["FUV young"] = "Lsun_FUV"
        self.table["FUV ionizing"] = "Lsun_FUV"
        self.table["Dust mass"] = "Msun"

        # Write the parameter table
        tables.write(self.table, self.parameter_table_path, format="ascii.ecsv")

# -----------------------------------------------------------------

def create_arguments(ski_path, input_path, output_path):

    """
    This function ...
    :param ski_path:
    :param input_path:
    :param output_path:
    :return:
    """

    # Create a new SkirtArguments object
    arguments = SkirtArguments()

    # The ski file pattern
    arguments.ski_pattern = ski_path
    arguments.recursive = False
    arguments.relative = False

    # Input and output
    arguments.input_path = input_path
    arguments.output_path = output_path

    # Parallelization settings
    arguments.parallel.threads = None
    arguments.parallel.processes = None

    return arguments

# -----------------------------------------------------------------
