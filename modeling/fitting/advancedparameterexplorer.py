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
from collections import defaultdict

# Import the relevant PTS classes and modules
from .parameterexploration import ParameterExplorer
from ...core.tools import filesystem, tables
from ...core.tools.logging import log
from ...core.launch.options import SchedulingOptions
from ...core.launch.parallelization import Parallelization
from ...core.basics.distribution import Distribution
from ...core.extract.timeline import TimeLineExtractor

# -----------------------------------------------------------------

class AdvancedParameterExplorer(ParameterExplorer):
    
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
        super(AdvancedParameterExplorer, self).__init__(config)

        # A dictionary with the scheduling options for the different remote hosts
        self.scheduling_options = dict()

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

        # Set options ...

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

        # 4. Set the combinations of parameter values
        self.set_parameters()

        # 5. Set the parallelization schemes for the different remote hosts
        self.set_parallelization()

        # 6. Estimate the runtimes for the different remote hosts
        self.estimate_runtimes()

        # 7. Launch the simulations for different parameter values
        self.simulate()

        # 8. Writing
        self.write()

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        pass

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
        #self.set_young_luminosity_range(young_luminosity)
        #self.set_ionizing_luminosity_range(ionizing_luminosity)
        #self.set_dust_mass_range(dust_mass)

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
            parallelization = Parallelization(cores, threads_per_core, processes)

            # Set the parallelization for this host
            self.launcher.set_parallelization_for_host(host.id, parallelization)

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        # Get the number of photon packages (per wavelength) for this batch of simulations
        current_packages = self.ski.packages()

        # A dictionary with the average runtime for the different remote hosts
        runtimes = dict()

        # Inform the user
        log.info("Estimating the runtimes based on the results of previous runs ...")

        # Debugging
        log.debug("Loading the table with the total runtimes of previous simulations ...")

        # Load the runtime table
        runtimes_table = tables.from_file(self.runtime_table_path, format="ascii.ecsv")

        # Get lists of the runtimes for each host, with the current configuration (packages, parallelization)
        runtimes_for_hosts = self.get_runtimes_hosts(runtimes_table)

        log.debug("Runtimes of previous simulations run with the same configuration (packages, parallelization) were found for the following remote hosts: '" + "', '".join(runtimes_for_hosts.keys()) + "'")

        # For each remote host, determine the most frequent runtime (runtimes are binned, most frequent value is center of most occupied bin)
        for host_id in runtimes_for_hosts:

            # Debugging
            parallelization_for_host = self.launcher.parallelization_for_host(host_id)
            log.debug("Determining the most frequent runtime for remote host '" + host_id + "' for the current "
                      "configuration of " + str(current_packages) + " photon packages and a parallelization scheme with"
                      + str(parallelization_for_host.cores) + " cores, " + str(parallelization_for_host.threads_per_core)
                      + " threads per core and " + str(parallelization_for_host.processes) + " processes ...")

            distribution = Distribution.from_values(runtimes_for_hosts[host_id], bins=25)
            runtimes[host_id] = distribution.most_frequent

            if log.is_debug(): distribution.plot(title="Distribution of runtimes for remote host '" + host_id + "'")

        # Overestimation factor
        overestimation = 1.2

        # Debugging
        log.debug("The factor for determining an overestimation of the runtime (for incorporating random effects) is " + str(overestimation))

        walltimes = dict()
        scheduling_hosts_without_runtime = []

        # Loop over the hosts that use a scheduling system to see whether we have a record of the total runtime for its
        # current configuration (packages, parallelization)
        for host in self.launcher.scheduler_hosts:

            if host.id in runtimes:
                walltimes[host.id] = runtimes[host.id] * overestimation
                log.debug("The walltime used for remote '" + host.id + "' is " + str(walltimes[host.id]) + " seconds")
            else: scheduling_hosts_without_runtime.append(host.id)

        # Debugging
        log.debug("No runtimes were found for the same configuration (packages, parallelization) for the following remote hosts (with scheduling system): '" + "', '".join(scheduling_hosts_without_runtime) + "'")

        # Remote hosts with scheduling system for which no runtimes could be found for the current configuration (packages, parallelization)
        if len(scheduling_hosts_without_runtime) > 0:

            serial_parallel_overhead_for_hosts = defaultdict(list)

            # Indices of the simulations in the runtime table with the current number of photon packages
            indices = tables.find_indices(runtimes_table, current_packages, "Packages")

            # Debugging
            log.debug(str(len(indices)) + " simulations were found that had the same number of photon packages as the current configuration")

            # Loop over the matching indices
            for index in indices:

                # The simulation name
                simulation_name = runtimes_table["Simulation name"][index]

                # The ID of the remote host on which the simulation was run
                host_id = runtimes_table["Host id"][index]

                # Get the parallelization properties for this particular simulation
                cores = runtimes_table["Cores"][index]
                threads_per_core = runtimes_table["Hyperthreads per core"][index]
                processes = runtimes_table["Processes"][index]

                # Determine the path to the timeline table file
                timeline_path = filesystem.join(self.fit_res_path, simulation_name, "timeline.dat")
                if not filesystem.is_file(timeline_path):
                    log.warning("The timeline table file does not exist for simulation '" + simulation_name + "'")
                    continue

                # Get the serial time, parallel time and overhead for this particular simulation
                extractor = TimeLineExtractor.open_table(timeline_path)
                serial = extractor.serial
                parallel = extractor.parallel
                overhead = extractor.overhead

                parallel_times_cores = parallel * cores
                overhead_per_core = overhead / cores

                serial_parallel_overhead_for_hosts[host_id].append((serial, parallel_times_cores, overhead_per_core))

            # Debugging
            log.debug("Found timeline information for simulations with the same number of photon packages as the current configuration for these remote hosts: '" + "', '".join(serial_parallel_overhead_for_hosts.keys()) + "'")

            # Loop over the hosts with scheduling system for which a runtime reference was not found for the current configuration (packages, parallelization)
            for host_id in scheduling_hosts_without_runtime:

                # Get the parallelization scheme that we have defined for this particular host
                parallelization_for_host = self.launcher.parallelization_for_host(host_id)
                cores_for_host = parallelization_for_host.cores

                if host_id in serial_parallel_overhead_for_hosts:

                    # Debugging
                    log.debug("Timeline information for remote host '" + host_id + "' was found")

                    runtimes = []
                    for serial, parallel_times_cores, overhead_per_core in serial_parallel_overhead_for_hosts[host_id]:
                        runtimes.append(serial + parallel_times_cores / cores_for_host + overhead_per_core * cores_for_host)

                    distribution = Distribution.from_values(runtimes, bins=25)

                    distribution.plot()

                    runtime = distribution.most_frequent
                    walltimes[host_id] = runtime * overestimation

                    # Debugging
                    log.debug("The walltime used for remote host '" + host_id + "' (parallelization scheme with " + str(cores_for_host) + " cores) is " + str(walltimes[host_id]) + " seconds")

                else:

                    # Debugging
                    log.debug("Timeline information for remote host '" + host_id + "' was not found, using information from other remote hosts ...")

                    runtimes = []
                    for other_host_id in serial_parallel_overhead_for_hosts:
                        for serial, parallel_times_cores, overhead_per_core in serial_parallel_overhead_for_hosts[other_host_id]:
                            runtimes.append(serial + parallel_times_cores / cores_for_host + overhead_per_core * cores_for_host)

                    distribution = Distribution.from_values(runtimes, bins=25)

                    if log.is_debug(): distribution.plot(title="Distribution of runtimes estimated for remote host '" + host_id + "' based on simulations run on other hosts")

                    runtime = distribution.most_frequent
                    walltimes[host_id] = runtime * overestimation

                    # Debugging
                    log.debug("The walltime used for remote host '" + host_id + "' (parallelization scheme with " + str(cores_for_host) + " cores) is " + str(walltimes[host_id]) + " seconds")

        # Create scheduling options
        for host_id in walltimes:
            self.scheduling_options[host_id] = SchedulingOptions()
            self.scheduling_options[host_id].walltime = walltimes[host_id]

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
        self.table["FUV young"].unit = "Lsun_FUV"
        self.table["FUV ionizing"].unit = "Lsun_FUV"
        self.table["Dust mass"].unit = "Msun"

        # Write the parameter table
        tables.write(self.table, self.parameter_table_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def get_runtimes_hosts(self, runtimes_table):

        """
        This function ...
        :param runtimes_table:
        :return:
        """

        # Get the number of photon packages (per wavelength) for this batch of simulations
        current_packages = self.ski.packages()

        # Keep a list of all the runtimes recorded for a certain remote host
        runtimes_for_hosts = defaultdict(list)

        # Loop over the entries in the runtime table
        # "Simulation name", "Host id", "Cluster name", "Cores", "Hyperthreads per core", "Processes", "Packages", "Runtime"
        for i in range(len(runtimes_table)):

            # Get the ID of the host and the cluster name for this particular simulation
            host_id = runtimes_table["Host id"][i]
            cluster_name = runtimes_table["Cluster name"][i]

            # Get the parallelization properties for this particular simulation
            cores = runtimes_table["Cores"][i]
            threads_per_core = runtimes_table["Hyperthreads per core"][i]
            processes = runtimes_table["Processes"][i]

            # Get the number of photon packages (per wavelength) used for this simulation
            packages_simulation = runtimes_table["Packages"][i]

            # Get the total runtime
            runtime = runtimes_table["Runtime"][i]

            # Get the parallelization scheme used for this simulation
            parallelization_simulation = Parallelization(cores, threads_per_core, processes)

            # Get the parallelization scheme that we have defined for this particular host
            parallelization_for_host = self.launcher.parallelization_for_host(host_id)

            if parallelization_for_host is None: continue

            # Check if the parallelization scheme of the simulations corresponds to the parallelization scheme
            # that is going to be used for the next batch of simulations launched on this host
            if parallelization_simulation == parallelization_for_host and packages_simulation == current_packages:
                # Add the runtime of the simulation to the list of runtimes for the host
                runtimes_for_hosts[host_id].append(runtime)

        return runtimes_for_hosts

    # -----------------------------------------------------------------

    def timeline_paths_for_host(self, runtimes_table, host_id):

        """
        This function ...
        :param runtimes_table:
        :param host_id:
        :return:
        """

        # Initialize a list to contain the paths to the timeline files
        paths = []

        # Loop over the entries in the runtimes table
        for i in range(len(runtimes_table)):

            # Get the simulation name
            simulation_name = runtimes_table["Simulation name"][i]

            # Get the ID of the host and the cluster name for this particular simulation
            host_id_simulation = runtimes_table["Host id"][i]
            cluster_name_simulation = runtimes_table["Cluster name"][i]

            # If the simulation host ID matches the specified host ID, determine the path to the extracted timeline info
            # for that simulation
            if host_id_simulation == host_id:

                timeline_table_path = filesystem.join(self.fit_res_path, simulation_name, "timeline.dat")

                # If the timeline file exists, add its path to the list
                if filesystem.is_file(timeline_table_path): paths.append(timeline_table_path)
                else: log.warning("The timeline table file does not exist for simulation '" + simulation_name + "'")

        # Return the list of timeline table paths
        return paths

    # -----------------------------------------------------------------

    def timeline_paths_for_other_hosts(self, runtimes_table, host_id):

        """
        This function ...
        :param runtimes_table:
        :param host_id:
        :return:
        """

        # Initialize a list to contain the paths to the timeline files
        paths = []

        # Loop over the entries in the runtimes table
        for i in range(len(runtimes_table)):

            # Get the simulation name
            simulation_name = runtimes_table["Simulation name"][i]

            # Get the ID of the host and the cluster name for this particular simulation
            host_id_simulation = runtimes_table["Host id"][i]
            cluster_name_simulation = runtimes_table["Cluster name"][i]

            # If the simulation host ID differs from the specified host ID, determine the path to the extracted timeline info
            # for that simulation
            if host_id_simulation != host_id:

                timeline_table_path = filesystem.join(self.fit_res_path, simulation_name, "timeline.dat")

                # If the timeline file exists, add its path to the list
                if filesystem.is_file(timeline_table_path): paths.append(timeline_table_path)
                else: log.warning("The timeline table file does not exist for simulation '" + simulation_name + "'")

        # Return the list of timeline table paths
        return paths

    # -----------------------------------------------------------------

    def timeline_paths(self, runtimes_table):

        """
        This function ...
        :param runtimes_table:
        :return:
        """

        # Initialize a dictionary
        paths = defaultdict(list)

        # Loop over the entries in the runtimes table
        for i in range(len(runtimes_table)):

            # Get the simulation name
            simulation_name = runtimes_table["Simulation name"][i]

            # Get the ID of the host and the cluster name for this particular simulation
            host_id_simulation = runtimes_table["Host id"][i]
            cluster_name_simulation = runtimes_table["Cluster name"][i]

            timeline_table_path = filesystem.join(self.fit_res_path, simulation_name, "timeline.dat")

            # If the timeline file exists, add its path to the list
            if filesystem.is_file(timeline_table_path): paths[host_id_simulation].append(timeline_table_path)
            else: log.warning("The timeline table file does not exist for simulation '" + simulation_name + "'")

        # Return the dictionary of timeline table paths for each host
        return paths

    # -----------------------------------------------------------------

    def timeline_paths_for_current_packages(self, runtimes_table):

        """
        This function ...
        :param runtimes_table:
        :return:
        """

        # Get the number of photon packages (per wavelength) for this batch of simulations
        current_packages = self.ski.packages()

        # Initialize a dictionary
        paths = defaultdict(list)

        # Loop over the entries in the runtimes table
        for i in range(len(runtimes_table)):

            # Get the simulation name
            simulation_name = runtimes_table["Simulation name"][i]

            # Get the ID of the host and the cluster name for this particular simulation
            host_id_simulation = runtimes_table["Host id"][i]
            cluster_name_simulation = runtimes_table["Cluster name"][i]

            # Get the number of photon packages (per wavelength) used for this simulation
            packages_simulation = runtimes_table["Packages"][i]

            # If the number of photon packages does not match, skip
            if packages_simulation != current_packages: continue

            timeline_table_path = filesystem.join(self.fit_res_path, simulation_name, "timeline.dat")

            # If the timeline file exists, add its path to the list
            if filesystem.is_file(timeline_table_path): paths[host_id_simulation].append(timeline_table_path)
            else: log.warning("The timeline table file does not exist for simulation '" + simulation_name + "'")

        # Return the dictionary of timeline table paths for each host
        return paths

# -----------------------------------------------------------------
