#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.runtimeestimator Contains the RuntimeEstimator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..simulation.parallelization import Parallelization
from ..basics.distribution import Distribution
from ..tools import tables
from ..launch.timing import TimingTable
from ..plot.distribution import DistributionPlotter
from ..basics.map import Map
from ..basics.configurable import Configurable
from ..simulation.skifile import SkiFile
from ..tools import introspection
from ..tools import filesystem as fs
from ..advanced.dustgridtool import DustGridTool

# -----------------------------------------------------------------

class RuntimeEstimator(object):

    """
    This class...
    """

    def __init__(self, timing_table):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        # Set the timing table
        self.timing_table = timing_table

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the timing table
        timing_table = TimingTable.from_file(path)

        # Create the RuntimeEstimator object
        return cls(timing_table)

    # -----------------------------------------------------------------

    def runtime_for(self, ski_file, parallelization, host_id, cluster_name=None, in_path=None, nwavelengths=None, ncells=None, fos=1.2, plot_path=None):

        """
        This function ...
        :param ski_file:
        :param parallelization:
        :param host_id:
        :param cluster_name:
        :param in_path:
        :param nwavelengths:
        :param ncells:
        :param fos: factor of safety
        :param plot_path:
        :return:
        """

        # Get the parameters that are relevant for timing
        parameters = timing_parameters(ski_file, parallelization, host_id, cluster_name, in_path, nwavelengths, ncells)

        # TODO: greatly expand the number of parameters that are used to estimate the runtime
        # such as: nwavelengths, self-absorption, transient heating, data parallel, ...

        # Get the list of runtimes for the specified host for the specified configuration of packages and parallelization
        previous_runtimes = self.previous_runtimes_for(parameters, parallelization)

        # Create the distribution plotter
        plotter = DistributionPlotter()

        # Check how many runtimes were found from previous simulations on the same host and with the same configuration
        if len(previous_runtimes) != 0:

            # Debugging
            log.debug("A total runtime was found for a total of " + str(len(previous_runtimes)) + " simulations that were run on the same remote host and with the same number of photon packages and parallelization")

            # Create a probability distribution from the recorded runtimes
            distribution = Distribution.from_values("Runtime", previous_runtimes, nbins=25)

            # If requested, plot the distribution
            if plot_path is not None:
                title = "Distribution of previously recorded runtimes"
                plotter.add_distribution(distribution, "Test")
                plotter.set_title(title)
                plotter.run(output_path=plot_path)

            # Return the most frequent (most probable) runtime, times the safety factor
            return distribution.most_frequent * fos

        # No previous runtimes were found for the specified host and configuration
        else:

            # Debugging
            log.debug("No simulations were found that had been run on the same remote host and with the same number of photon packages and parallelization")

            # Calculate estimated runtimes based on simulations run on the specified remote host and with the specified number of photon pacakages, but regardless of the parallelization scheme
            estimated_runtimes = self.estimated_runtimes_for(parameters, parallelization)

            # Check how many runtimes could be estimated based on previous simulations on the same host and with the same number of photon packages
            if len(estimated_runtimes) != 0:

                # Debugging
                log.debug("On the basis of " + str(len(estimated_runtimes)) + " previously run simulations with the same number of photon packages, " + str(len(estimated_runtimes)) + " runtimes could be estimated")

                # Create a probability distribution from the estimated runtimes
                distribution = Distribution.from_values("Runtime", estimated_runtimes, nbins=25)

                # If requested, plot the distribution
                if plot_path is not None:
                    title = "Distribution of runtimes estimated based on simulations run on the specified remote host with " \
                             "the same number of photon packages, with various parallelization schemes"
                    plotter.add_distribution(distribution, "Test")
                    plotter.set_title(title)
                    plotter.run(output_path=plot_path)

                # Return the most probable runtime, times the safety factor
                return distribution.most_frequent * fos

            else:

                # Calculate estimated runtimes based on simulations run on no-matter-which remote host with no-matter-which parallelization scheme, but with the specified number of photon packages
                estimated_runtimes = self.estimated_runtimes_for_all_hosts(parameters, parallelization)

                # Check how many runtimes could be estimated based on previous simulations on other hosts
                if len(estimated_runtimes) != 0: #raise RuntimeError("The runtime could not be estimated: no reference with the same number of photon packages")

                    # Create the probability distribution of estimated runtimes
                    distribution = Distribution.from_values("Runtime", estimated_runtimes, nbins=25)

                    # If requested, plot the distribution
                    if plot_path is not None:
                        title = "Distribution of runtimes estimated based on " \
                                  "simulations with the same number of photon packages " \
                                  "on various hosts and with various parallelization " \
                                  "schemes"
                        plotter.add_distribution(distribution, "Test")
                        plotter.set_title(title)
                        plotter.run(output_path=plot_path)

                    # Return the most probable runtime, times the safety factor
                    return distribution.most_frequent * fos

                # If not a single simulation had the same number of photon packages, the runtime could not be estimated (currently?)
                else:

                    estimated_runtimes = self.estimated_runtimes_for_all_hosts_all_npackages(parameters, parallelization)

                    # Create the probability distribution of estimated runtimes
                    distribution = Distribution.from_values("Runtime", estimated_runtimes, nbins=25)

                    # If requested, plot the distribution
                    if plot_path is not None:
                        title = "Distribution of runtimes estimated based on " \
                                "simulations with a different number of photon packages " \
                                "on various hosts and with various parallelization " \
                                "schemes"
                        plotter.add_distribution(distribution, "Test")
                        plotter.set_title(title)
                        plotter.run(output_path=plot_path)

                    # Return the most probable runtime, times the safety factor
                    return distribution.most_frequent * fos

    # -----------------------------------------------------------------

    def previous_runtimes_for(self, parameters, parallelization):

        """
        This function ...
        :param parameters:
        :param parallelization:
        :return:
        """

        # Get host ID and packages (TODO: this set is way to limited)
        host_id = parameters.host_id
        packages = parameters.npackages

        # Initialize a list to contain the runtimes found for the specified host and configuration
        runtimes = []

        # Loop over the entries in the timing table
        # Columns:
        # "Simulation name" / "Submission time"
        # "Host id"
        # "Cluster name"
        # "Cores"
        # "Hyperthreads per core"
        # "Processes"
        # "Packages"
        # "Total runtime"
        # "Serial runtime"
        # "Parallel runtime"
        # "Runtime overhead"
        for i in range(len(self.timing_table)):

            # Get the ID of the host and the cluster name for this particular simulation
            host_id_sim = self.timing_table["Host id"][i]
            cluster_name_sim = self.timing_table["Cluster name"][i]

            # If the host ID of the simulation does not correspond to the specified host ID, skip this entry
            if host_id_sim != host_id: continue

            # Get the parallelization properties for this particular simulation
            parallelization_sim = self.parallelization_for_entry(i)

            # Get the number of photon packages (per wavelength) used for this simulation
            packages_sim = self.timing_table["Packages"][i]

            # Check if the parallelization scheme of the simulation corresponds to the parallelization scheme
            # that is going to be used for the next batch of simulations launched on this host
            if parallelization_sim == parallelization and packages_sim == packages:

                # Get the total runtime for the current simulation
                runtime_sim = self.timing_table["Total runtime"][i]

                # Add the runtime to the list
                runtimes.append(runtime_sim)

        # Return the list of recorded runtimes
        return runtimes

    # -----------------------------------------------------------------

    def estimated_runtimes_for(self, parameters, parallelization):

        """
        This function ...
        :param parameters:
        :param parallelization:
        :return:
        """

        # Get host ID and npackages (TODO: this set is way to limited)
        host_id = parameters.host_id
        packages = parameters.npackages

        # Indices of the simulations in the timing table from the specified host and for the specified number of photon packages
        indices_configuration = tables.find_indices(self.timing_table, [host_id, packages], ["Host id", "Packages"])

        # Debugging
        log.debug(str(len(indices_configuration)) + " simulations were found that were run on the specified host and had the same number of photon packages as the specified amount (regardless of parallelization scheme)")

        # Return the estimated runtimes for the matching entries
        return self.estimated_runtimes_from_entries(indices_configuration, parallelization)

    # -----------------------------------------------------------------

    def estimated_runtimes_for_all_hosts(self, parameters, parallelization):

        """
        This function ...
        :param parameters:
        :param parallelization:
        :return:
        """

        # Get the number of packages
        packages = parameters.npackages

        # Indices of the simulations in the runtime table with the current number of photon packages
        indices_configuration = tables.find_indices(self.timing_table, packages, "Packages")

        # Debugging
        log.debug(str(len(indices_configuration)) + " simulations were found that had the same number of photon packages as the current configuration (regardless of the parallelization scheme or remote host)")

        # Return the estimated runtimes for the matching entries
        return self.estimated_runtimes_from_entries(indices_configuration, parallelization)

    # -----------------------------------------------------------------

    def estimated_runtimes_for_all_hosts_all_npackages(self, parameters, parallelization):

        """
        This function ...
        :param parameters:
        :param parallelization:
        :return:
        """

        # Get the number of packages
        packages = parameters.npackages

        # Get the number of dust cells
        ncells = parameters.ncells

        # Initialize a list to contain the runtimes estimated from previous simulations
        estimated_runtimes = []

        # Loop over the entire timing table
        for i in range(len(self.timing_table)):

            packages_i = self.timing_table["Packages"][i]
            packages_ratio_i = float(packages) / float(packages_i)

            ncells_i = self.timing_table["Dust cells"][i]
            ncells_ratio_i = float(ncells) / float(ncells_i)

            # print("PACKAGES RATIO", packages_ratio_i)
            # print("NCELLS RATIO", ncells_ratio_i)

            # Get the parallelization scheme for this simulation
            parallelization_sim = self.parallelization_for_entry(i)

            setup_time = self.timing_table["Setup time"][i] * ncells_ratio_i # TODO: incorporate the ncells better (number of tree levels, other parameters?)
            writing_time = self.timing_table["Writing time"][i]
            intermediate_time = self.timing_table["Intermediate time"][i]

            stellar_emission_time = self.timing_table["Stellar emission time"][i] * packages_ratio_i
            spectra_calculation_time = self.timing_table["Spectra calculation time"][i]
            dust_emission_time = self.timing_table["Dust emission time"][i] * packages_ratio_i

            communication_time = self.timing_table["Communication time"][i]
            waiting_time = self.timing_table["Waiting time"][i]

            ## THIS IS IDENTICAL AS IN THE FUNCTION BELOW

            # Get the serial runtime, parallel runtime and runtime overhead
            serial = setup_time + writing_time + intermediate_time
            parallel = stellar_emission_time + spectra_calculation_time + dust_emission_time
            overhead = communication_time + waiting_time

            # TODO: the steps below can be more advanced (cores is not necessarily the total number of threads
            # (hyperthreading), hyperthreading gives 30% performance boost?)
            parallel_times_cores = parallel * parallelization_sim.cores
            overhead_per_core = overhead / parallelization_sim.cores

            # Estimate the runtime
            runtime = serial + parallel_times_cores / parallelization.cores + overhead_per_core * parallelization.cores

            # Add the estimated runtime to the list
            estimated_runtimes.append(runtime)

        # Return the list of estimated runtimes
        return estimated_runtimes

    # -----------------------------------------------------------------

    def estimated_runtimes_from_entries(self, indices, parallelization):

        """
        This function ...
        :param self:
        :param indices:
        :param parallelization:
        :return:
        """

        # Initialize a list to contain the runtimes estimated from previous simulations
        estimated_runtimes = []

        # Loop over the matching indices
        for index in indices:

            # Get the parallelization scheme for this simulation
            parallelization_sim = self.parallelization_for_entry(index)

            # Get the serial runtime, parallel runtime and runtime overhead
            serial = self.timing_table["Setup time"][index] + self.timing_table["Writing time"][index] + self.timing_table["Intermediate time"][index]
            parallel = self.timing_table["Stellar emission time"][index] + self.timing_table["Spectra calculation time"][index] + self.timing_table["Dust emission time"][index]
            overhead = self.timing_table["Communication time"][index] + self.timing_table["Waiting time"][index]

            # TODO: the steps below can be more advanced (cores is not necessarily the total number of threads
            # (hyperthreading), hyperthreading gives 30% performance boost?)
            parallel_times_cores = parallel * parallelization_sim.cores
            overhead_per_core = overhead / parallelization_sim.cores

            # Estimate the runtime
            runtime = serial + parallel_times_cores / parallelization.cores + overhead_per_core * parallelization.cores

            # Add the estimated runtime to the list
            estimated_runtimes.append(runtime)

        # Return the list of estimated runtimes
        return estimated_runtimes

    # -----------------------------------------------------------------

    def parallelization_for_entry(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get the number of cores, threads per core and the number of processes
        cores = self.timing_table["Cores"][index]
        threads_per_core = self.timing_table["Threads per core"][index]
        processes = self.timing_table["Processes"][index]

        # Create and return a Parallelization instance
        return Parallelization(cores, threads_per_core, processes)

# -----------------------------------------------------------------

def timing_parameters(ski_file, parallelization, host_id, cluster_name=None, in_path=None, nwavelengths=None, ncells=None):

    """
    This function ...
    :param ski_file:
    :param parallelization:
    :param host_id:
    :param cluster_name:
    :param in_path:
    :param nwavelengths:
    :param ncells:
    :return:
    """

    #values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths, packages,
    #          cells, selfabsorption, transient_heating, data_parallel, total_runtime, setup_time, stellar_time,
    #          spectra_time, dust_time, writing_time, waiting_time, communication_time, intermediate_time]

    parameters = Map()

    # Remote host
    parameters.host_id = host_id
    parameters.cluster_name = cluster_name

    # Parallelization
    parameters.cores = parallelization.cores
    parameters.threads_per_core = parallelization.threads_per_core
    parameters.processes = parallelization.processes

    # Ski file parameters
    try: parameters.nwavelengths = ski_file.nwavelengths()
    except ValueError:
        if nwavelengths is not None: parameters.nwavelengths = nwavelengths
        elif in_path is not None: parameters.nwavelengths = ski_file.nwavelengthsfile(in_path)
        else: raise ValueError("Cannot determine the number of wavelengths: either the input path should be specified (so that the wavelengths file can be found) or the number of wavelengths should be specified explicitly")
    parameters.npackages = ski_file.packages()
    try: parameters.ncells = ski_file.ncells()
    except ValueError:

        if ncells is None: raise ValueError("Cannot determine the number of dust cells: if a tree dust grid structure is used, the approximate number of cells should be specified explicitly")

        # Set the number of dust cells
        parameters.ncells = ncells

        #parameters.ncells = None # the number of dust cells cannot be predicted if using a tree
        # dust grid, but we esimate runtimes based on ski files with the same model anyway, so the number of dust cells
        # can be expected to be the same (if the dust grid parameters are identical, that is...)

        #parameters.grid_type = None
        #parameters.rel_scale = None
        #parameters.min_level = None
        #parameters.max_mass_fraction = None

    parameters.selfabsorption = ski_file.dustselfabsorption()
    parameters.transient_heating = ski_file.transientheating()

    # Other
    #parameters.data_parallel = data_parallel
    parameters.data_parallel = parallelization.data_parallel

    # Return the parameters
    return parameters

# -----------------------------------------------------------------
