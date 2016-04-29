#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.advancedparameterexplorer Contains the AdvancedParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .parameterexploration import ParameterExplorer
from ...core.tools import tables
from ...core.tools.logging import log
from ...core.launch.options import SchedulingOptions
from ...core.launch.parallelization import Parallelization
from ...core.launch.runtime import RuntimeEstimator
from ...core.basics.distribution import Distribution

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

        # The probability distributions for the different fit parameters
        self.distributions = dict()

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

        # Set the number of simulations to launch in the batch
        if arguments.simulations is not None: explorer.config.simulations = arguments.simulations

        # Set the remote host IDs
        if arguments.remotes is not None: explorer.config.remotes = arguments.remotes

        # Set the limits of the FUV luminosity of the young stellar population
        if arguments.young is not None:
            explorer.config.young_stars.min = arguments.young[0]
            explorer.config.young_stars_max = arguments.young[1]

        # Set the limits of the FUV luminosity of the ionizing stellar population
        if arguments.ionizing is not None:
            explorer.config.ionizing_stars.min = arguments.ionizing[0]
            explorer.config.ionizing_stars.max = arguments.ionizing[1]

        # Set the limits of the dust mass
        if arguments.dust is not None:
            explorer.config.dust.min = arguments.dust[0]
            explorer.config.dust.max = arguments.dust[1]

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

        # Load the probability distributions for the different parameters
        self.load_distributions()

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

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the probability distributions for the different fit parameters ...")

        # Loop over the different fit parameters
        for parameter_name in self.parameter_names:

            # Load the probability distribution
            distribution = Distribution.from_file(self.distribution_table_paths[parameter_name])

            # Set the distribution
            self.distributions[parameter_name] = distribution

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Picking random parameter values based on the probability distributions ...")

        # Draw parameters values for the specified number of simulations
        for _ in range(self.config.simulations):

            # Draw a random FUV luminosity of the young stellar population
            young_luminosity = self.distributions["FUV young"].random()

            # Draw a random FUV luminosity of the ionizing stellar population
            ionizing_luminosity = self.distributions["FUV ionizing"].random()

            # Draw a random dust mass
            dust_mass = self.distributions["Dust mass"].random()

            # Add the combination of parameter values to the list
            combination = (young_luminosity, ionizing_luminosity, dust_mass)
            self.parameters.append(combination)

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

        # Inform the user
        log.info("Estimating the runtimes based on the results of previous runs ...")

        # Get the number of photon packages (per wavelength) for this batch of simulations
        current_packages = self.ski.packages()

        # Debugging
        log.debug("Loading the table with the total runtimes of previous simulations ...")

        # Load the timing table
        timing_table = tables.from_file(self.timing_table_path, format="ascii.ecsv")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator(timing_table)

        # Initialize a dictionary to contain the estimated walltimes for the different hosts with scheduling system
        walltimes = dict()

        # Loop over the hosts which use a scheduling system and estimate the walltime
        for host_id in self.launcher.scheduler_host_ids:

            # Get the parallelization scheme that we have defined for this remote host
            parallelization = self.launcher.parallelization_for_host(host_id)

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = estimator.runtime_for(host_id, current_packages, parallelization)

            # Set the estimated walltime
            walltimes[host_id] = runtime

        # Create and set scheduling options for each host that uses a scheduling system
        for host_id in walltimes: self.scheduling_options[host_id] = SchedulingOptions.from_dict({"walltime": walltimes[host_id]})

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
