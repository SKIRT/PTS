#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.evolver Contains the ModelEvolver class, an abstract base class for GeneticModelEvolver
#  and BasicModelEvolver

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from .launcher import FittingModelLauncher
from ...core.tools import time
from ...core.tools.logging import log
from ...core.launch.options import SchedulingOptions
from ...core.launch.parallelization import Parallelization
from ...core.launch.runtime import RuntimeEstimator
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ModelEvolver(FittingModelLauncher):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelEvolver, self).__init__(config)

        # -- Attributes --

        # The probability distributions for the different fit parameters
        self.distributions = dict()

    # -----------------------------------------------------------------

    @abstractmethod
    def load_input(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def set_parameters(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    def update_animations(self, young_luminosity, ionizing_luminosity, dust_mass):

        """
        This function ...
        :param young_luminosity:
        :param ionizing_luminosity:
        :param dust_mass:
        :return:
        """

        # Add the point (and thus a frame) to the animation of parameter points
        self.scatter_animation.add_point(young_luminosity, ionizing_luminosity, dust_mass)

        # Update the distribution animations
        if self.number_of_models > 1:

            # Add a frame to the animation of the distribution of the FUV luminosity of young starss
            self.fuv_young_animation.add_value(young_luminosity)

            # Add a frame to the animation of the distribution of the FUV luminosity of ionizing stars
            self.fuv_ionizing_animation.add_value(ionizing_luminosity)

            # Add a frame to the animation of the distribution of the dust mass
            self.dust_mass_animation.add_value(dust_mass)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the corresponding system).
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for the remote host(s) that use a scheduling system ...")

        # Loop over the IDs of the hosts used by the batch launcher that use a scheduling system
        for host in self.launcher.scheduler_hosts:

            # Debugging
            log.debug("Setting the parallelization scheme for host '" + host.id + "' ...")

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

            # Debugging
            log.debug("Parallelization scheme: " + str(parallelization))

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

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator.from_file(self.timing_table_path)

        # Initialize a dictionary to contain the estimated walltimes for the different hosts with scheduling system
        walltimes = dict()

        # Loop over the hosts which use a scheduling system and estimate the walltime
        for host_id in self.launcher.scheduler_host_ids:

            # Debugging
            log.debug("Estimating the runtime for host '" + host_id + "' ...")

            # Get the parallelization scheme that we have defined for this remote host
            parallelization = self.launcher.parallelization_for_host(host_id)

            # Visualisation of the distribution of estimated runtimes
            if self.config.visualise: plot_path = fs.join(self.visualisation_path, time.unique_name("advancedparameterexploration_runtime_"+host_id) + ".pdf")
            else: plot_path = None

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = estimator.runtime_for(host_id, current_packages, parallelization, plot_path=plot_path)

            # Debugging
            log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

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

        # Write the animations
        if self.config.visualise: self.write_animations()

# -----------------------------------------------------------------
