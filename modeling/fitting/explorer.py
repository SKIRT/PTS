#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.explorer Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import abstractmethod

# Import the relevant PTS classes and modules
from .launcher import FittingModelLauncher
from ...core.tools.logging import log
from ...core.launch.parallelization import Parallelization
from ...magic.animation.scatter import ScatterAnimation

# -----------------------------------------------------------------

class ParameterExplorer(FittingModelLauncher):
    
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

        # -- Attributes --

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExplorer, self).setup()

        # Get the names of the filters for which we have photometry
        self.filter_names = self.get_observed_filter_names()

        # Set options for the BatchLauncher: basic options
        self.launcher.config.shared_input = True  # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = self.config.remotes  # the remote hosts on which to run the simulations

        # Set options for the BatchLauncher: simulation analysis options
        self.launcher.config.analysis.extraction.path = self.fit_res_path
        self.launcher.config.analysis.misc.path = self.fit_res_path # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.plotting.path = self.fit_plot_path # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.extraction.timeline = True # extract the simulation timeline
        self.launcher.config.analysis.plotting.seds = True  # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_path # the path to the reference SED (for plotting the simulated SED against the reference points)
        #self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_dustpedia_path # the path to the DustPedia SED
        #self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        #self.launcher.config.analysis.misc.images = True  # Make observed images
        self.launcher.config.analysis.misc.observation_filters = self.filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.plotting.format = "png" # plot in PNG format so that an animation can be made from the fit SEDs
        self.launcher.config.analysis.timing_table_path = self.timing_table_path # The path to the timing table file
        self.launcher.config.analysis.memory_table_path = self.memory_table_path # The path to the memory table file

        # Set remote for the 'extra' simulations
        self.launcher.config.extra_remote = "nancy"

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Load the ski file
        self.load_ski()

    # -----------------------------------------------------------------

    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the animations ...")

        # Initialize the scatter animation
        self.scatter_animation = ScatterAnimation([self.config.young_stars.min, self.config.young_stars.max],
                                                  [self.config.ionizing_stars.min, self.config.ionizing_stars.max],
                                                  [self.config.dust.min, self.config.dust.max])
        self.scatter_animation.x_label = "FUV luminosity of young stars"
        self.scatter_animation.y_label = "FUV luminosity of ionizing stars"
        self.scatter_animation.z_label = "Dust mass"

        # No distribution animations for the Explorer classes
        self.fuv_young_animation = None
        self.fuv_ionizing_animation = None
        self.dust_mass_animation = None

    # -----------------------------------------------------------------

    @abstractmethod
    def set_parameters(self):

        """
        This function ...
        :return:
        """

        return

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

        # not necessary

        return

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
