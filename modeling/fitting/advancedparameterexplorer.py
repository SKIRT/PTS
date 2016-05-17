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

# Import standard modules
import io
import imageio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde

# Import the relevant PTS classes and modules
from .parameterexploration import ParameterExplorer
from ...core.tools import tables, time
from ...core.tools.logging import log
from ...core.launch.options import SchedulingOptions
from ...core.launch.parallelization import Parallelization
from ...core.launch.runtime import RuntimeEstimator
from ...core.basics.distribution import Distribution
from ...core.basics.animatedgif import AnimatedGif
from ...core.tools import filesystem as fs
from ...core.plot.distribution import DistributionPlotter

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

        # -- Attributes --

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

        # Make visualisations
        explorer.config.visualise = arguments.visualise

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
        #self.set_parameters()

        # 5. Set the parallelization schemes for the different remote hosts
        self.set_parallelization()

        # 6. Estimate the runtimes for the different remote hosts
        self.estimate_runtimes()

        exit()

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

            # Normalize the distribution
            distribution.normalize(value=1.0, method="max")

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

        # Create animations
        if self.config.visualise:
            points_animation = AnimatedGif()
            fuv_young_animation = AnimatedGif()
            fuv_ionizing_animation = AnimatedGif()
            dust_mass_animation = AnimatedGif()
        else:
            points_animation = None
            fuv_young_animation = None
            fuv_ionizing_animation = None
            dust_mass_animation = None

        # Debugging
        if log.is_debug() and False:

            # Young stars
            x_limits = [self.config.young_stars.min, self.config.young_stars.max]
            #print(self.distributions["FUV young"].cumulative_smooth(x_limits[0], x_limits[1]))
            self.distributions["FUV young"].plot_smooth(x_limits=x_limits, title="Probability distribution from which FUV luminosities of young stars will be drawn")
            self.distributions["FUV young"].plot_smooth(x_limits=x_limits, title="Probability distribution from which FUV luminosities of young stars will be drawn (in log scale)")
            self.distributions["FUV young"].plot_cumulative_smooth(x_limits=x_limits, title="Cumulative distribution of FUV luminosities of young stars")

            # Ionizing stars
            x_limits = [self.config.ionizing_stars.min, self.config.ionizing_stars.max]
            #print(self.distributions["FUV ionizing"].cumulative_smooth(x_limits[0], x_limits[1]))
            self.distributions["FUV ionizing"].plot_smooth(x_limits=x_limits, title="Probability distribution from which FUV luminosities of ionizing stars will be drawn")
            self.distributions["FUV ionizing"].plot_smooth(x_limits=x_limits, title="Probability distribution from which FUV luminosities of ionizing stars will be drawn (in log scale)")
            self.distributions["FUV ionizing"].plot_cumulative_smooth(x_limits=x_limits, title="Cumulative distribution of FUV luminosities of ionizing stars")

            # Dust mass
            x_limits = [self.config.dust.min, self.config.dust.max]
            #print(self.distributions["Dust mass"].cumulative_smooth(x_limits[0], x_limits[1]))
            self.distributions["Dust mass"].plot_smooth(x_limits=x_limits, title="Probability distribution from which dust masses will be drawn")
            self.distributions["Dust mass"].plot_smooth(x_limits=x_limits, title="Probability distribution from which dust masses will be drawn (in log scale)")
            self.distributions["Dust mass"].plot_cumulative_smooth(x_limits=x_limits, title="Cumulative distribution of dust masses")

        # Draw parameters values for the specified number of simulations
        #for counter in range(self.config.simulations):
        for counter in range(100):

            # Debugging
            log.debug("Calculating random parameter set " + str(counter+1) + " of " + str(self.config.simulations) + " ...")

            # Draw a random FUV luminosity of the young stellar population
            young_luminosity = self.distributions["FUV young"].random(self.config.young_stars.min, self.config.young_stars.max)

            # Draw a random FUV luminosity of the ionizing stellar population
            ionizing_luminosity = self.distributions["FUV ionizing"].random(self.config.ionizing_stars.min, self.config.ionizing_stars.max)

            # Draw a random dust mass
            dust_mass = self.distributions["Dust mass"].random(self.config.dust.min, self.config.dust.max)

            # Add the parameter values to the dictionary
            self.parameters["FUV young"].append(young_luminosity)
            self.parameters["FUV ionizing"].append(ionizing_luminosity)
            self.parameters["Dust mass"].append(dust_mass)

            # Add a frame to the animation of parameter points
            if points_animation is not None and self.number_of_models > 1:

                buf = io.BytesIO()

                fig = plt.figure(figsize=(15,15))

                #ax = Axes3D(fig)

                # Add first subplot
                ax = fig.add_subplot(2, 2, 1, projection='3d')

                ax.scatter(self.parameters["FUV young"], self.parameters["FUV ionizing"], self.parameters["Dust mass"])
                ax.set_xlim([self.config.young_stars.min, self.config.young_stars.max])
                ax.set_ylim([self.config.ionizing_stars.min, self.config.ionizing_stars.max])
                ax.set_zlim([self.config.dust.min, self.config.dust.max])
                ax.set_xlabel("FUV luminosity of young stars")
                ax.set_ylabel("FUV luminosity of ionizing stars")
                ax.set_zlabel("Dust mass")
                # To draw projected points against the axis planes:
                #ax.plot(self.parameters["FUV young"], self.parameters["Dust mass"], 'r+', zdir='y', zs=self.config.ionizing_stars.max)
                #ax.plot(self.parameters["FUV ionizing"], self.parameters["Dust mass"], 'g+', zdir='x', zs=self.config.young_stars.min)
                #ax.plot(self.parameters["FUV young"], self.parameters["FUV ionizing"], 'k+', zdir='z', zs=self.config.dust.min)

                # Add second subplot
                ax = fig.add_subplot(2, 2, 2)

                # Density plot of FUV young vs. FUV ionizing
                x = np.array(self.parameters["FUV young"])
                y = np.array(self.parameters["FUV ionizing"])
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)
                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='')
                ax.set_xlabel("FUV luminosity of young stars")
                ax.set_ylabel("FUV luminosity of ionizing stars")
                ax.set_xlim([self.config.young_stars.min, self.config.young_stars.max])
                ax.set_ylim([self.config.ionizing_stars.min, self.config.ionizing_stars.max])

                # Add third subplot
                ax = fig.add_subplot(2, 2, 3)

                # Density plot of FUV young vs. dust mass
                x = np.array(self.parameters["FUV young"])
                y = np.array(self.parameters["Dust mass"])
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)
                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='')
                ax.set_xlabel("FUV luminosity of young stars")
                ax.set_ylabel("Dust mass")
                ax.set_xlim([self.config.young_stars.min, self.config.young_stars.max])
                ax.set_ylim([self.config.dust.min, self.config.dust.max])

                # Add fourth subplot
                ax = fig.add_subplot(2, 2, 4)

                # Density plot of FUV ionizing vs. dust mass
                x = np.array(self.parameters["FUV ionizing"])
                y = np.array(self.parameters["Dust mass"])
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)
                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='')
                ax.set_xlabel("FUV luminosity of ionizing stars")
                ax.set_ylabel("Dust mass")
                ax.set_xlim([self.config.ionizing_stars.min, self.config.ionizing_stars.max])
                ax.set_ylim([self.config.dust.min, self.config.dust.max])

                plt.tight_layout()

                plt.savefig(buf, format="png")
                plt.close()
                buf.seek(0)
                im = imageio.imread(buf)
                buf.close()
                points_animation.add_frame(im)

            plotter = DistributionPlotter()

            # Add a frame to the animation of the distribution of the FUV luminosity of young starss
            if fuv_young_animation is not None and self.number_of_models > 1:

                new_distribution = Distribution.from_values(self.parameters["FUV young"])
                new_distribution.normalize(1.0, method="max")
                buf = io.BytesIO()
                plotter.add_distribution(self.distributions["FUV young"], "Previous models")
                plotter.add_distribution(new_distribution, "New models")
                plotter.set_variable_name("FUV luminosity of young stars")
                plotter.run(buf, format="png", min_value=self.config.young_stars.min, max_value=self.config.young_stars.max, max_count=1., logscale=True)
                #new_distribution.plot(title="FUV luminosity of young stars",
                #                      x_limits=[self.config.dust.min, self.config.dust.max], path=buf, format="png")
                buf.seek(0)
                im = imageio.imread(buf)
                buf.close()
                fuv_young_animation.add_frame(im)

            # Add a frame to the animation of the distribution of the FUV luminosity of ionizing stars
            if fuv_ionizing_animation is not None and self.number_of_models > 1:

                new_distribution = Distribution.from_values(self.parameters["FUV ionizing"])
                new_distribution.normalize(1.0, method="max")
                buf = io.BytesIO()
                plotter.clear()
                plotter.add_distribution(self.distributions["FUV ionizing"], "Previous models")
                plotter.add_distribution(new_distribution, "New models")
                plotter.set_variable_name("FUV luminosity of ionizing stars")
                plotter.run(buf, format="png", min_value=self.config.ionizing_stars.min, max_value=self.config.ionizing_stars.max, max_count=1., logscale=True)
                #new_distribution.plot(title="FUV luminosity of ionizing stars",
                #                      x_limits=[self.config.dust.min, self.config.dust.max], path=buf, format="png")
                buf.seek(0)
                im = imageio.imread(buf)
                buf.close()
                fuv_ionizing_animation.add_frame(im)

            # Add a frame to the animation of the distribution of the dust mass
            if dust_mass_animation is not None and self.number_of_models > 1:

                new_distribution = Distribution.from_values(self.parameters["Dust mass"])
                new_distribution.normalize(1.0, method="max")
                buf = io.BytesIO()
                plotter.clear()
                plotter.add_distribution(self.distributions["Dust mass"], "Previous models")
                plotter.add_distribution(new_distribution, "New models")
                plotter.set_variable_name("Dust mass")
                plotter.run(buf, format="png", min_value=self.config.dust.min, max_value=self.config.dust.max, max_count=1., logscale=True)
                #new_distribution.plot(title="Dust mass", x_limits=[self.config.dust.min, self.config.dust.max],
                #                      path=buf, format="png", other=self.distributions["Dust mass"])
                buf.seek(0)
                im = imageio.imread(buf)
                buf.close()
                dust_mass_animation.add_frame(im)

        # Save the animation of the parameter values as points in a 3D plot
        if points_animation is not None:
            path = fs.join(self.visualisation_path, time.unique_name("advancedparameterexploration") + ".gif")
            points_animation.save(path)

        # Save the animation of the distribution of values for the FUV luminosity of the young stars
        if fuv_young_animation is not None:
            path = fs.join(self.visualisation_path, time.unique_name("advancedparameterexploration_fuvyoung") + ".gif")
            fuv_young_animation.save(path)

        # Save the animation of the distribution of values for the FUV luminosity of the ionizing stars
        if fuv_ionizing_animation is not None:
            path = fs.join(self.visualisation_path, time.unique_name("advancedparameterexploration_fuvionizing") + ".gif")
            fuv_ionizing_animation.save(path)

        # Save the animation of the distribution of values for the dust mass
        if dust_mass_animation is not None:
            path = fs.join(self.visualisation_path, time.unique_name("advancedparameterexploration_dustmass") + ".gif")
            dust_mass_animation.save(path)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the correponding system).
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
