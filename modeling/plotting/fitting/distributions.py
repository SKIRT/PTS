#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.distributions

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.simulation.wavelengthgrid import WavelengthGrid
from ....core.basics.distribution import Distribution
from ....core.plot.distribution import DistributionPlotter

# -----------------------------------------------------------------

class DistributionsPLotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DistributionsPLotter, self).__init__(*args, **kwargs)

        # The prior and posterior probability distributions
        self.prior_distributions = dict()  # indexed on generation, then on the free parameter labels
        self.posterior_distributions = dict()  # indexed on the free parameter labels

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DistributionsPLotter, self).setup(**kwargs)

        # Directories for prior and posterior distributions
        self.plot_fitting_prior_distributions_path = fs.create_directory_in(self.plot_fitting_distributions_path,
                                                                            "prior")
        self.plot_fitting_posterior_distributions_path = fs.create_directory_in(self.plot_fitting_distributions_path,
                                                                                "posterior")

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the distributions ...")

        # Load prior distributions
        self.load_prior_distributions()

        # Load the posterior distributions
        self.load_posterior_distributions()

    # -----------------------------------------------------------------

    def load_prior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prior probability distributions for the different free parameters for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Initialize a dictionary
            distributions_dict = dict()

            # Load the parameters table
            parameters_table = self.parameters_table_for_generation(generation_name)

            # Loop over the different parameters
            for label in self.free_parameter_labels:

                # Test if the values are not all equal
                if np.min(parameters_table[label]) == np.max(parameters_table[label]):

                    # Give a warning and skip this parameter (we can't make a distribution of values that are all the same)
                    log.warning("The prior values for the '" + label + "' parameter of the '" + generation_name + "' generation are all identical")
                    continue

                # Create a distribution from the list of parameter values
                distribution = Distribution.from_values("Parameter value", parameters_table[label])

                # Add the distribution to the dictionary
                distributions_dict[label] = distribution

            # Add the dictionary for this generation
            self.prior_distributions[generation_name] = distributions_dict

    # -----------------------------------------------------------------

    def load_posterior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the probability distributions for the different fit parameters ...")

        # Loop over the different fit parameters
        for parameter_name in self.free_parameter_labels:

            # Check if the posterior distribution has been calculated
            if not self.has_distribution(parameter_name): continue

            # Load and set the distribution
            self.posterior_distributions[parameter_name] = self.get_parameter_distribution(parameter_name)

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Plot prior distributions
        self.plot_prior_distributions()

        # Plot posterior distributions
        self.plot_posterior_distributions()

    # -----------------------------------------------------------------

    def plot_prior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the posterior parameter distributions ...")

        # Create a distribution plotter
        plotter = DistributionPlotter()

        # Loop over the generations
        for generation_name in self.prior_distributions:

            # Create directory for this generation
            generation_path = fs.create_directory_in(self.plot_fitting_prior_distributions_path, generation_name)

            # Get the parameter ranges for this generation
            ranges = self.parameter_ranges_for_generation(generation_name)

            # Directories
            linear_path = fs.create_directory_in(generation_path, "linear")
            linear_smooth_path = fs.create_directory_in(generation_path, "linear smooth")
            log_path = fs.create_directory_in(generation_path, "log")
            log_smooth_path = fs.create_directory_in(generation_path, "log smooth")

            # Loop over the distributions (for the different parameters)
            for label in self.prior_distributions[generation_name]:
                # Get the limits
                limits = [ranges[label].min, ranges[label].max]

                # Show the limits
                # print(self.prior_distributions[generation_name][label].min_value)
                # print(self.prior_distributions[generation_name][label].max_value)

                # Define the title
                title = "Prior probability distribution of " + self.parameter_descriptions[label]

                # Add the distribution
                plotter.add_distribution(self.prior_distributions[generation_name][label], label)
                plotter.title = title

                # Plot linear
                path = fs.join(linear_path, label + ".pdf")
                plotter.run(output_path=path)

                plotter.clear(clear_distributions=False)

                # Plot linear smooth
                path = fs.join(linear_smooth_path, label + ".pdf")
                plotter.run(output_path=path, add_smooth=True)

                plotter.clear(clear_distributions=False)

                # Plot log
                path = fs.join(log_path, label + ".pdf")
                plotter.run(output_path=path, logscale=True)

                plotter.clear(clear_distributions=False)

                # Plot log smooth
                path = fs.join(log_smooth_path, label + ".pdf")
                plotter.run(output_path=path, add_smooth=True, logscale=True)

                # Clear
                plotter.clear()

    # -----------------------------------------------------------------

    def plot_posterior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the posterior parameter distributions ...")

        # Loop over the parameter labels
        for label in self.posterior_distributions:
            # Path
            linear_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_linear.pdf")
            log_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_log.pdf")
            cumulative_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_cumulative.pdf")

            # Get the limits
            x_limits = [self.free_parameter_ranges[label].min, self.free_parameter_ranges[label].max]

            # Plot the distributions
            self.posterior_distributions[label].plot_smooth(x_limits=x_limits,
                                                            title="Posterior probability distribution of " +
                                                                  self.parameter_descriptions[label],
                                                            path=linear_path)
            self.posterior_distributions[label].plot_smooth(x_limits=x_limits,
                                                            title="Posterior probability distribution of " +
                                                                  self.parameter_descriptions[
                                                                      label] + " (in log scale)", path=log_path)
            self.posterior_distributions[label].plot_cumulative_smooth(x_limits=x_limits,
                                                                       title="Cumulative distribution of " +
                                                                             self.parameter_descriptions[label],
                                                                       path=cumulative_path)

# -----------------------------------------------------------------
