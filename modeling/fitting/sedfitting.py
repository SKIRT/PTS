#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.sedfitting Contains the SEDFitter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import imageio
from matplotlib import pyplot as plt

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import time
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.basics.distribution import Distribution
from ...core.basics.animation import Animation
from .tables import ModelProbabilitiesTable, ParameterProbabilitiesTable
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class SEDFitter(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDFitter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The model probabilities tables
        self.model_probabilities = dict()

        # The parameter probabilities tables
        self.parameter_probabilities = dict()

        # The model parameter table
        self.parameter_tables = dict()

        # The tables with the probability distributions for the different fit parameters
        self.distributions = dict()

        # The animation
        self.animation = None

        # The directory with the probability tables for all finished generations
        self.prob_generations_path = None

        # The dictionary with ...
        self.prob_generations_table_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the parameters of the best models for each generation
        self.get_best_parameters()

        # 3. Calculate the probabilities
        self.calculate_probabilities()

        # 4. Calculate the probability distributions
        self.create_distributions()

        # 5. Make an animation of the fitting procedure
        if self.config.visualise: self.animate()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    @property
    def has_finished_generations(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.has_finished_generations

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        if self.config.unfinished: return self.fitting_run.generation_names
        else: return self.fitting_run.finished_generations

    # -----------------------------------------------------------------

    @property
    def ngenerations(self):

        """
        This function ...
        :return:
        """

        return len(self.generation_names)

    # -----------------------------------------------------------------

    @property
    def has_generations(self):

        """
        This function ...
        :return:
        """

        return self.ngenerations > 0

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDFitter, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # The directory with the probability tables for all finished generations
        self.prob_generations_path = fs.create_directory_in(self.fitting_run.prob_path, "generations")

        # Check if there are finished generations
        if not self.has_finished_generations and not self.config.unfinished: raise RuntimeError("There are no finished generations")
        if not self.has_generations: raise RuntimeError("There are no generations")

        # For each finished generation, determine the path to the probability table
        for generation_name in self.generation_names:

            # Determine probabilities table path
            path = fs.join(self.prob_generations_path, generation_name + ".dat")
            self.prob_generations_table_paths[generation_name] = path

    # -----------------------------------------------------------------

    @property
    def only_finished(self):

        """
        This function ...
        :return:
        """

        return not self.config.unfinished

    # -----------------------------------------------------------------

    def get_best_parameters(self):

        """"
        This function ...
        """

        # Inform the user
        log.info("Getting the parameter values of the best model for the finished generations (if not already done) ...")

        # Loop over the finished generations
        for generation_name in self.generation_names:

            # Check if the generation is already in the best parameters table
            if generation_name in self.fitting_run.best_parameters_table.generation_names: continue

            # Otherwise, add the best parameter values
            values, chi_squared = self.fitting_run.best_parameter_values_for_generation(generation_name, return_chi_squared=True, only_finished=self.only_finished)

            # Add an entry to the best parameters table file
            self.fitting_run.best_parameters_table.add_entry(generation_name, values, chi_squared)

    # -----------------------------------------------------------------

    def calculate_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities ...")

        # Calculate the probability tables
        self.calculate_model_probabilities()

        # Calcualte the combined probability tables (all finished generations)
        self.calculate_parameter_probabilities()

    # -----------------------------------------------------------------

    def calculate_model_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the model probabilities for the finished generations (if not already done) ...")

        # Loop over the finished generations
        for generation_name in self.generation_names:

            # Check whether the probabilities table is already present for this generation
            if fs.is_file(self.prob_generations_table_paths[generation_name]):

                # Debugging
                log.debug("Loading the model probabilities table for generation " + generation_name + " ...")

                # Load the probabilities table
                probabilities_table = ModelProbabilitiesTable.from_file(self.prob_generations_table_paths[generation_name])

                # Add to the dictionary
                self.model_probabilities[generation_name] = probabilities_table

            # Otherwise, calculate the probabilities based on the chi squared table
            else:

                # Load the parameter table
                parameter_table = self.fitting_run.parameters_table_for_generation(generation_name)

                # Load the chi squared table
                chi_squared_table = self.fitting_run.chi_squared_table_for_generation(generation_name)

                # Sort the table for decreasing chi squared value
                chi_squared_table.sort("Chi squared")
                chi_squared_table.reverse()

                # Get the chi squared values
                chi_squared_values = chi_squared_table["Chi squared"]

                # Calculate the probability for each model
                probabilities = np.exp(-0.5 * chi_squared_values)

                # Create the probabilities table
                probabilities_table = ModelProbabilitiesTable(parameters=self.fitting_run.free_parameter_labels, units=self.fitting_run.parameter_units)

                # Add the entries to the model probabilities table
                for i in range(len(chi_squared_table)):

                    # Get the simulation name
                    simulation_name = chi_squared_table["Simulation name"][i]

                    # Get a dictionary with the parameter values for this simulation
                    parameter_values = parameter_table.parameter_values_for_simulation(simulation_name)

                    # Add an entry to the table
                    probabilities_table.add_entry(simulation_name, parameter_values, probabilities[i])

                # Save the model probabilities table
                probabilities_table.saveto(self.prob_generations_table_paths[generation_name])

                # Add to the dictionary
                self.model_probabilities[generation_name] = probabilities_table

    # -----------------------------------------------------------------

    def calculate_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities of the different parameter values ...")

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Create a set for the unique values
            unique_values = set()

            # Loop over the generations to extract all unique values for the parameter
            for generation_name in self.model_probabilities:

                # Loop over the values of this parameter for this generation, and expand the set accordingly
                for value in self.model_probabilities[generation_name][label]: unique_values.add(value)

            # Get a (sorted) list of all the unique values for this parameter
            unique_values = sorted(list(unique_values))

            # Initialize a ParameterProbabilitiesTable instance for this parameter
            table = ParameterProbabilitiesTable()

            # Add an entry for each unique parameter value that has been encountered
            for value in unique_values:

                # Add the probabilities from all models that have this value
                individual_probabilities = []

                # Loop over the generations
                for generation_name in self.model_probabilities:
                    simulation_indices = self.model_probabilities[generation_name][label] == value
                    nsimulations_for_value = np.sum(simulation_indices)
                    individual_probabilities += list(self.model_probabilities[generation_name]["Probability"][simulation_indices])

                # Combine the individual probabilities
                combined_probability = np.sum(np.array(individual_probabilities))

                # Add an entry to the table
                table.add_entry(value, combined_probability)

            # Set the table
            self.parameter_probabilities[label] = table

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the probability distributions ...")

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Convert the probability lists into NumPy arrays and normalize them
            normalized_probabilities = np.array(self.parameter_probabilities[label]["Probability"]) / sum(self.parameter_probabilities[label]["Probability"])

            # Create the probability distributions for the different parameters
            self.distributions[label] = Distribution.from_probabilities(normalized_probabilities, self.parameter_probabilities[label]["Value"], label)

    # -----------------------------------------------------------------

    def animate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating an animation of the SED fitting procedure ...")

        # Create an Animation instance
        self.animation = Animation()

        # Loop over the generations
        for generation_name in self.generation_names:

            # Get the model probabilities table for this generation
            model_probabilities_table = self.model_probabilities[generation_name]

            # Loop over the simulations in the model probability table
            for i in range(len(model_probabilities_table)):

                # Get the name of the simulation
                simulation_name = model_probabilities_table["Simulation name"][i]

                # Determine the path to the corresponding SED plot file
                path = fs.join(self.fitting_run.generations_path, generation_name, simulation_name, "plot", "sed.png")

                # Load the image (as a NumPy array)
                image = imageio.imread(path)

                # Add the image to the animation
                self.animation.add_frame(image)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the parameter probabilities
        self.write_parameter_probabilities()

        # Write the ski file of the best simulation
        self.write_best_parameters()

        # Write the probability distributions in table format
        self.write_distributions()

        # Plot the probability distributions as histograms
        self.plot_distributions()

        # Write the animated GIF
        if self.config.visualise: self.write_animation()

    # -----------------------------------------------------------------

    def write_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter probabilities ...")

        # Loop over the probability tables for the different free parameter
        for label in self.fitting_run.free_parameter_labels:

            # Save the table
            self.parameter_probabilities[label].saveto(self.fitting_run.get_parameter_probabilities_path(label))

    # -----------------------------------------------------------------

    def write_best_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the best model parameters table ...")

        # Save the best parameters table
        self.fitting_run.best_parameters_table.save()

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the probability distributions ...")

        # Loop over the entries in the 'probabilities' table
        for label in self.distributions:

            # Debugging
            log.debug("Writing the probability distribution of the " + self.fitting_run.parameter_descriptions[label] + " ...")

            # Write the table of probabilities for this parameter
            self.distributions[label].saveto(self.fitting_run.get_parameter_distribution_path(label))

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Loop over the different fit parameters
        for parameter_name in self.distributions:

            # Debugging
            log.debug("Plotting the probability distribution of the " + self.fitting_run.parameter_descriptions[parameter_name] + " ...")

            # Get the probability distribution for this parameter
            distribution = self.distributions[parameter_name]
            description = self.fitting_run.parameter_descriptions[parameter_name]

            # Create a plot file for the probability distribution
            path = fs.join(self.fitting_run.prob_distributions_path, parameter_name + ".pdf")
            try: distribution.plot(title="Probability of the " + description, path=path, logscale=False, xlogscale=True)
            except ValueError: log.warning("Could not create the distribution plot for parameter '" + parameter_name + "'")

    # -----------------------------------------------------------------

    def write_animation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the animation ...")

        # Determine the path to the new animation
        path = fs.join(self.visualisation_path, time.unique_name("sedfitter") + ".gif")

        # Write the animation
        self.animation.saveto(path)

# -----------------------------------------------------------------

#if False:
#    "THIS DOES NOT WORK!"
#    "Need to interpolate somehow because the 3 parameter vectors have different lengths..."
#    fig = plt.figure(figsize=(10, 4))
#    plotContours(params1, MdustProb, FyoungProb, FionizedProb)
#    fig.savefig(outpath + "plotChi2Contours.pdf", format='pdf')

# -----------------------------------------------------------------

# PLOT CONTOURS: http://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

# -----------------------------------------------------------------

def plotContours(params, MdustProb, FyoungProb, FionizedProb):

    MdustScale     = 1.e7 # in 1e7 Msun
    LyoungScale    = 3.846e26 / (1.425e21*0.153) * 1e8 # in 1e8 Lsun
    LionizingScale = 3.53e14 # in Msun/yr

    MdustBestFit    = 69079833.3333     / MdustScale
    FyoungBestFit   = 1.69488344353e+15 / LyoungScale
    FionizedBestFit = 3.53e+14          / LionizingScale

    p1 = MdustProb[None,:]*FionizedProb[:,None]
    p2 = FyoungProb[None,:]*FionizedProb[:,None]
    p3 = MdustProb[None,:]*FyoungProb[:,None]

    locplot = [[0.07,0.15,0.305,0.81],[0.375,0.15,0.305,0.81],[0.68,0.15,0.305,0.81]]

    fig_a = plt.axes(locplot[0])
    fig_a.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_a.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0]) / MdustScale
    y = np.array(params[2]) / LionizingScale
    fig_a.imshow(p1, cmap='gray', interpolation=None,
               origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_a.set_aspect('auto')

    fig_b = plt.axes(locplot[1])
    fig_b.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_b.set_xlabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    x = np.array(params[1]) / LyoungScale
    y = np.array(params[2]) / LionizingScale
    fig_b.imshow(p2, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_b.set_aspect('auto')

    fig_c = plt.axes(locplot[2])
    fig_c.set_ylabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    fig_c.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0])/MdustScale
    y = np.array(params[1])/LyoungScale
    fig_c.imshow(p3, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_c.set_aspect('auto')

# -----------------------------------------------------------------
