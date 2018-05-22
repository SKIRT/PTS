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
from ...core.plot.distribution import plot_distributions

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

        # The model probabilities (per generation)
        self.model_probabilities = dict()

        # The parameter probabilities (per generation)
        self.parameter_probabilities = dict()

        # All generations parameter probabilities
        self.parameter_probabilities_all = dict()

        # The model parameter table
        self.parameter_tables = dict()

        # The tables with the probability distributions for the different fit parameters
        self.distributions = dict()

        # The animation
        self.animation = None

        # The directory with the probability tables for all finished generations
        self.prob_generations_path = None

        # The dictionary with ...
        #self.prob_generations_table_paths = dict()
        self.prob_generations_paths = dict()

    # -----------------------------------------------------------------

    @property
    def do_create_distributions(self):

        """
        This function ...
        :return:
        """

        return not self.config.per_generation

    # -----------------------------------------------------------------

    @property
    def do_animate(self):

        """
        This function ...
        :return:
        """

        return self.config.visualise

    # -----------------------------------------------------------------

    @property
    def do_writing(self):

        """
        This function ...
        :return:
        """

        return self.config.write

    # -----------------------------------------------------------------

    @property
    def do_plotting(self):

        """
        This function ...
        :return:
        """

        return self.config.plot

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the chi squared values
        self.get_chi_squared()

        # 2. Get the parameters of the best models for each generation
        self.get_best_parameters()

        # 3. Calculate the probabilities
        self.calculate_probabilities()

        # 4. Calculate the probability distributions
        if self.do_create_distributions: self.create_distributions()

        # 5. Make an animation of the fitting procedure
        if self.do_animate: self.animate()

        # 6. Writing
        if self.do_writing: self.write()

        # 7. Plot
        if self.do_plotting: self.plot()

    # -----------------------------------------------------------------

    @property
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_labels

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

        # For each finished generation, determine the path to the probabilities directory
        for generation_name in self.generation_names: self.prob_generations_paths[generation_name] = fs.create_directory_in(self.prob_generations_path, generation_name)

        # Rerun?
        if self.config.rerun_all:
            fs.clear_directory(self.parameters_directory_path)
            fs.clear_directory(self.distributions_directory_path)
            for generation_name in self.generation_names:
                generation_path = self.prob_generations_paths[generation_name]
                fs.clear_directory(generation_path)
        elif self.config.rerun is not None:
            fs.clear_directory(self.parameters_directory_path)
            fs.clear_directory(self.distributions_directory_path)
            generation_path = self.prob_generations_paths[self.config.rerun]
            fs.clear_directory(generation_path)

    # -----------------------------------------------------------------

    @property
    def only_finished(self):

        """
        This function ...
        :return:
        """

        return not self.config.unfinished

    # -----------------------------------------------------------------

    def get_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Calculate the chi squared values
        if self.config.recalculate_chisquared: self.calculate_chi_squared()

        # Load the chi squared values
        else: self.load_chi_squared()

    # -----------------------------------------------------------------

    def load_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the chi squared values ...")

        # Loop over the generations
        for generation_name in self.generation_names: pass

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared values ...")

        # Open the chi squared table
        chi_squared = generation.chi_squared_table

        # -----------------------------------------------------------------

        # Only remove certain simulations
        if config.simulations is not None:

            # Loop over the simulations, remove the chi squared entry
            for simulation_name in simulation_names:

                # Does simulation exist in table?
                if not chi_squared.has_simulation(simulation_name): continue

                # Remove entry from table
                chi_squared.remove_simulation(simulation_name)

        # Remove all simulations
        else: chi_squared.remove_all_simulations()

        # -----------------------------------------------------------------

        # Save the chi squared table
        chi_squared.save()

        nzeros = 0

        # Loop over the simulations, run the SED fit model analyser
        for simulation_name in simulation_names:

            # Debugging
            log.debug("Refitting the '" + simulation_name + "' simulation ...")

            # Get the simulation object
            simulation = generation.get_simulation_or_basic(simulation_name)

            # Get simulation misc directory
            misc_path = generation.get_simulation_misc_path(simulation_name)

            # Get the mock sed
            mock_sed = generation.get_mock_sed(simulation_name)

            # Create the fit model analyser
            analyser = SEDFitModelAnalyser()

            # Set options
            analyser.config.write = False

            # Run the analyser
            analyser.run(simulation=simulation, fitting_run=fitting_run, mock_sed=mock_sed)

            # Get the chi squared
            chisq = analyser.chi_squared

            # Get the differences table
            differences = analyser.differences
            differences.sort()

            # Plot the chi squared terms?
            #if config.plot_terms:

            # Calculate the probability
            probability = np.exp(-0.5 * chisq)
            if probability == 0: nzeros += 1

            # Show chi squared
            log.debug("The chi squared value for simulation '" + simulation_name + "' is " + str(chisq) + " and the probability is " + str(probability))

            # Debugging
            log.debug("Adding to the chi squared table ...")

            # Set the chi squared value
            chi_squared.add_entry(simulation_name, chisq)

            # Save the chi squared table
            chi_squared.save()

            # Debugging
            log.debug("Writing the differences table ...")

            # Set the differences path
            differences_path = generation.get_simulation_misc_differences_path(simulation_name)

            # Save the differences table
            differences.saveto(differences_path)

        # Show number of zeros
        if nzeros > 5: log.warning(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probabilities of zero")
        else: log.debug(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probability of zero")

    # -----------------------------------------------------------------

    def get_best_parameters(self):

        """"
        This function ...
        """

        # Inform the user
        log.info("Getting the parameter values of the best model for the generations (if not already done) ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Check if the generation is already in the best parameters table:
            # NO: Refitting can be performed
            #if generation_name in self.fitting_run.best_parameters_table.generation_names: continue

            # Debugging
            log.debug("Getting the best parameter values for generation '" + generation_name + "' ...")

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

        # Get the probability tables
        self.get_model_probabilities()

        # Calculate the parameter probabilities for each generation separately
        self.calculate_parameter_probabilities()

        # Calculate the combined probabilities for all generations
        if not self.config.per_generation: self.calculate_all_parameter_probabilities()

    # -----------------------------------------------------------------

    def get_model_probabilities_table_path_for_generation(self, generation_name):

        """
        Thisf unction ...
        :param generation_name:
        :return:
        """

        return fs.join(self.prob_generations_paths[generation_name], "models.dat")

    # -----------------------------------------------------------------

    def has_model_probabilities_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        path = self.get_model_probabilities_table_path_for_generation(generation_name)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table_path_for_generation(self, generation_name, parameter_label):

        """
        This function ...
        :param generation_name:
        :param parameter_label:
        :return:
        """

        return fs.join(self.prob_generations_paths[generation_name], parameter_label + ".dat")

    # -----------------------------------------------------------------

    def has_parameter_probabilities_table_path_for_generation(self, generation_name, parameter_label):

        """
        This function ...
        :param generation_name:
        :param parameter_label:
        :return:
        """

        path = self.get_parameter_probabilities_table_path_for_generation(generation_name, parameter_label)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def get_model_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model probabilities for the generations ...")

        # Loop over the finished generations
        for generation_name in self.generation_names:

            # Check whether the probabilities table is already present for this generation
            if self.has_model_probabilities_table_for_generation(generation_name): self.load_model_probabilities_generation(generation_name)

            # Otherwise, calculate the probabilities based on the chi squared table
            else: self.calculate_model_probabilities_generation(generation_name)

    # -----------------------------------------------------------------

    def load_model_probabilities_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Debugging
        log.debug("Loading the model probabilities table for generation '" + generation_name + "' ...")

        # Determine the path
        path = self.get_model_probabilities_table_path_for_generation(generation_name)

        # Load the probabilities table
        probabilities_table = ModelProbabilitiesTable.from_file(path)

        # Add to the dictionary
        self.model_probabilities[generation_name] = probabilities_table

    # -----------------------------------------------------------------

    def get_simulation_names_parameters_and_chi_squared_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Load the parameter table
        parameter_table = self.fitting_run.parameters_table_for_generation(generation_name)

        # Load the chi squared table
        chi_squared_table = self.fitting_run.chi_squared_table_for_generation(generation_name)

        # Sort the table for decreasing chi squared value
        chi_squared_table.sort("Chi squared")
        chi_squared_table.reverse()

        # Get the chi squared values
        chi_squared_values = list(chi_squared_table["Chi squared"])

        # Initialize lists
        simulation_names = []
        parameter_values = []

        # Loop over the simulations
        for i in range(len(chi_squared_table)):

            # Get the simulation name
            simulation_name = chi_squared_table["Simulation name"][i]

            # Get a dictionary with the parameter values for this simulation
            values = parameter_table.parameter_values_for_simulation(simulation_name)

            # Add to lists
            simulation_names.append(simulation_name)
            parameter_values.append(values)

        # Return
        return simulation_names, parameter_values, chi_squared_values

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.parameter_units

    # -----------------------------------------------------------------

    def calculate_model_probabilities_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Debugging
        log.debug("Calculating model probabilities for generation '" + generation_name + "' ...")

        # Get simulation names with parameter values and chi squared values
        simulation_names, parameter_values, chi_squared_values = self.get_simulation_names_parameters_and_chi_squared_for_generation(generation_name)
        nsimulations = len(simulation_names)

        # Convert chi squared values to probabilities
        probabilities_table = chi_squared_to_probabilities(chi_squared_values, simulation_names, parameter_values, self.free_parameter_labels, self.parameter_units)

        # Check the probabilities
        if probabilities_table.all_zero:

            # Give warning
            log.warning("All probabilities for the '" + generation_name + "' generation are zero")

            # Get the minimum chi squared and subtract it from all chi squared values
            min_chi_squared = min(chi_squared_values)
            new_chi_squared_values = [value - min_chi_squared + 1 for value in chi_squared_values]

            # Try converting new chi squared values to probabilities
            probabilities_table = chi_squared_to_probabilities(new_chi_squared_values, simulation_names, parameter_values,
                                                               self.free_parameter_labels, self.parameter_units)

            # Check again
            if probabilities_table.all_zero: raise RuntimeError("All probabilities for the '" + generation_name + "' are still zero")
            if probabilities_table.all_zero_but_one: raise RuntimeError("All probabilities for the '" + generation_name + "' are still zero")
            if probabilities_table.has_zeros: log.warning("Some probabilities for the '" + generation_name + "' generation are zeros (" + str(probabilities_table.nzeros) + " out of " + str(nsimulations) + ")")

        # Some zeros
        elif probabilities_table.has_zeros: log.warning("Some probabilities for the '" + generation_name + "' generation are zeros (" + str(probabilities_table.nzeros) + " out of " + str(nsimulations) + ")")

        # Save the model probabilities table
        table_path = self.get_model_probabilities_table_path_for_generation(generation_name)
        probabilities_table.saveto(table_path)

        # Add to the dictionary
        self.model_probabilities[generation_name] = probabilities_table

    # -----------------------------------------------------------------

    def calculate_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities of the different parameter values for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Initialize dictionary for this generation
            self.parameter_probabilities[generation_name] = dict()

            # Loop over the free parameters
            for label in self.free_parameter_labels:

                # Create a set for the unique values
                unique_values = set()

                # Initialize a ParameterProbabilitiesTable instance for this parameter
                table = ParameterProbabilitiesTable()

                # Loop over the values of this parameter for this generation, and expand the set accordingly
                for value in self.model_probabilities[generation_name][label]: unique_values.add(value)

                # Get a (sorted) list of all the unique values for this parameter
                unique_values = sorted(list(unique_values))

                # Add an entry for each unique parameter value that has been encountered
                for value in unique_values:

                    # Get an array of the probabilities of all the models that have this unique value
                    simulation_indices = self.model_probabilities[generation_name][label] == value
                    #nsimulations_for_value = np.sum(simulation_indices)
                    #individual_probabilities += list(self.model_probabilities[generation_name]["Probability"][simulation_indices])
                    individual_probabilities = self.model_probabilities[generation_name]["Probability"][simulation_indices]

                    # Combine the individual probabilities
                    combined_probability = np.sum(np.asarray(individual_probabilities))

                    # Add an entry to the table
                    table.add_entry(value, combined_probability)

                # Set the table
                self.parameter_probabilities[generation_name][label] = table

    # -----------------------------------------------------------------

    def calculate_all_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities of the different parameter values for all generations ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

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
            self.parameter_probabilities_all[label] = table

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the probability distributions ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Debugging
            log.debug("Creating distribution for the '" + label + "' parameter ...")

            # Convert the probability lists into NumPy arrays and normalize them
            normalized_probabilities = np.array(self.parameter_probabilities_all[label]["Probability"]) / sum(self.parameter_probabilities_all[label]["Probability"])

            # Create the probability distributions for the different parameters
            self.distributions[label] = Distribution.from_probabilities(label, normalized_probabilities, self.parameter_probabilities_all[label]["Value"])

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

    def get_description(self, parameter_label):

        """
        This function ...
        :param parameter_label:
        :return:
        """

        return self.fitting_run.parameter_descriptions[parameter_label]

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

        # Write all generations parameter probabilities
        if not self.config.per_generation: self.write_all_parameter_probabilities()

        # Write the best parameters table
        self.write_best_parameters()

        # Write the probability distributions in table format
        if not self.config.per_generation: self.write_distributions()

        # Write the animated GIF
        if self.config.visualise: self.write_animation()

    # -----------------------------------------------------------------

    def write_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter probabilities for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Loop over the free parameters
            for label in self.free_parameter_labels:

                # Get the parameter probabilities
                probabilities = self.parameter_probabilities[generation_name][label]

                # Determine the path
                path = self.get_parameter_probabilities_table_path_for_generation(generation_name, label)

                # Save the table
                probabilities.saveto(path)

    # -----------------------------------------------------------------

    @property
    def parameters_directory_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.prob_parameters_path

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table_path(self, parameter_label):

        """
        This function ...
        :param parameter_label:
        :return:
        """

        return self.fitting_run.get_parameter_probabilities_path(parameter_label)

    # -----------------------------------------------------------------

    def write_all_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter probabilities from all generations ...")

        # Loop over the probability tables for the different free parameter
        for label in self.fitting_run.free_parameter_labels:

            # Debugging
            log.debug("Writing the parameter probabilities for the " + self.get_description(label) + " ...")

            # Get path
            path = self.get_parameter_probabilities_table_path(label)

            # Save the table
            self.parameter_probabilities_all[label].saveto(path)

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

    @property
    def distributions_directory_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.prob_distributions_path

    # -----------------------------------------------------------------

    def get_distribution_table_path(self, parameter_label):

        """
        This function ...
        :param parameter_label:
        :return:
        """

        return self.fitting_run.get_parameter_distribution_path(parameter_label)

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
            log.debug("Writing the probability distribution of the " + self.get_description(label) + " ...")

            # Get path
            path = self.get_distribution_table_path(label)

            # Write the table of probabilities for this parameter
            self.distributions[label].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the parameter probabilities per generation
        self.plot_probabilities()

        # Plot the probability distributions as histograms
        if not self.config.per_generation: self.plot_distributions()

    # -----------------------------------------------------------------

    def plot_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the parameter probabilities for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Determine the path
            path = fs.join(self.prob_generations_paths[generation_name], "probabilities.pdf")

            # Make distributions
            distributions = dict()
            for label in self.free_parameter_labels:
                distribution = Distribution.from_probabilities(label, self.parameter_probabilities[generation_name][label]["Probability"], self.parameter_probabilities[generation_name][label]["Value"])
                distribution.normalize(method="sum")
                distributions[label] = distribution

            # Plot in different panels
            plot_distributions(distributions, panels=True, extrema=True, frequencies=True, path=path, logscale=True)

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Determine the path
        path = fs.join(self.fitting_run.prob_distributions_path, "distributions.pdf")

        # Plot in different panels
        plot_distributions(self.distributions, panels=True, extrema=True, frequencies=True, path=path, logscale=True)

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

def chi_squared_table_to_probabilities(table, parameter_table):

    """
    This function ...
    :param table:
    :param parameter_table:
    :return:
    """

    # Get the free parameter labels and the parameter units
    parameter_labels = parameter_table.parameter_labels
    parameter_units = parameter_table.parameter_units

    # Sort the table for decreasing chi squared value
    table.sort("Chi squared")
    table.reverse()

    # Get the chi squared values
    chi_squared_values = list(table["Chi squared"])

    # Initialize lists
    simulation_names = []
    parameter_values = []

    # Loop over the simulations
    for i in range(len(table)):

        # Get the simulation name
        simulation_name = table["Simulation name"][i]

        # Get a dictionary with the parameter values for this simulation
        values = parameter_table.parameter_values_for_simulation(simulation_name)

        # Add to lists
        simulation_names.append(simulation_name)
        parameter_values.append(values)

    # Calculate probabilities
    return chi_squared_to_probabilities(chi_squared_values, simulation_names, parameter_values, parameter_labels, parameter_units)

# -----------------------------------------------------------------

def chi_squared_to_probabilities(chi_squared_values, simulation_names, parameter_values, parameter_labels, parameter_units):

    """
    This function ...
    :param chi_squared_values:
    :param simulation_names:
    :param parameter_values:
    :param parameter_labels:
    :param parameter_units:
    :return:
    """

    # Get number of simulations
    nsimulations = len(simulation_names)

    # Calculate the probability for each model
    probabilities = np.exp(-0.5 * np.asarray(chi_squared_values))

    # Create the probabilities table
    probabilities_table = ModelProbabilitiesTable(parameters=parameter_labels, units=parameter_units)

    # Add the entries to the model probabilities table
    for i in range(nsimulations):

        # Get the simulation name
        simulation_name = simulation_names[i]

        # Get a dictionary with the parameter values for this simulation
        values = parameter_values[i]

        # Get the probability
        probability = probabilities[i]

        # Add an entry to the table
        probabilities_table.add_entry(simulation_name, values, probability)

    # Return the probability table
    return probabilities_table

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
