#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.generation_refit Refit the models (calculate the chis squared) for the simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.fitting.modelanalyser import SEDFitModelAnalyser
from pts.core.basics.log import log
from pts.core.basics.distribution import Distribution
from pts.core.plot.distribution import plot_distribution
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.sedfitting import chi_squared_table_to_probabilities

# -----------------------------------------------------------------

# Load the modeling environment and the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Fitting run
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation
definition.add_required("generation", "string", "name of the generation")

# Simulation names
definition.add_positional_optional("simulations", "string_list", "names of the simulations to refit")

# Plotting
definition.add_flag("plot_terms", "plot the chi squared terms", False)
definition.add_flag("plot_ranks", "plot chi squared as a function of rank")
definition.add_flag("plot_distribution", "plot distribution of chi squared values")
definition.add_optional("max_chisquared", "positive_real", "maximum chi squared to show")
definition.add_optional("min_probability", "positive_real", "minimum probability to show")
definition.add_flag("save_plots", "save the plots in appropriate directories")

# Create the configuration
config = parse_arguments("generation_refit", definition, "Refit the models (calculate the chis squared) for the simulations of a generation")

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.run)
generation = fitting_run.get_generation(config.generation)

# Set the simulation names
if config.simulations is not None: simulation_names = config.simulations
else: simulation_names = generation.simulation_names
nsimulations = len(simulation_names)

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

nzeros = 0

# -----------------------------------------------------------------

# Loop over the simulations, run the SED fit model analyser
for simulation_name in simulation_names:

    # Debugging
    log.debug("Refitting the '" + simulation_name + "' simulation ...")

    # Get the simulation object
    simulation = generation.get_simulation_or_basic(simulation_name)

    # Get simulation misc directory
    misc_path = generation.get_simulation_misc_path(simulation_name)

    # Get the mock sed
    mock_sed = generation.get_simulation_mock_sed(simulation_name)

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
    if config.plot_terms:

        # Determine plot path
        if config.save_plots: path = fs.join(misc_path, "chi_squared_terms.pdf")
        else: path = None

        # Create distribution over the wavelength range
        distribution = Distribution.from_columns("Wavelength", differences.wavelengths(add_unit=False), differences.chi_squared_terms(), unit="micron", y_name="Weighed squared difference")
        plot_distribution(distribution, logscale=True, statistics=False, color="xkcd:orange", path=path)

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

# -----------------------------------------------------------------

# Show number of zeros
if nzeros > 5: log.warning(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probabilities of zero")
else: log.debug(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probability of zero")

# -----------------------------------------------------------------

# Get the parameters table
parameters_table = generation.parameters_table

# -----------------------------------------------------------------

# Create the probabilities table
probabilities_table = chi_squared_table_to_probabilities(chi_squared, parameters_table)
probabilities = probabilities_table.probabilities
#nsimulations = len(probabilities)

# Get the number of zero probabilities
nzeros = probabilities_table.nzeros
#print(nzeros)

# -----------------------------------------------------------------

# Plot ranks?
if config.plot_ranks:

    # Determine the plot path
    if config.save_plots: path = fs.join(generation.path, "chi_squared_simulations.pdf")
    else: path = None

    # Plot histogram of the chi squared values w.r.t. simulation rank
    distribution = Distribution.by_rank("Simulation rank", chi_squared.chi_squared_values, y_name="Chi squared")
    plot_distribution(distribution, statistics=False, path=path)

    # Determine the plot path
    if config.save_plots: path = fs.join(generation.path, "probability_simulations.pdf")
    else: path = None

    # Plot histogram of the probability values w.r.t. simulation rank
    distribution = Distribution.by_rank("Simulation rank", probabilities, y_name="Probability")
    plot_distribution(distribution, statistics=False, color="red", logfrequency=True, path=path)

# -----------------------------------------------------------------

# Plot distribution?
if config.plot_distribution:

    # Determine the plot path
    if config.save_plots: path = fs.join(generation.path, "chi_squared_distribution.pdf")
    else: path = None

    # Make distribution of chi squared values and plot histogram
    distribution2 = Distribution.from_values("Chi squared", chi_squared.chi_squared_values, nbins=50, clip_above=config.max_chisquared, density=False)
    plot_distribution(distribution2, path=path)

    # Determine the plot path
    if config.save_plots: path = fs.join(generation.path, "probability_distribution.pdf")
    else: path = None

    # Make distribution of probability values and plot histogram
    distribution2 = Distribution.from_values("Probability", probabilities, nbins=50, ignore_value=0, clip_below=config.min_probability, density=False)
    plot_distribution(distribution2, color="red", path=path)

# -----------------------------------------------------------------
