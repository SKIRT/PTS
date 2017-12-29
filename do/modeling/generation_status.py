#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_status View the status of the simulations of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.distribution import Distribution
from pts.core.tools import formatting as fmt
from pts.core.tools import strings
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Extra
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_flag("plot_memory", "plot memory usage")
definition.add_flag("plot_chisquared", "plot chi squared")

# Get configuration
config = parse_arguments("generation_status", definition, "View the status of the simulations of a certain generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Get the chi squared table
chi_squared = generation.chi_squared_table

# -----------------------------------------------------------------

# Check
if not generation.has_assignment_table: raise RuntimeError("No assignment for this generation")

# -----------------------------------------------------------------

# Get timing and memory table
timing = fitting_run.timing_table
memory = fitting_run.memory_table

# -----------------------------------------------------------------

# Get number of simulations
nsimulations = generation.nsimulations
nretrieved = generation.nretrieved_simulations
nanalysed = generation.nanalysed_simulations
#print(nsimulations, nretrieved, nanalysed)
fraction_analysed = float(nanalysed) / nsimulations

print("")
print("Total number of simulations: " + str(nsimulations))
print("Number of analysed simulations: " + str(nanalysed) + " (" + tostr(fraction_analysed*100, round=True, ndigits=2) + "%)")
print("")

# -----------------------------------------------------------------

total_times = []
setup_times = []
stellar_times = []
spectra_times = []
dust_times = []
writing_times = []
waiting_times = []
communication_times = []
intermediate_times = []

# -----------------------------------------------------------------

# Get the simulations
for simulation in generation.simulations:

    simulation_name = simulation.name

    if simulation.analysed:

        # Get chi squared
        chisq = chi_squared.chi_squared_for(simulation_name)
        ndecimal = 1
        ndigits = 7
        print(" - " + fmt.green + simulation_name + ": " + strings.number(chisq, ndecimal, ndigits, fill=" "))

    else: print(" - " + fmt.red + simulation_name + fmt.reset)

    # Get timing
    if timing.has_simulation(simulation_name):

        #parameters = timing.ski_parameters_for_simulation(simulation_name, which="all")
        #print(parameters)

        # Get timings
        total_time = timing.total_time_for_simulation(simulation_name)
        setup_time = timing.setup_time_for_simulation(simulation_name)
        stellar_time = timing.stellar_emission_time_for_simulation(simulation_name)
        spectra_time = timing.spectra_calculation_time_for_simulation(simulation_name)
        dust_time = timing.dust_emission_time_for_simulation(simulation_name)
        writing_time = timing.writing_time_for_simulation(simulation_name)
        waiting_time = timing.waiting_time_for_simulation(simulation_name)
        communication_time = timing.communication_time_for_simulation(simulation_name)
        intermediate_time = timing.intermediate_time_for_simulation(simulation_name)

        # Add timings
        total_times.append(total_time)
        setup_times.append(setup_time)
        stellar_times.append(stellar_time)
        spectra_times.append(spectra_time)
        dust_times.append(dust_time)
        writing_times.append(writing_time)
        waiting_times.append(waiting_time)
        communication_times.append(communication_time)
        intermediate_times.append(intermediate_time)

# -----------------------------------------------------------------

# Plot?
if config.plot_runtimes:

    total = Distribution.from_values(total_times, unit="min")
    setup = Distribution.from_values(setup_times, unit="min")
    stellar = Distribution.from_values(stellar_times, unit="min")
    spectra = Distribution.from_values(spectra_times, unit="min")
    dust = Distribution.from_values(dust_times, unit="min")
    writing = Distribution.from_values(writing_times, unit="min")
    waiting = Distribution.from_values(waiting_times, unit="min")
    communication = Distribution.from_values(communication_times, unit="min")
    intermediate = Distribution.from_values(intermediate_times, unit="min")

    total.plot(title="Total")
    setup.plot(title="Setup")
    stellar.plot(title="Stellar")
    spectra.plot(title="Spectra")
    dust.plot(title="Dust")
    writing.plot(title="Writing")
    waiting.plot(title="Waiting")
    communication.plot(title="Communication")
    intermediate.plot(title="Intermediate")

# -----------------------------------------------------------------