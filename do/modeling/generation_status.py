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
from pts.core.simulation.remote import SKIRTRemote
from pts.core.basics.log import log
from pts.core.tools import numbers
from pts.core.plot.distribution import plot_distribution

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

# Plotting
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_flag("plot_memory", "plot memory usage")
definition.add_flag("plot_chisquared", "plot chi squared")

# Show
definition.add_flag("show_runtimes", "show runtimes")
definition.add_flag("show_memory", "show memory")

# Show parameters
definition.add_flag("parameters", "show the parameter values")
definition.add_flag("extra", "show extra info")

# Flags
definition.add_flag("offline", "offline modus")
definition.add_flag("lazy", "lazy modus")

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

# Get the maximum chi squared value
min_chi_squared = chi_squared.best_chi_squared
max_chi_squared = chi_squared.worst_chi_squared
max_chi_squared_magnitude = numbers.order_of_magnitude(max_chi_squared)

# Set number of digits for chi squared values
chisq_ndecimal = 1
chisq_ndigits = max_chi_squared_magnitude + 1 + chisq_ndecimal

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

# Determine number of retrieved and analysed simulations
if config.lazy:
    nretrieved = generation.nretrieved_simulations_basic
    nanalysed = generation.nanalysed_simulations_basic
else:
    nretrieved = generation.nretrieved_simulations
    nanalysed = generation.nanalysed_simulations

# -----------------------------------------------------------------

# Determine fraction analysed
fraction_retrieved = float(nretrieved) / nsimulations
fraction_analysed = float(nanalysed) / nsimulations

# Show
print("")
print(fmt.bold + "Total number of simulations: " + fmt.reset + str(nsimulations))
print(fmt.bold + "Number of retrieved simulations: " + fmt.reset + str(nretrieved) + " (" + tostr(fraction_retrieved*100, round=True, ndigits=2) + "%)")
print(fmt.bold + "Number of analysed simulations: " + fmt.reset + str(nanalysed) + " (" + tostr(fraction_analysed*100, round=True, ndigits=2) + "%)")
print("")

# -----------------------------------------------------------------

print(fmt.bold + "Best chi squared: " + fmt.reset + tostr(min_chi_squared))
print(fmt.bold + "Worst chi squared: " + fmt.reset + tostr(max_chi_squared))
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

total_memories = []
setup_memories = []
stellar_memories = []
spectra_memories = []
dust_memories = []
writing_memories = []

# -----------------------------------------------------------------

remotes = dict()
states = dict()
if not config.offline:

    # Inform the user
    log.info("Loading the remotes ...")

    # Create the remote instances
    for host_id in generation.host_ids: remotes[host_id] = SKIRTRemote(host_id=host_id)

    # Get screen states
    for host_id in remotes: states[host_id] = remotes[host_id].screen_states()

# -----------------------------------------------------------------

# Get the parameters table
parameters = generation.parameters_table

# Get parameter units
parameter_units = [parameters.unit_for(label) for label in fitting_run.free_parameter_labels]

# -----------------------------------------------------------------

# Number of digits for parameter values
parameters_ndigits = 3

# -----------------------------------------------------------------

# Print header
if config.parameters:

    if config.extra: nspaces = 52
    else: nspaces = 40

    # Print
    print(" " * nspaces + "\t" + "\t".join(fitting_run.free_parameter_labels))
    print(" " * nspaces + "\t" + "\t".join(["[" + tostr(unit) + "]   " for unit in parameter_units]))

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get the simulation
    if config.lazy or not generation.has_simulation(simulation_name):

        # No simulation object
        simulation = None

        # Get host ID and simulation ID from generation assignment table
        host_id = generation.get_host_id(simulation_name)
        simulation_id = generation.get_simulation_id(simulation_name)

    # Simulation file is not present anymore
    else:

        # Load the simulation
        simulation = generation.get_simulation(simulation_name)

        # Get properties
        host_id = simulation.host_id
        simulation_id = simulation.id

    # Get the parameter values
    if config.parameters: parameter_values = parameters.parameter_values_for_simulation(simulation_name)
    else: parameter_values = None

    # Get parameter string
    if parameter_values is not None:
        parameter_strings = [tostr(parameter_values[label].value, scientific=True, decimal_places=3) + "  " for label in fitting_run.free_parameter_labels]
        parameters_string = "\t".join(parameter_strings)
    else: parameters_string = ""

    # Get extra info
    if config.extra:
        ndigits = 3
        id_string = strings.integer(simulation_id, ndigits)
        host_string = strings.to_length(host_id, 5)
        extra_string = " (" + host_string + " " + id_string + ")"
    else: extra_string = ""

    # No simulation object
    if simulation is None:

        # Has chi squared
        if generation.is_analysed(simulation_name):

            chisq = chi_squared.chi_squared_for(simulation_name)
            print(" - " + fmt.green + simulation_name + extra_string + ": " + strings.number(chisq, chisq_ndecimal, chisq_ndigits, fill=" ") + "\t" + parameters_string + fmt.reset)

        # Has miscellaneous output, but no chi squared
        elif generation.has_misc_output(simulation_name): print(" - " + fmt.yellow + simulation_name + extra_string + ": no chisq" + "\t" + parameters_string + fmt.reset)

        # Has plotting output
        elif generation.has_plotting_output(simulation_name): print(" - " + fmt.yellow + simulation_name + extra_string + ": plotted" + "\t" + parameters_string + fmt.reset)

        # Has extraction output
        elif generation.has_extraction_output(simulation_name): print(" - " + fmt.yellow + simulation_name + extra_string + ": extracted" + "\t" + parameters_string + fmt.reset)

        # Has simulation output
        elif generation.is_retrieved(simulation_name): print(" - " + fmt.yellow + simulation_name + extra_string + ": retrieved" + "\t" + parameters_string + fmt.reset)

        # No simulation output
        else: print(" - " + fmt.red + simulation_name + extra_string + ": unknown" + "\t"  + parameters_string + fmt.reset)

    # Simulation object
    else:

        # Already analysed
        if simulation.analysed:

            # Get chi squared
            chisq = chi_squared.chi_squared_for(simulation_name)
            print(" - " + fmt.green + simulation_name + extra_string +  ": " + strings.number(chisq, chisq_ndecimal, chisq_ndigits, fill=" ") + "\t" + parameters_string + fmt.reset)

        elif simulation.retrieved: print(" - " + fmt.yellow + simulation_name + ": not analysed" + fmt.reset)
        else:

            # Not yet retrieved, what is the status?
            if host_id in remotes: simulation_status = remotes[host_id].get_simulation_status(simulation, screen_states=states[host_id])
            else: simulation_status = " unknown"

            # Show
            if simulation_status == "finished": print(" - " + fmt.yellow + simulation_name + extra_string + ": " + simulation_status + "\t" + parameters_string + fmt.reset)
            else: print(" - " + fmt.red + simulation_name + extra_string + ": " + simulation_status + "\t" + parameters_string + fmt.reset)

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

    # Get memory
    if memory.has_simulation(simulation_name):

        # parameters = memory.ski_parameters_for_simulation(simulation_name, which="all")
        # print(parameters)

        # Get memory info
        total_memory = memory.total_memory_for_simulation(simulation_name)
        setup_memory = memory.setup_memory_for_simulation(simulation_name)
        stellar_memory = memory.stellar_memory_for_simulation(simulation_name)
        spectra_memory = memory.spectra_memory_for_simulation(simulation_name)
        dust_memory = memory.dust_memory_for_simulation(simulation_name)
        writing_memory = memory.writing_memory_for_simulation(simulation_name)

        # Add memory info
        total_memories.append(total_memory)
        setup_memories.append(setup_memory)
        stellar_memories.append(stellar_memory)
        spectra_memories.append(spectra_memory)
        dust_memories.append(dust_memory)
        writing_memories.append(writing_memory)

# -----------------------------------------------------------------

# Make scalar
total_times = [time.to("min").value for time in total_times]
setup_times = [time.to("min").value for time in setup_times]
stellar_times = [time.to("min").value for time in stellar_times]
spectra_times = [time.to("min").value for time in spectra_times]
dust_times = [time.to("min").value for time in dust_times]
writing_times = [time.to("min").value for time in writing_times]
waiting_times = [time.to("min").value for time in waiting_times]
communication_times = [time.to("min").value for time in communication_times]
intermediate_times = [time.to("min").value for time in intermediate_times]

# Make scalar
total_memories = [memory.to("GB").value for memory in total_memories]
setup_memories = [memory.to("GB").value for memory in setup_memories]
stellar_memories = [memory.to("GB").value for memory in stellar_memories]
spectra_memories = [memory.to("GB").value for memory in spectra_memories]
dust_memories = [memory.to("GB").value for memory in dust_memories]
writing_memories = [memory.to("GB").value for memory in writing_memories]

# -----------------------------------------------------------------

# Get number of measurements
ntotal_times = len(total_times)
nsetup_times = len(setup_times)
nstellar_times = len(stellar_times)
nspectra_times = len(spectra_times)
ndust_times = len(dust_times)
nwriting_times = len(writing_times)
nwaiting_times = len(waiting_times)
ncommunication_times = len(communication_times)
nintermediate_times = len(intermediate_times)
ntotal_memories = len(total_memories)
nsetup_memories = len(setup_memories)
nstellar_memories = len(stellar_memories)
nspectra_memories = len(spectra_memories)
ndust_memories = len(dust_memories)
nwriting_memories = len(writing_memories)

# -----------------------------------------------------------------

# Sigma clip runtimes
total_times, total_time_noutliers = numbers.sigma_clip(total_times, return_nmasked=True)
setup_times, setup_time_noutliers = numbers.sigma_clip(setup_times, return_nmasked=True)
stellar_times, stellar_time_noutliers = numbers.sigma_clip(stellar_times, return_nmasked=True)
spectra_times, spectra_time_noutliers = numbers.sigma_clip(spectra_times, return_nmasked=True)
dust_times, dust_time_noutliers = numbers.sigma_clip(dust_times, return_nmasked=True)
writing_times, writing_time_noutliers = numbers.sigma_clip(writing_times, return_nmasked=True)
waiting_times, waiting_time_noutliers = numbers.sigma_clip(waiting_times, return_nmasked=True)
communication_times, communication_time_noutliers = numbers.sigma_clip(communication_times, return_nmasked=True)
intermediate_times, intermediate_time_noutliers = numbers.sigma_clip(intermediate_times, return_nmasked=True)

# Sigma clip memory usages
total_memories, total_memory_noutliers = numbers.sigma_clip(total_memories, return_nmasked=True)
setup_memories, setup_memory_noutliers = numbers.sigma_clip(setup_memories, return_nmasked=True)
stellar_memories, stellar_memory_noutliers = numbers.sigma_clip(stellar_memories, return_nmasked=True)
spectra_memories, spectra_memory_noutliers = numbers.sigma_clip(spectra_memories, return_nmasked=True)
dust_memories, dust_memory_noutliers = numbers.sigma_clip(dust_memories, return_nmasked=True)
writing_memories, writing_memory_noutliers = numbers.sigma_clip(writing_memories, return_nmasked=True)

# -----------------------------------------------------------------

# Show?
if config.show_runtimes:

    # Total
    total = numbers.arithmetic_mean(*total_times)
    total_err = numbers.standard_deviation(*total_times, mean=total)

    # Setup
    setup = numbers.arithmetic_mean(*setup_times)
    setup_err = numbers.standard_deviation(*setup_times, mean=setup)

    # Stellar
    stellar = numbers.arithmetic_mean(*stellar_times)
    stellar_err = numbers.standard_deviation(*stellar_times, mean=stellar)

    # Spectra
    spectra = numbers.arithmetic_mean(*spectra_times)
    spectra_err = numbers.standard_deviation(*spectra_times, mean=spectra)

    # Dust
    dust = numbers.arithmetic_mean(*dust_times)
    dust_err = numbers.standard_deviation(*dust_times, mean=dust)

    # Writing
    writing = numbers.arithmetic_mean(*writing_times)
    writing_err = numbers.standard_deviation(*writing_times, mean=writing)

    # Waiting
    waiting = numbers.arithmetic_mean(*waiting_times)
    waiting_err = numbers.standard_deviation(*waiting_times, mean=waiting)

    # Communication
    communication = numbers.arithmetic_mean(*communication_times)
    communication_err = numbers.standard_deviation(*communication_times, mean=communication)

    # Intermediate
    intermediate = numbers.arithmetic_mean(*intermediate_times)
    intermediate_err = numbers.standard_deviation(*intermediate_times, mean=intermediate)

    # Show
    print("")
    print(fmt.bold + "Runtimes:" + fmt.reset)
    print("")
    print(" - Total time: (" + tostr(total, round=True, ndigits=3) + " ± " + tostr(total_err, round=True, ndigits=3) + ") minutes [" + str(total_time_noutliers) + " outliers out of " + str(ntotal_times) + " data points]")
    print(" - Setup time: (" + tostr(setup, round=True, ndigits=3) + " ± " + tostr(setup_err, round=True, ndigits=3) + ") minutes [" + str(setup_time_noutliers) + " outliers out of " + str(nsetup_times) + " data points]")
    print(" - Stellar time: (" + tostr(stellar, round=True, ndigits=3) + " ± " + tostr(stellar_err, round=True, ndigits=3) + ") minutes [" + str(stellar_time_noutliers) + " outliers out of " + str(nstellar_times) + " data points]")
    print(" - Spectra time: (" + tostr(spectra, round=True, ndigits=3) + " ± " + tostr(spectra_err, round=True, ndigits=3) + ") minutes [" + str(spectra_time_noutliers) + " outliers out of " + str(nspectra_times) + " data points]")
    print(" - Dust time: (" + tostr(dust, round=True, ndigits=3) + " ± " + tostr(dust_err, round=True, ndigits=3) + ") minutes [" + str(dust_time_noutliers) + " outliers out of " + str(ndust_times) + " data points]")
    print(" - Writing time: (" + tostr(writing, round=True, ndigits=3) + " ± " + tostr(writing_err, round=True, ndigits=3) + ") minutes [" + str(writing_time_noutliers) + " outliers out of " + str(nwriting_times) + " data points]")
    print(" - Waiting time: (" + tostr(waiting, round=True, ndigits=3) + " ± " + tostr(waiting_err, round=True, ndigits=3) + ") minutes [" + str(waiting_time_noutliers) + " outliers out of " + str(nwaiting_times) + " data points]")
    print(" - Communication time: (" + tostr(communication, round=True, ndigits=3) + " ± " + tostr(communication_err, round=True, ndigits=3) + ") minutes [" + str(communication_time_noutliers) + " outliers out of " + str(ncommunication_times) + " data points]")
    print(" - Intermediate time: (" + tostr(intermediate, round=True, ndigits=3) + " ± " + tostr(intermediate_err, round=True, ndigits=3) + ") minutes [" + str(intermediate_time_noutliers) + " outliers out of " + str(nintermediate_times) + " data points]")

# -----------------------------------------------------------------

# Plot?
if config.plot_runtimes:

    # Create distributions
    total = Distribution.from_values("Runtime", total_times, unit="min")
    setup = Distribution.from_values("Runtime", setup_times, unit="min")
    stellar = Distribution.from_values("Runtime", stellar_times, unit="min")
    spectra = Distribution.from_values("Runtime", spectra_times, unit="min")
    dust = Distribution.from_values("Runtime", dust_times, unit="min")
    writing = Distribution.from_values("Runtime", writing_times, unit="min")
    waiting = Distribution.from_values("Runtime", waiting_times, unit="min")
    communication = Distribution.from_values("Runtime", communication_times, unit="min")
    intermediate = Distribution.from_values("Runtime", intermediate_times, unit="min")

    # Plot
    plot_distribution(total, title="Total")
    plot_distribution(setup, title="Setup")
    plot_distribution(stellar, title="Stellar")
    plot_distribution(spectra, title="Spectra")
    plot_distribution(dust, title="Dust")
    plot_distribution(writing, title="Writing")
    plot_distribution(waiting, title="Waiting")
    plot_distribution(communication, title="Communication")
    plot_distribution(intermediate, title="Intermediate")

# -----------------------------------------------------------------

# Show
if config.show_memory:

    # Total
    total = numbers.arithmetic_mean(*total_memories)
    total_err = numbers.standard_deviation(*total_memories, mean=total)

    # Setup
    setup = numbers.arithmetic_mean(*setup_memories)
    setup_err = numbers.standard_deviation(*setup_memories, mean=setup)

    # Stellar
    stellar = numbers.arithmetic_mean(*stellar_memories)
    stellar_err = numbers.standard_deviation(*stellar_memories, mean=stellar)

    # Spectra
    spectra = numbers.arithmetic_mean(*spectra_memories)
    spectra_err = numbers.standard_deviation(*spectra_memories, mean=spectra)

    # Dust
    dust = numbers.arithmetic_mean(*dust_memories)
    dust_err = numbers.standard_deviation(*dust_memories, mean=dust)

    # Writing
    writing = numbers.arithmetic_mean(*writing_memories)
    writing_err = numbers.standard_deviation(*writing_memories, mean=writing)

    # Show
    print("")
    print(fmt.bold + "Memory usage:" + fmt.reset)
    print("")
    print(" - Total memory: (" + tostr(total, round=True, ndigits=3) + " ± " + tostr(total_err, round=True, ndigits=3) + ") GB [" + str(total_memory_noutliers) + " outliers out of " + str(ntotal_memories) + " data points]")
    print(" - Setup memory: (" + tostr(setup, round=True, ndigits=3) + " ± " + tostr(setup_err, round=True, ndigits=3) + ") GB [" + str(setup_memory_noutliers) + " outliers out of " + str(nsetup_memories) + " data points]")
    print(" - Stellar memory: (" + tostr(stellar, round=True, ndigits=3) + " ± " + tostr(stellar_err, round=True, ndigits=3) + ") GB [" + str(stellar_memory_noutliers) + " outliers out of " + str(nstellar_memories) + " data points]")
    print(" - Spectra memory: (" + tostr(spectra, round=True, ndigits=3) + " ± " + tostr(spectra_err, round=True, ndigits=3) + ") GB [" + str(spectra_memory_noutliers) + " outliers out of " + str(nspectra_memories) + " data points]")
    print(" - Dust memory: (" + tostr(dust, round=True, ndigits=3) + " ± " + tostr(dust_err, round=True, ndigits=3) + ") GB [" + str(dust_memory_noutliers) + " outliers out of " + str(ndust_memories) + " data points]")
    print(" - Writing memory: (" + tostr(writing, round=True, ndigits=3) + " ± " + tostr(writing_err, round=True, ndigits=3) + ") GB [" + str(writing_memory_noutliers) + " outliers out of " + str(nwriting_memories) + " data points]")

# -----------------------------------------------------------------

# Plot?
if config.plot_memory:

    # Create distributions
    total = Distribution.from_values("Memory usage", total_memories, unit="GB")
    setup = Distribution.from_values("Memory usage", setup_memories, unit="GB")
    stellar = Distribution.from_values("Memory usage", stellar_memories, unit="GB")
    spectra = Distribution.from_values("Memory usage", spectra_memories, unit="GB")
    dust = Distribution.from_values("Memory usage", dust_memories, unit="GB")
    writing = Distribution.from_values("Memory usage", writing_memories, unit="GB")

    # Plot
    plot_distribution(total, title="Total")
    plot_distribution(setup, title="Setup")
    plot_distribution(stellar, title="Stellar")
    plot_distribution(spectra, title="Spectra")
    plot_distribution(dust, title="Dust")
    plot_distribution(writing, title="Writing")

# -----------------------------------------------------------------

# Plot chi squared?
if config.plot_chisquared: chi_squared.distribution.plot(title="Chi squared values", xlogscale=True)

# -----------------------------------------------------------------
