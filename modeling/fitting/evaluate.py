#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.evaluate Contains the evaluate function.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.launch.launcher import SKIRTLauncher
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...core.tools.stringify import stringify
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools.filelock import FileLock
from ...core.tools.stringify import tostr
from ...evolve.optimize.parameters import get_parameters_from_genome
from ...core.tools import introspection

# -----------------------------------------------------------------

def get_parameter_unit(label, fitting_run):

    """
    This function ...
    :param label:
    :param fitting_run:
    :return:
    """

    # unit = Unit(fitting_run.parameter_units[label])

    if label in fitting_run.parameter_units:
        if fitting_run.parameter_units[label] is None: return None
        elif fitting_run.parameter_units[label] == "": return None
        else: return fitting_run.parameter_units[label]
    else: return None

# -----------------------------------------------------------------

def get_parameter_values_from_generator(parameters, index, fitting_run):

    """
    This function ...
    :param parameters:
    :param index:
    :param fitting_run:
    :return:
    """

    # Set the parameter values as a dictionary for this individual model
    parameter_values = dict()
    for label in fitting_run.free_parameter_labels:

        # Get the value for this model from the generator and get the unit defined for this parameter
        value = parameters[label][index]

        # Get unit (if any)
        unit = get_parameter_unit(label, fitting_run)

        # Set value with unit
        if unit is not None: parameter_values[label] = value * unit

        # Set dimensionless value
        else: parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_scalar_parameter_values(genome, fitting_run):

    """
    This function ...
    :param genome:
    :param fitting_run:
    :return:
    """

    # Initialize a dictionary for the parameter values
    parameter_values = dict()

    # Loop over all the genes (parameters)
    for index in range(len(genome)):

        # Get the parameter value
        value = genome[index]

        # Get the laebl of the parameter
        label = fitting_run.free_parameter_labels[index]

        # Set tin the dictionary
        parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_parameter_values_for_named_individual(parameters, name, fitting_run):

    """
    This function ...
    :param parameters: 
    :param name:
    :param fitting_run:
    :return: 
    """

    # Set the parameter values as a dictionary for this individual model
    parameter_values = dict()
    for label in fitting_run.free_parameter_labels:

        # Get the value for this model from the generator
        value = parameters[label][name]

        # Get unit (if any)
        unit = get_parameter_unit(label, fitting_run)

        # Set value with unit
        if unit is not None: parameter_values[label] = value * unit

        # Set dimensionless value
        else: parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_parameter_values_from_genome(genome, fitting_run, minima, maxima, nbits, parameter_scales, gray=False):

    """
    This function ...
    :param genome:
    :param fitting_run:
    :param minima:
    :param maxima:
    :param nbits:
    :param parameter_scales:
    :param gray:
    :return:
    """

    # Get the raw parameter values from the genome
    parameters = get_parameters_from_genome(genome, minima, maxima, nbits, parameter_scales, gray=gray)

    values = dict()

    # Loop over the parameters
    for label_index, label in enumerate(fitting_run.free_parameter_labels):

        # Add unit
        if label in fitting_run.parameter_units: value = parameters[label_index] * fitting_run.parameter_units[label]
        else: value = parameters[label_index]

        # Set the value with unit
        values[label] = value

    # Return the parameter values
    return values

# -----------------------------------------------------------------

def make_test_definition(simulation_name, ski, parameter_values, object_name, simulation_input, scientific=False,
                         fancy=False, ndigits=None):

    """
    This function ...
    :param simulation_name:
    :param ski:
    :param parameter_values:
    :param object_name:
    :param simulation_input:
    :param scientific:
    :param fancy:
    :param ndigits:
    :return:
    """

    # Debugging
    log.debug("Adjusting ski file for the following model parameters:")
    for label in parameter_values: log.debug(" - " + label + ": " + tostr(parameter_values[label], scientific=scientific, fancy=fancy, ndigits=ndigits[label]))

    # Set the parameter values in the ski file template
    ski.set_labeled_values(parameter_values)

    # Create a directory for this simulation
    simulation_path = introspection.create_temp_dir(simulation_name)

    # Create an output directory for this simulation
    simulation_output_path = fs.create_directory_in(simulation_path, "out")

    # Put the ski file with adjusted parameters into the simulation directory
    ski_path = fs.join(simulation_path, object_name + ".ski")
    ski.saveto(ski_path)

    # Create the SKIRT simulation definition
    definition = SingleSimulationDefinition(ski_path, simulation_output_path, simulation_input, name=simulation_name)

    # Return the definition
    return definition

# -----------------------------------------------------------------

def prepare_simulation(simulation_name, ski, parameter_values, object_name, simulation_input, generation_path,
                       scientific=False, fancy=False, ndigits=None):

    """
    This function ...
    :param simulation_name:
    :param ski:
    :param parameter_values:
    :param object_name:
    :param simulation_input:
    :param generation_path:
    :param scientific:
    :param fancy:
    :param ndigits:
    :return:
    """

    # Debugging
    log.debug("Adjusting ski file for the following model parameters:")
    for label in parameter_values: log.debug(" - " + label + ": " + tostr(parameter_values[label], scientific=scientific, fancy=fancy, ndigits=ndigits[label]))

    # Set the parameter values in the ski file template
    ski.set_labeled_values(parameter_values)

    # Create a directory for this simulation
    simulation_path = fs.create_directory_in(generation_path, simulation_name)

    # Create an output directory for this simulation
    simulation_output_path = fs.create_directory_in(simulation_path, "out")

    # Put the ski file with adjusted parameters into the simulation directory
    ski_path = fs.join(simulation_path, object_name + ".ski")
    ski.saveto(ski_path)

    # Create the SKIRT simulation definition
    definition = SingleSimulationDefinition(ski_path, simulation_output_path, simulation_input, name=simulation_name)

    # Return the simulation definition
    return definition

# -----------------------------------------------------------------

def evaluate(genome, **kwargs):

    """
    This function ...
    :param genome:
    :return:
    """

    # Get the fitting run
    fitting_run = kwargs.pop("fitting_run")

    # Get the parameter values
    parameter_values = get_parameter_values_from_genome(genome, fitting_run)

    # Generate simulation name
    simulation_name = generate_simulation_name()

    # Prepare simulation directories, ski file, and return the simulation definition
    definition = prepare_simulation(simulation_name, ski_template, parameter_values, object_name, input_paths, generation.path)
    #simulation_name = definition.name

    # Debugging
    log.debug("Launching a simulation to the queue with:")
    log.debug(" - name: " + simulation_name)
    log.debug(" - input: " + stringify(input_paths)[1])
    log.debug(" - ski path: " + definition.ski_path)
    log.debug(" - output path: " + definition.output_path)

    # Create the SKIRT launcher
    launcher = SKIRTLauncher()

    # Run
    launcher.run(definition=definition)

    # Acquire lock
    with FileLock(parameters_table_path):

        # Open the table
        # ...

        # Add an entry to the parameters table
        parameters_table.add_entry(simulation_name, parameter_values)

        # Save the parameters table
        parameters_table.save()

# -----------------------------------------------------------------
