#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.select Contains functions to select a model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ..fitting.run import get_ski_file_for_simulation
from ..config.parameters import parameter_descriptions, default_units, parsing_types_for_parameter_types
from ..config.parameters import types as parameter_types
from ...core.basics.configuration import prompt_string, prompt_string_list, prompt_variable
from ...core.simulation.skifile import SkiFile
from ...core.tools import introspection

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_modeling_ski_templates_path(), "labeled_template.ski")

# -----------------------------------------------------------------

def select_from_model_suite(model_suite, adapt=True, name=None):

    """
    This function ...
    :param model_suite:
    :param adapt:
    :param name:
    :return:
    """

    # Ask for the model name
    if name is None: model_name = prompt_string("model", "name of the model", choices=model_suite.model_names, required=True)
    elif name in model_suite.model_names: model_name = name
    else: raise ValueError("Model name '" + name + "' does not exist")

    # Load the labeled ski template file
    ski = SkiFile(template_ski_path)
    #labels_before = ski.labels

    # Get paths for each label
    label_paths = ski.get_labels_and_paths()
    #print(label_paths["position_angle"])

    # Load the model
    definition = model_suite.get_model_definition(model_name)

    # Add the components to the ski file and to the input map paths dictionary
    input_paths = dict()
    model_suite.add_model_components(model_name, ski, input_paths)

    # Re-add the labels
    ski.add_labels(label_paths, allow_change=False)

    # Set the correct values for the instruments from labeled values elsewhere in the ski file
    ski.fix_labels("instrumentSystem")

    #labels_after = ski.labels
    # Check
    #if not sequences.same_contents(labels_before, labels_after):
    #    labels_string = ", ".join(sequences.difference(labels_before, labels_after))
    #    log.warning("The parameter labels '" + labels_string + "' have been lost in the ski file after setting the components")

    # Get the parameter values
    if adapt: parameter_values = prompt_parameters(ski) # only adapted parameters
    else: parameter_values = get_default_parameters(ski) # all parameters

    # Return
    return model_name, ski, definition, input_paths, parameter_values

# -----------------------------------------------------------------

def get_default_parameters(ski):

    """
    This function ...
    :param ski:
    :return:
    """

    # Get all parameter labels
    parameter_labels = parameter_descriptions.keys()

    # Get the default parameter values
    default_values = dict()
    for label in parameter_labels:

        # Get the value from the ski file
        value = ski.get_labeled_value(label)

        # Cannot find the labeled value
        if value is None:
            log.warning("Could not find the parameter with label '" + label + "' in the ski file")
            continue

        # Set the default value
        default_values[label] = value

    # Return the default values
    return default_values

# -----------------------------------------------------------------

def prompt_parameters(ski):

    """
    This function ...
    :param ski:
    :return:
    """

    # Get parameter labels to adapt
    parameter_labels = prompt_string_list("parameters", "names of the parameters to adapt the value", choices=parameter_descriptions)

    # Get the default parameter values
    default_values = dict()
    for label in parameter_labels:

        # Get the value from the ski file
        value = ski.get_labeled_value(label)

        # Set the default value
        default_values[label] = value

    # Prompt for parameter values
    parameter_values = dict()
    for label in parameter_labels:

        # Get the parameter type
        parameter_type = parameter_types[label]

        # Get the parsing type
        ptype = parsing_types_for_parameter_types[parameter_type]

        # Get the default unit
        unit = default_units[label]

        # Ask for the value
        value = prompt_variable(label, ptype, "value for the " + parameter_descriptions[label], default=default_values[label], required=False) # can be optional

        # Set the value
        parameter_values[label] = value.to(unit)

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def select_from_fitting_context(fitting_context):

    """
    This function ...
    :param fitting_context:
    :return:
    """

    modeling_path = fitting_context.modeling_path

    # Prompt for fitting run
    run_id = prompt_string("fitting_run", "name of the fitting run", choices=fitting_context.fitting_run_names)

    # Load the fitting run
    fitting_run = fitting_context.load_fitting_run(run_id)

    # Prompt for the generation
    generation_name = prompt_generation(fitting_run)

    # Get the parameter values
    simulation_name, parameter_values, chi_squared = get_parameters(fitting_run, generation_name)

    # Load the ski file
    if simulation_name is not None: ski = get_ski_file_for_simulation(modeling_path, run_id, generation_name, simulation_name)
    else:

        ski = fitting_run.ski_template

        # Set parameter values (ALTHOUGH PROBABLY UNNECESSARY)
        ski.set_labeled_values(parameter_values)

    # Load the input paths
    input_map_paths = fitting_run.input_map_paths
    input_paths = input_map_paths

    # Return
    return run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values

# -----------------------------------------------------------------

def prompt_generation(fitting_run):

    """
    This function ...
    :param fitting_run:
    :return:
    """

    # Inform the user
    log.info("Prompting for the generation ...")

    # There are finished generations
    if fitting_run.has_finished_generations:

        # Set the default option for the generation name
        last_generation_name = fitting_run.last_finished_generation
        if last_generation_name is None: last_generation_name = "first_guess"

        # Set the choices for the generationn name
        generation_names = ["first_guess"] + fitting_run.finished_generations

        # Get the generation name
        generation_name = prompt_string("generation", "name of the (finished) generation", default=last_generation_name, choices=generation_names)

    # No finished generations: use first guess
    else:

        # Give warning
        log.warning("There are no finished generations for fitting run '" + fitting_run.name + "'. Using the first guess model values.")

        # Select first guess (basically the model definition)
        generation_name = "first_guess"

    # Set the generation name to None
    if generation_name == "first_guess": generation_name = None

    # Return the genreation name
    return generation_name

# -----------------------------------------------------------------

def get_parameters(fitting_run, generation_name):

    """
    This function ...
    :param fitting_run:
    :param generation_name:
    :return:
    """

    # Inform the user
    log.info("Getting the parameter values from the fitting run ...")

    # If the first guess model should be used
    if generation_name is None:

        # Inform the user
        log.info("Using the parameter values from the initial guess model ...")

        # Get the parameter values
        parameter_values = fitting_run.first_guess_parameter_values
        simulation_name = None
        chi_squared = None

    # If the best simulation of a generation has to be used
    else:

        # Inform the user
        log.info("Using the parameter values from the best model in the '" + generation_name + "' generation ...")

        # Get the parameter values
        #parameter_values = fitting_run.best_parameter_values_for_generation(generation_name)
        simulation_name, parameter_values, chi_squared = fitting_run.best_simulation_name_parameter_values_and_chi_squared_for_generation(generation_name, only_finished=True)

    # Return the parameter values
    return simulation_name, parameter_values, chi_squared

# -----------------------------------------------------------------

def select_from_analysis_context(analysis_context):

    """
    This function ...
    :param analysis_context:
    :return:
    """

    modeling_path = analysis_context.modeling_path

    # Prompt for analysis run
    analysis_run_name = prompt_string("analysis_run", "name of the analysis run", choices=analysis_context.all_analysis_run_names) # cached and local

    # TODO: this function still has to be further implemented

# -----------------------------------------------------------------
