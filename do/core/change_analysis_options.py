#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.change_analysis_options Change certain analysis options for a single or multiple simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.tools import introspection
from pts.do.core.change_simulation_settings import save_simulation, save_simulations, get_common_ptype
from pts.core.launch.options import get_analysis_property_names_and_descriptions, get_analysis_section_names_and_descriptions, get_analysis_property_names_and_descriptions_for_section
from pts.core.tools import sequences
from pts.core.tools.stringify import stringify, tostr
from pts.core.basics.configuration import prompt_proceed, prompt_variable
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("ids", "integer_list", "simulation IDs")
definition.add_positional_optional("matching", "string", "only adapt settings with a name matching this string", suggestions=["remote"])

# -----------------------------------------------------------------

# Select certain properties
definition.add_optional("contains", "string", "only adapt properties containing this string in their name")
definition.add_optional("not_contains", "string", "don't adapt properties containing this string in their name")
definition.add_optional("exact_name", "string", "only adapt properties with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't adapt properties with this exact string as their name")
definition.add_optional("startswith", "string", "only adapt properties whose name starts with this string")
definition.add_optional("endswith", "string", "only adapt properties whose name starts with this string")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("change_analysis_options", definition, description="Change certain analysis options for a simulation")

# -----------------------------------------------------------------

# Check settings
if config.matching is not None:
    if config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
    config.contains = config.matching

# -----------------------------------------------------------------

# No IDs specified?
if config.ids is None: config.ids = introspection.simulation_ids_for_host(config.remote)

# -----------------------------------------------------------------

nids = len(config.ids)
has_single_id = nids == 1

# -----------------------------------------------------------------

def get_analysis_value_for_simulation(simulation, name, section=None):

    """
    This function ...
    :param simulation:
    :param name:
    :param section:
    :return:
    """

    if section is not None: return simulation.analysis[section][name]
    else: return simulation.analysis[name]

# -----------------------------------------------------------------

def get_analysis_values_for_simulations(simulations, name, section=None):

    """
    This function ...
    :param simulations:
    :param name:
    :param section:
    :return:
    """

    values = dict()

    # Loop over the simulations
    for simulation_id in simulations:

        # Get value
        value = get_analysis_value_for_simulation(simulations[simulation_id], name, section=section)

        # Set value
        values[simulation_id] = value

    # Return the values
    return values

# -----------------------------------------------------------------

def set_analysis_value_for_simulation(simulation, name, value, section=None):

    """
    This function ...
    :param simulation:
    :param name:
    :param value:
    :param section:
    :return:
    """

    # Debugging
    log.debug("Changing the value of '" + name + "' to '" + tostr(value) + "' ...")

    # Set the new value
    if section is not None: simulation.analysis[section][name] = value
    else: simulation.analysis[name] = value

# -----------------------------------------------------------------

def set_analysis_value_for_simulations(simulations, name, value, section=None):

    """
    This function ...
    :param simulations:
    :param name:
    :param section:
    :return:
    """

    # Debugging
    log.debug("Changing the value of '" + name + "' to '" + tostr(value) + " for all simulations ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Set value
        set_analysis_value_for_simulation(simulations[simulation_id], name, value, section=section)

# -----------------------------------------------------------------

def set_analysis_value_for_simulations_prompt(simulations, name, value, changed=None, section=None):

    """
    This function ...
    :param simulations:
    :param name:
    :param value:
    :param changed:
    :param section:
    :return:
    """

    # Initialize changed dictinoary
    if changed is None: changed = {simulation_id: False for simulation_id in simulations}

    # Debugging
    log.debug("Changing the values of '" + name + "' for each simulation ...")

    # Loop over the simulations, set the value
    for simulation_id in simulations:

        # Get simulation name
        simulation = simulations[simulation_id]
        simulation_name = simulation.name

        # Get current value
        current_value = get_analysis_value_for_simulation(simulation, name, section=section)

        # Check
        if current_value == value: continue

        # Prompt
        if prompt_proceed("Replace the value of '" + name + "' from '" + tostr(current_value) + "' to '" + tostr(value) + " for simulation '" + simulation_name + "'?"):
            set_analysis_value_for_simulation(simulation, name, value, section=section)
            changed[simulation_id] = True

    # Return the changed dictionary
    return changed

# -----------------------------------------------------------------

# Adapt a single simulation
if has_single_id:

    # Open the simulation object
    single_id = config.ids[0]
    simulation = get_simulation_for_host(config.remote, single_id)

    # Update
    simulation.update_analysis_options()

    # Check whether analysis options are defined
    simulation.analysis.prompt_properties(contains=config.contains, not_contains=config.not_contains, exact_name=config.exact_name, exact_not_name=config.exact_not_name, startswith=config.startswith, endswith=config.endswith)

    # Save the simulation
    save_simulation(simulation)

# -----------------------------------------------------------------

# Adapt multiple simulations
else:

    # Load the simulations and put them in a dictionary
    simulations = OrderedDict()
    for simulation_id in config.ids: simulations[simulation_id] = get_simulation_for_host(config.remote, simulation_id)

    # Update analysis options in each simulation
    for simulation_id in config.ids: simulations[simulation_id].update_analysis_options()

    # Create a dictionary to contain a flag for each simulation that tells whether it has changed
    changed = dict()
    for simulation_id in simulations: changed[simulation_id] = False

    # Get properties
    properties = get_analysis_property_names_and_descriptions()

    # Loop over the properties
    for name in properties:

        # Checks
        if config.contains is not None and config.contains not in name: continue
        if config.not_contains is not None and config.not_contains in name: continue
        if config.exact_name is not None and name != config.exact_name: continue
        if config.exact_not_name is not None and name == config.exact_not_name: continue
        if config.startswith is not None and not name.startswith(config.startswith): continue
        if config.endswith is not None and not name.endswith(config.endswith): continue

        # Get description
        description = properties[name]

        # Get the analysis options for all the simulations
        values = get_analysis_values_for_simulations(simulations, name)

        # Get unique values
        unique_values = sequences.unique_values(values.values())

        # Only one unique value
        if len(unique_values) == 1:

            default = unique_values[0]
            ptype, string = stringify(default)
            choices = None
            suggestions = None

        # Multiple unique values
        else:

            # Prompt to change this property
            change = prompt_proceed("Change the analysis option '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
            if not change: continue
            default = None
            ptype = get_common_ptype(unique_values)
            choices = None
            suggestions = unique_values

        # Ask for the new value
        value = prompt_variable(name, ptype, description, default=default, required=True, choices=choices, suggestions=suggestions)

        # Each simulation had the same value: adapt each simulation's analysis options simultaneously without prompting to proceed
        if len(unique_values) == 1:

            if value != default:
                set_analysis_value_for_simulations(simulations, name, value)
                for simulation_id in simulations: changed[simulation_id] = True

        # Different simulations had different values
        else: changed = set_analysis_value_for_simulations_prompt(simulations, name, value, changed=changed)

    # Get sections
    sections = get_analysis_section_names_and_descriptions()

    # Loop over the sections
    for section_name in sections:

        # Get section description
        section_description = sections[section_name]

        # Debug
        log.debug("Entering section '" + section_name + "' ...")

        # Get properties
        properties = get_analysis_property_names_and_descriptions_for_section(section_name)

        # Loop over the properties in this section
        for name in properties:

            # Checks
            if config.contains is not None and config.contains not in name: continue
            if config.not_contains is not None and config.not_contains in name: continue
            if config.exact_name is not None and name != config.exact_name: continue
            if config.exact_not_name is not None and name == config.exact_not_name: continue
            if config.startswith is not None and not name.startswith(config.startswith): continue
            if config.endswith is not None and not name.endswith(config.endswith): continue

            # Get description
            description = properties[name]

            # Get the analysis options for all the simulations
            values = get_analysis_values_for_simulations(simulations, name, section=section_name)

            # Get unique values
            unique_values = sequences.unique_values(values.values())

            # Only one unique value
            if len(unique_values) == 1:

                default = unique_values[0]
                ptype, string = stringify(default)
                choices = None
                suggestions = None

            # Multiple unique values
            else:

                # Prompt to change this property
                change = prompt_proceed("Change the analysis option '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
                if not change: continue
                default = None
                ptype = get_common_ptype(unique_values)
                choices = None
                suggestions = unique_values

            # Ask for the new value
            value = prompt_variable(name, ptype, description, default=default, required=True, choices=choices, suggestions=suggestions)

            # Each simulation had the same value: adapt each simulation's analysis options simultaneously without prompting to proceed
            if len(unique_values) == 1:

                if value != default:
                    set_analysis_value_for_simulations(simulations, name, value, section=section_name)
                    for simulation_id in simulations: changed[simulation_id] = True

                # Different simulations had different values
                else: changed = set_analysis_value_for_simulations_prompt(simulations, name, value, changed=changed, section=section_name)

    # Save the simulations
    save_simulations(simulations, changed=changed)

# -----------------------------------------------------------------
