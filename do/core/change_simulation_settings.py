#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.change_simulation_settings Change certain settings of a simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_variable, prompt_proceed
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.basics.log import log
from pts.core.tools.stringify import stringify, tostr
from pts.core.tools import sequences
from pts.core.basics.configuration import parent_type
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("ids", "integer_list", "simulation IDs (if none specified, all simulation IDs will be used)")
definition.add_positional_optional("matching", "string", "only adapt settings with a name matching this string")

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
config = parse_arguments("change_simulation_settings", definition, description="Change certain settings of a simulation")

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

properties = OrderedDict()
properties["ski_path"] = "local ski file path"
properties["input_path"] = "local input path"
properties["output_path"] = "local output path"
properties["base_path"] = "local simulation base path"
properties["name"] = "simulation name"
properties["host_id"] = "remote host ID"
properties["cluster_name"] = "remote cluster name"
properties["id"] = "remote host simulation ID"
properties["remote_ski_path"] = "remote ski file path"
properties["remote_simulation_path"] = "remote simulation path"
properties["remote_input_path"] = "remote input path"
properties["remote_output_path"] = "remote output path"
properties["submitted_at"] = "simulation submission timestamp"
properties["retrieve_types"] = "retrieve output file types"
properties["remove_remote_input"] = "remove remote input directory"
properties["remove_remote_output"] = "remove remote output directory"
properties["remove_remote_simulation_directory"] = "remove remote simulation directory"
properties["remove_local_output"] = "remove local output after analysis"
properties["retrieved"] = "retrieved flag"

# -----------------------------------------------------------------

def get_value_for_simulation(simulation, name):

    """
    Thisn function ...
    :param simulation:
    :param name:
    :return:
    """

    # Get the current value
    return getattr(simulation, name)

# -----------------------------------------------------------------

def set_value_for_simulation(simulation, name, value):

    """
    This function ...
    :param simulation:
    :param name:
    :param value:
    :return:
    """

    # Debugging
    log.debug("Changing the value of '" + name + "' to '" + tostr(value) + "' ...")

    # Set the new value
    setattr(simulation, name, value)

# -----------------------------------------------------------------

def get_values_for_simulations(simulations, name):

    """
    This function ...
    :param simulations:
    :param name:
    :return:
    """

    values = dict()

    # Loop over the simulations
    for simulation_id in simulations:

        # Get value
        value = get_value_for_simulation(simulations[simulation_id], name)

        # Set value
        values[simulation_id] = value

    # Return the values
    return values

# -----------------------------------------------------------------

def set_value_for_simulations(simulations, name, value):

    """
    This function ...
    :param simulations:
    :param name:
    :return:
    """

    # Debugging
    log.debug("Changing the value of '" + name + "' to '" + tostr(value) + " for all simulations ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Set value
        set_value_for_simulation(simulations[simulation_id], name, value)

# -----------------------------------------------------------------

def set_value_for_simulations_prompt(simulations, name, value):

    """
    This function ...
    :param simulations:
    :param name:
    :param value:
    :return:
    """

    # Debugging
    log.debug("Changing the values of '" + name + "' for each simulation ...")

    # Loop over the simulations, set the value
    for simulation_id in simulations:

        # Get simulation name
        simulation = simulations[simulation_id]
        simulation_name = simulation.name

        # Get current value
        current_value = get_value_for_simulation(simulation, name)

        # Check
        if current_value == value: continue

        # Prompt
        if prompt_proceed("Replace the value of '" + name + "' from '" + tostr(current_value) + "' to '" + tostr(value) + " for simulation '" + simulation_name + "'?"):
            set_value_for_simulation(simulation, name, value)

# -----------------------------------------------------------------

def get_common_ptype(values):

    """
    Thisj function ...
    :param values:
    :return:
    """

    ptype = None
    ptypes = set()

    for value in values:

        parsetype, val = stringify(value)

        if ptype is None: ptype = parsetype
        elif ptype != parsetype: ptype = "mixed"

        # Add the parse type
        ptypes.add(parsetype)

    ptypes = list(ptypes)

    if len(ptypes) == 1: ptype = ptypes[0]
    elif sequences.all_equal(ptypes): ptype = ptypes[0]
    else:

        # Investigate the different ptypes
        parent_types = [parent_type(type_name) for type_name in ptypes]

        # Check
        for i in range(len(parent_types)):
            if parent_types[i] is None: log.warning("Could not determine the parent type for '" + ptypes[i] + "'. All parent types: " + str(parent_types))
        if sequences.all_equal(parent_types) and parent_types[0] is not None: ptype = parent_types[0]
        elif ptype == "mixed": log.warning("Could not determine a common type for '" + stringify(parent_types)[1] + "'")

    # Return the type
    return ptype

# -----------------------------------------------------------------

# Adapt a single simulation
if has_single_id:

    # Open the simulation object
    single_id = config.ids[0]
    simulation = get_simulation_for_host(config.remote, single_id)

    # Loop over the properties
    for name in properties:

        # Checks
        if config.contains is not None and config.contains not in name: continue
        if config.not_contains is not None and config.not_contains in name: continue
        if config.exact_name is not None and name != config.exact_name: continue
        if config.exact_not_name is not None and name == config.exact_not_name: continue
        if config.startswith is not None and not name.startswith(config.startswith): continue
        if config.endswith is not None and not name.endswith(config.endswith): continue

        # Get the description
        description = properties[name]

        # Get the current value
        default = get_value_for_simulation(simulation, name)

        # There is a current value
        if default is not None: ptype, pstring = stringify(default)

        # No current value
        else: ptype = "any"

        # Check ptype
        if ptype is None: ptype = "any"

        # Ask for the new value
        value = prompt_variable(name, ptype, description, default=default, required=True)
        if default is None and value == "": continue

        # Set the property
        if value != default: set_value_for_simulation(simulation, name, value)

    # Save the simulation
    simulation.save()

# -----------------------------------------------------------------

# Adapt multiple simulations
else:

    # Load the simulations and put them in a dictionary
    simulations = OrderedDict()
    for simulation_id in config.ids: simulations[simulation_id] = get_simulation_for_host(config.remote, simulation_id)

    # Loop over the properties
    for name in properties:

        # Checks
        if config.contains is not None and config.contains not in name: continue
        if config.not_contains is not None and config.not_contains in name: continue
        if config.exact_name is not None and name != config.exact_name: continue
        if config.exact_not_name is not None and name == config.exact_not_name: continue
        if config.startswith is not None and not name.startswith(config.startswith): continue
        if config.endswith is not None and not name.endswith(config.endswith): continue

        # Get the description
        description = properties[name]

        # Get the values for all the simulations
        values = get_values_for_simulations(simulations, name)

        # Get unique values
        unique_values = sequences.unique_values(values.values())

        if len(unique_values) == 1:
            default = unique_values[0]
            ptype, string = stringify(default)
            choices = None
            suggestions = None
        else:
            default = None
            ptype = get_common_ptype(unique_values)
            choices = None
            suggestions = unique_values

        # Ask for the new value
        value = prompt_variable(name, ptype, description, default=default, required=True, choices=choices, suggestions=suggestions)
        #if default is None and value == "": continue

        # Each simulation had the same value: adapt each simulation simultaneously without prompting to proceed
        if len(unique_values) == 1:
            if value != default: set_value_for_simulations(simulations, name, value)

        # Different simulations had different values
        else: set_value_for_simulations_prompt(simulations, name, value)

# -----------------------------------------------------------------
