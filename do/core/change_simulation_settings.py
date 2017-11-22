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
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_variable
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.basics.log import log
from pts.core.tools.stringify import stringify, tostr

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_required("id", "positive_integer", "simulation ID")
definition.add_positional_optional("matching", "string", "only adapt settings with a name matching this string")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("change_simulation_settings", definition, description="Change certain settings of a simulation")

# -----------------------------------------------------------------

# Open the simulation object
simulation = get_simulation_for_host(config.remote, config.id)

# -----------------------------------------------------------------

contains = config.matching
not_contains = None
exact_name = None
exact_not_name = None

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
#properties["handle"] = "execution handle"
properties["retrieved"] = "retrieved flag"

# -----------------------------------------------------------------

# Loop over the properties
for name in properties:

    # Checks
    if contains is not None and contains not in name: continue
    if not_contains is not None and not_contains in name: continue
    if exact_name is not None and name != exact_name: continue
    if exact_not_name is not None and name == exact_not_name: continue

    # Get the description
    description = properties[name]

    # Get the current value
    default = getattr(simulation, name)

    # There is a current value
    if default is not None: ptype, pstring = stringify(default)

    # No current value
    else: ptype = "any"

    # Check ptype
    if ptype is None: ptype = "any"

    # Ask for the new value
    value = prompt_variable(name, ptype, description, default=default, required=True)

    # Set the property
    if value != default:

        # Debugging
        log.debug("Changing the value of '" + name + "' to '" + tostr(value) + "' ...")

        # Set the new value
        setattr(simulation, name, value)

# -----------------------------------------------------------------

# Save the simulation
simulation.save()

# -----------------------------------------------------------------
