#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.adapter Contains the SimulationAdapter and AnalysisAdapter classes.

# -----------------------------------------------------------------

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..launch.options import get_analysis_property_names_and_descriptions, get_analysis_section_names_and_descriptions, get_analysis_property_names_and_descriptions_for_section
from ..basics.configuration import prompt_proceed, prompt_variable
from ..tools.utils import lazyproperty
from ..simulation.remote import get_simulation_for_host
from ..basics.log import log
from ..tools.stringify import stringify, tostr
from ..tools import sequences
from ..basics.configuration import parent_type
from ..tools import filesystem as fs
from ..simulation.remote import get_simulations_for_host
from ..tools import introspection

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

class SimulationAdapter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SimulationAdapter, self).__init__(*args, **kwargs)

        # The simulations
        self.simulations = OrderedDict()

        # Flags
        self.changed = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Adapt the settings
        self.adapt()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationAdapter, self).setup(**kwargs)

        # Check settings
        if self.config.matching is not None:
            if self.config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
            self.config.contains = self.config.matching

        # Set simulation names
        if self.config.from_directories:
            if self.config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
            self.config.names = fs.directories_in_path(returns="name")

        # Load the simulations
        self.load_simulations()

    # -----------------------------------------------------------------

    @property
    def simulation_ids(self):

        """
        This function ...
        :return:
        """

        return self.simulations.keys()

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.simulations)

    # -----------------------------------------------------------------

    @property
    def has_single_simulation(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations == 1

    # -----------------------------------------------------------------

    @property
    def single_simulation(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        return self.simulations[self.simulation_ids[0]]

    # -----------------------------------------------------------------

    @property
    def single_changed(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        return self.changed[self.simulation_ids[0]]

    # -----------------------------------------------------------------

    @single_changed.setter
    def single_changed(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        self.changed[self.simulation_ids[0]] = value

    # -----------------------------------------------------------------

    @lazyproperty
    def all_simulations(self):

        """
        This function ...
        :return:
        """

        return get_simulations_for_host(self.config.remote, as_dict=True)

    # -----------------------------------------------------------------

    @property
    def all_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.all_simulations.keys()

    # -----------------------------------------------------------------

    def load_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulations ...")

        # Set simulation IDs
        if self.config.names is None and self.config.ids is None: self.config.ids = introspection.simulation_ids_for_host(self.config.remote)

        # From names
        if self.config.names is not None:

            # Loop over the names
            for name in self.config.names:
                simulation_id = self.all_simulations[name].id
                self.simulations[simulation_id] = self.all_simulations[name]

        # From IDS
        elif self.config.ids is not None:

            # Load the simulations and put them in the dictionary
            for simulation_id in self.config.ids: self.simulations[simulation_id] = get_simulation_for_host(self.config.remote, simulation_id)

        # Nothing
        else: raise ValueError("Names or IDs must be specified")

        # Set changed flags
        for simulation_id in self.simulation_ids: self.changed[simulation_id] = False

    # -----------------------------------------------------------------

    def adapt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting ...")

        # Single simulation?
        if self.has_single_simulation: self.adapt_single()

        # Multiple simulations
        else: self.adapt_multi()

    # -----------------------------------------------------------------

    def adapt_single(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting single simulation ...")

        # Loop over the properties
        for name in self.property_names:

            # Get the description
            description = properties[name]

            # Get the current value
            default = get_value_for_simulation(self.single_simulation, name)

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
            if value == default: continue

            # Change
            set_value_for_simulation(self.single_simulation, name, value)
            self.single_changed = True

    # -----------------------------------------------------------------

    @lazyproperty
    def property_names(self):

        """
        This function ...
        :return:
        """

        names = []

        # Loop over the properties
        for name in properties:

            # Checks
            if self.config.contains is not None and self.config.contains not in name: continue
            if self.config.not_contains is not None and self.config.not_contains in name: continue
            if self.config.exact_name is not None and name != self.config.exact_name: continue
            if self.config.exact_not_name is not None and name == self.config.exact_not_name: continue
            if self.config.startswith is not None and not name.startswith(self.config.startswith): continue
            if self.config.endswith is not None and not name.endswith(self.config.endswith): continue

            # Add the name
            names.append(name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    def adapt_multi(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting multiple simulations ...")

        # Loop over the properties
        for name in self.property_names:

            # Get the description
            description = properties[name]

            # Get the values for all the simulations
            values = get_values_for_simulations(self.simulations, name)

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
                change = prompt_proceed("Change the property '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
                if not change: continue
                default = None
                ptype = get_common_ptype(unique_values)
                choices = None
                suggestions = unique_values

            # Ask for the new value
            value = prompt_variable(name, ptype, description, default=default, required=True, choices=choices, suggestions=suggestions)
            # if default is None and value == "": continue

            # Each simulation had the same value: adapt each simulation simultaneously without prompting to proceed
            if len(unique_values) == 1:
                if value != default:
                    set_value_for_simulations(self.simulations, name, value)
                    for simulation_id in self.simulation_ids: self.changed[simulation_id] = True

            # Different simulations had different values
            else: self.changed = set_value_for_simulations_prompt(self.simulations, name, value, changed=self.changed)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the simulations
        self.write_simulations()

    # -----------------------------------------------------------------

    def write_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulations ...")

        # Single
        if self.has_single_simulation: save_simulation(self.single_simulation)

        # Multi
        else: save_simulations(self.simulations, changed=self.changed)

# -----------------------------------------------------------------

class AnalysisAdapter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AnalysisAdapter, self).__init__(*args, **kwargs)

        # The simulations
        self.simulations = OrderedDict()

        # Flags
        self.changed = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Adapt the settings
        self.adapt()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisAdapter, self).setup(**kwargs)

        # Check settings
        if self.config.matching is not None:
            if self.config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
            self.config.contains = self.config.matching

        # Set simulation names
        if self.config.from_directories:
            if self.config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
            self.config.names = fs.directories_in_path(returns="name")

        # Load the simulations
        self.load_simulations()

    # -----------------------------------------------------------------

    @property
    def simulation_ids(self):

        """
        This function ...
        :return:
        """

        return self.simulations.keys()

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        This function ...
        """

        return len(self.simulations)

    # -----------------------------------------------------------------

    @property
    def has_single_simulation(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations == 1

    # -----------------------------------------------------------------

    @property
    def single_simulation(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        return self.simulations[self.simulation_ids[0]]

    # -----------------------------------------------------------------

    @property
    def single_changed(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        return self.changed[self.simulation_ids[0]]

    # -----------------------------------------------------------------

    @single_changed.setter
    def single_changed(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        self.changed[self.simulation_ids[0]] = value

    # -----------------------------------------------------------------

    @lazyproperty
    def all_simulations(self):

        """
        This function ...
        :return:
        """

        return get_simulations_for_host(self.config.remote, as_dict=True)

    # -----------------------------------------------------------------

    @property
    def all_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.all_simulations.keys()

    # -----------------------------------------------------------------

    def load_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulations ...")

        # Set simulation IDs
        if self.config.names is None and self.config.ids is None: self.config.ids = introspection.simulation_ids_for_host(self.config.remote)

        # From names
        if self.config.names is not None:

            # Loop over the names
            for name in self.config.names:
                simulation_id = self.all_simulations[name].id
                self.simulations[simulation_id] = self.all_simulations[name]

        # From IDS
        elif self.config.ids is not None:

            # Load the simulations and put them in the dictionary
            for simulation_id in self.config.ids: self.simulations[simulation_id] = get_simulation_for_host(self.config.remote, simulation_id)

        # Nothing
        else: raise ValueError("Names or IDs must be specified")

        # Set changed flags
        for simulation_id in self.simulation_ids: self.changed[simulation_id] = False

    # -----------------------------------------------------------------

    @lazyproperty
    def properties(self):

        """
        This function ...
        :return:
        """

        return get_analysis_property_names_and_descriptions()

    # -----------------------------------------------------------------

    @lazyproperty
    def property_names(self):

        """
        This function ...
        :return:
        """

        # Initialize list for the property names
        names = []

        # Loop over the properties
        for name in self.properties:

            # Checks
            if self.config.contains is not None and self.config.contains not in name: continue
            if self.config.not_contains is not None and self.config.not_contains in name: continue
            if self.config.exact_name is not None and name != self.config.exact_name: continue
            if self.config.exact_not_name is not None and name == self.config.exact_not_name: continue
            if self.config.startswith is not None and not name.startswith(self.config.startswith): continue
            if self.config.endswith is not None and not name.endswith(self.config.endswith): continue

            # Add the name
            names.append(name)

        # Return
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def sections(self):

        """
        This function ...
        :return:
        """

        # Get sections
        return get_analysis_section_names_and_descriptions()

    # -----------------------------------------------------------------

    @lazyproperty
    def section_names(self):

        """
        This function ...
        :return:
        """

        return self.sections.keys()

    # -----------------------------------------------------------------

    def adapt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting ...")

        # Single simulation
        if self.has_single_simulation: self.adapt_single()

        # Multiple simulations
        else: self.adapt_multi()

    # -----------------------------------------------------------------

    def adapt_single(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting single simulation ...")

        # Check whether analysis options are defined
        has_changed = self.single_simulation.analysis.prompt_properties(contains=self.config.contains, not_contains=self.config.not_contains,
                                              exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                              startswith=self.config.startswith, endswith=self.config.endswith)

        # Set changed flag
        self.single_changed = has_changed

    # -----------------------------------------------------------------

    def adapt_multi(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting multiple simulations ...")

        # Properties
        self.adapt_properties()

        # Sections
        self.adapt_sections()

    # -----------------------------------------------------------------

    def adapt_properties(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting properties ...")

        # Loop over the properties
        for name in self.property_names:

            # Get description
            description = self.properties[name]

            # Get the analysis options for all the simulations
            values = get_analysis_values_for_simulations(self.simulations, name)

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
                    set_analysis_value_for_simulations(self.simulations, name, value)
                    for simulation_id in self.simulations: self.changed[simulation_id] = True

            # Different simulations had different values
            else: self.changed = set_analysis_value_for_simulations_prompt(self.simulations, name, value, changed=self.changed)

    # -----------------------------------------------------------------

    def adapt_sections(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adapting sections ...")

        # Loop over the sections
        for section_name in self.section_names:

            # Get section description
            section_description = self.sections[section_name]

            # Debug
            log.debug("Entering section '" + section_name + "' ...")

            # Get properties
            properties = get_analysis_property_names_and_descriptions_for_section(section_name)

            # Loop over the properties in this section
            for name in properties:

                # Checks
                if self.config.contains is not None and self.config.contains not in name: continue
                if self.config.not_contains is not None and self.config.not_contains in name: continue
                if self.config.exact_name is not None and name != self.config.exact_name: continue
                if self.config.exact_not_name is not None and name == self.config.exact_not_name: continue
                if self.config.startswith is not None and not name.startswith(self.config.startswith): continue
                if self.config.endswith is not None and not name.endswith(self.config.endswith): continue

                # Get description
                description = properties[name]

                # Get the analysis options for all the simulations
                values = get_analysis_values_for_simulations(self.simulations, name, section=section_name)

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
                        set_analysis_value_for_simulations(self.simulations, name, value, section=section_name)
                        for simulation_id in self.simulations: self.changed[simulation_id] = True

                    # Different simulations had different values
                    else: self.changed = set_analysis_value_for_simulations_prompt(self.simulations, name, value, changed=self.changed, section=section_name)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the simulations
        self.write_simulations()

    # -----------------------------------------------------------------

    def write_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulation(s) ...")

        # Save the simulations
        if self.has_single_simulation: save_simulation(self.single_simulation)

        # Multiple simulations
        else: save_simulations(self.simulations, changed=self.changed)

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

def set_value_for_simulations_prompt(simulations, name, value, changed=None):

    """
    This function ...
    :param simulations:
    :param name:
    :param value:
    :param changed:
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
        current_value = get_value_for_simulation(simulation, name)

        # Check
        if current_value == value: continue

        # Prompt
        if prompt_proceed("Replace the value of '" + name + "' from '" + tostr(current_value) + "' to '" + tostr(value) + " for simulation '" + simulation_name + "'?"):
            set_value_for_simulation(simulation, name, value)
            changed[simulation_id] = True

    # Return the changed dictionary
    return changed

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

def save_simulation(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    # Debugging
    log.debug("Saving the simulation ...")

    # Save the simulation
    simulation.save()

# -----------------------------------------------------------------

def save_simulations(simulations, changed=None):

    """
    This function ...
    :param simulations:
    :param changed:
    :return:
    """

    # Debugging
    log.debug("Saving the simulations ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Not changed?
        if changed is not None and not changed[simulation_id]: continue

        # Save
        simulations[simulation_id].save()

# -----------------------------------------------------------------
