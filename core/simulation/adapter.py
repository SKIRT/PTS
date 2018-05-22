#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.adapter Contains the SimulationAdapter and AnalysisAdapter classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

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
from ..tools import filesystem as fs
from ..simulation.remote import get_simulations_for_host
from ..tools import introspection
from .shower import properties, get_common_ptype, get_values_for_simulations, get_value_for_simulation, get_analysis_values_for_simulations, get_analysis_value_for_simulation
from ..tools import types
from ..basics import containers

# -----------------------------------------------------------------

def adapt_simulation(simulation, config=None):

    """
    This function ...
    :param simulation:
    :param config:
    :return:
    """

    # Create adapter
    adapter = SimulationAdapter(config=config)

    # Run
    adapter.run(simulation=simulation)

# -----------------------------------------------------------------

def adapt_simulations(*simulations, **kwargs):

    """
    This function ...
    :param simulations:
    :param kwargs:
    :return:
    """

    # Create adapter
    adapter = SimulationAdapter(kwargs)

    # Run
    adapter.run(simulations=simulations)

# -----------------------------------------------------------------

def adapt_analysis(simulation, config=None):

    """
    This function ...
    :param simulation:
    :param config:
    :return:
    """

    # Create adapter
    adapter = AnalysisAdapter(config=config)

    # Run
    adapter.run(simulation=simulation)

# -----------------------------------------------------------------

def adapt_analysis_simulations(*simulations, **kwargs):

    """
    This function ...
    :param simulations:
    :param kwargs:
    :return:
    """

    # Create adapter
    adapter = AnalysisAdapter(kwargs)

    # Run
    adapter.run(simulations=simulations)

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Get simulations
        if "simulations" in kwargs:
            simulations = kwargs.pop("simulations")
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute="id")
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")
        elif "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.simulations[simulation.id] = simulation

        # Load the simulations
        if not self.has_simulations: self.load_simulations()

    # -----------------------------------------------------------------

    @property
    def has_simulations(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations > 0

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
                if name not in self.all_simulations: raise ValueError("Simulation '" + name + "' is not a simulation of host '" + self.config.remote + "'")
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

            # Check types
            if self.config.types is not None and ptype not in self.config.types: continue

            # Replace?
            if ptype == "string" and self.config.replace_string is not None and default is not None and self.config.replace_string[0] in default:

                #print(default)
                value = default.replace(self.config.replace_string[0], self.config.replace_string[1])
                #print(value)

            # No replacements: skip
            elif self.config.only_replacements: continue

            # No replacements: prompt for new value
            else:

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

            # Check types
            if self.config.types is not None and ptype not in self.config.types: continue

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Get simulations
        if "simulations" in kwargs:
            simulations = kwargs.pop("simulations")
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute="id")
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")
        elif "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.simulations[simulation.id] = simulation

        # Load the simulations
        if not self.has_simulations: self.load_simulations()

    # -----------------------------------------------------------------

    @property
    def has_simulations(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations > 0

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

        # Update analysis options in each simulation
        if self.config.update:
            for simulation_id in self.simulation_ids:
                self.simulations[simulation_id].update_analysis_options()
                #print(self.simulations[simulation_id].analysis.misc)

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
                                              startswith=self.config.startswith, endswith=self.config.endswith, replace_string=self.config.replace_string,
                                              types=self.config.types, only_replacements=self.config.only_replacements)

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

    def get_ptype(self, name, section_name=None):

        """
        This function ...
        :param name:
        :param section_name:
        :return:
        """

        if section_name is not None: return self.simulations[self.simulation_ids[0]].analysis[section_name].get_ptype(name)
        else: return self.simulations[self.simulation_ids[0]].analysis.get_ptype(name)

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

            # Get the ptype
            ptype = self.get_ptype(name)

            # Check types
            if self.config.types is not None and ptype not in self.config.types: continue

            # Only one unique value
            if len(unique_values) == 1:

                default = unique_values[0]
                #ptype, string = stringify(default)
                choices = None
                suggestions = None

            # Multiple unique values
            else:

                # Prompt to change this property
                change = prompt_proceed("Change the analysis option '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
                if not change: continue
                default = None
                #ptype = get_common_ptype(unique_values)
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
                #print(unique_values)

                # Get the ptype
                ptype = self.get_ptype(name, section_name)

                # Check types
                if self.config.types is not None and ptype not in self.config.types: continue

                # Only one unique value
                if len(unique_values) == 1:

                    default = unique_values[0]
                    #ptype, string = stringify(default)
                    choices = None
                    suggestions = None

                # Multiple unique values
                else:

                    # Prompt to change this property
                    change = prompt_proceed("Change the analysis option '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
                    if not change: continue
                    default = None
                    #ptype = get_common_ptype(unique_values)
                    choices = None
                    suggestions = unique_values

                # Ask for the new value
                #print(name, ptype, description, default, choices, suggestions)
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
    if section is not None: log.debug("in section '" + section + "'")
    if section is not None: log.debug("Original value: '" + tostr(simulation.analysis[section][name]) + "'")
    else: log.debug("Original value: '" + tostr(simulation.analysis[name]) + "'")

    #print(simulation.analysis[section][name])
    #print(value)

    # Set the new value
    if section is not None:
        section = simulation.analysis[section]
        #print(section)
        section[name] = value
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
    if section is not None: log.debug("in section '" + section + "'")

    # Loop over the simulations
    for simulation_id in simulations:

        # Debugging
        #log.debug("Original value: '" + get_analysis_value_for_simulation(simulations[simulation_id], name, section=section) + "'")

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
    if section is not None: log.debug("in section '" + section + "'")

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
        if prompt_proceed("Replace the value of '" + name + "' from '" + tostr(current_value) + "' to '" + tostr(value) + "' for simulation '" + simulation_name + "'?"):
            set_analysis_value_for_simulation(simulation, name, value, section=section)
            changed[simulation_id] = True

    # Return the changed dictionary
    return changed

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
    log.debug("Original value: '" + tostr(getattr(simulation, name)) + "' ...")

    # Set the new value
    setattr(simulation, name, value)

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
        if prompt_proceed("Replace the value of '" + name + "' from '" + tostr(current_value) + "' to '" + tostr(value) + "' for simulation '" + simulation_name + "'?"):
            set_value_for_simulation(simulation, name, value)
            changed[simulation_id] = True

    # Return the changed dictionary
    return changed

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
