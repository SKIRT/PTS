#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.shower Contains the SimulationShower and AnalysisShower classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..launch.options import get_analysis_property_names_and_descriptions, get_analysis_section_names_and_descriptions, get_analysis_property_names_and_descriptions_for_section
from ..tools.utils import lazyproperty
from ..simulation.remote import get_simulation_for_host
from ..basics.log import log
from ..tools.stringify import stringify, stringify_list_fancy
from ..tools import sequences
from ..basics.configuration import parent_type
from ..tools import filesystem as fs
from ..simulation.remote import get_simulations_for_host
from ..tools import introspection
from ..tools import types
from ..basics import containers
from ..tools import formatting as fmt
from ..basics.configuration import prompt_choices, prompt_yn

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

def select_simulation_settings(successive=False, contains=None, not_contains=None, exact_name=None, exact_not_name=None,
                               startswith=None, endswith=None):

    """
    This function ...
    :param successive:
    :param contains:
    :param not_contains:
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :return:
    """

    # Initialize list for the selected settings
    settings = []

    # Select options in succession
    if successive:

        # Loop over the properties
        for name in properties:

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue
            if startswith is not None and not name.startswith(startswith): continue
            if endswith is not None and not name.endswith(endswith): continue

            # Select?
            description = properties[name]
            if not prompt_yn(name, "select '" + description + "'?", default=False): continue

            # Add the setting
            settings.append(name)

    # Show all options at the same time
    else: settings = prompt_choices("settings", "simulation settings", properties)

    # Return the choosen settings
    return settings

# -----------------------------------------------------------------

def select_analysis_options(successive=False, hierarchic=True, contains=None, not_contains=None, exact_name=None,
                            exact_not_name=None, startswith=None, endswith=None):

    """
    This function ...
    :param successive:
    :param hierarchic:
    :param contains:
    :param not_contains:
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :return:
    """

    # Initialize list for the selected options
    options = []

    # Hierachic: first properties, then sections
    if hierarchic:

        ## Properties

        # Get properties
        all_properties = get_analysis_property_names_and_descriptions()

        # Subset of properties
        properties = OrderedDict()

        # Loop over the properties
        for name in all_properties:

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue
            if startswith is not None and not name.startswith(startswith): continue
            if endswith is not None and not name.endswith(endswith): continue

            # Add
            properties[name] = all_properties[name]

        # Show
        print("")
        print(fmt.yellow + fmt.bold + "PROPERTIES" + fmt.reset)
        print("")

        # Select properties in succession
        if successive:

            # Loop over the properties
            for name in properties:

                # Select?
                description = properties[name]
                if not prompt_yn(name, "select '" + description + "'?", default=False): continue

                # Add the property
                options.append(name)

        # Show all properties at the same time
        else: options += prompt_choices("properties", "analysis properties", properties)

        ## Sections

        # Get sections
        sections = get_analysis_section_names_and_descriptions()

        # Loop over the sections
        for section_name in sections:

            # Get section description
            section_description = sections[section_name]

            # Show section
            print("")
            print(fmt.yellow + fmt.bold + section_name.upper() + fmt.reset_bold + ": " + section_description)
            print("")

            # Debug
            log.debug("Entering section '" + section_name + "' ...")

            # Get properties
            all_section_properties = get_analysis_property_names_and_descriptions_for_section(section_name)

            # Make selection
            section_properties = OrderedDict()

            # Loop over the properties
            for name in all_section_properties:

                # Checks
                if contains is not None and contains not in name: continue
                if not_contains is not None and not_contains in name: continue
                if exact_name is not None and name != exact_name: continue
                if exact_not_name is not None and name == exact_not_name: continue
                if startswith is not None and not name.startswith(startswith): continue
                if endswith is not None and not name.endswith(endswith): continue

                # Add
                section_properties[name] = all_section_properties[name]

            # Any?
            if len(section_properties) == 0: continue

            # Select properties in succession
            if successive:

                # Loop over the properties
                for name in section_properties:

                    # Select
                    description = section_properties[name]
                    if not prompt_yn(name, "select '" + description + "' in '" + section_name + "'?", default=False): continue

                    # Add the property
                    options.append((section_name, name))

            # Show all properties at the same time
            else:

                # Select
                section_options = prompt_choices("properties", section_name + " properties", section_properties)

                # Add
                for name in section_options: options.append((section_name, name))

    # Not hierarchic
    else:

        # Get properties
        all_properties = get_analysis_property_names_and_descriptions()

        # Make subset
        properties = OrderedDict()

        # Loop over the properties
        for name in all_properties:

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue
            if startswith is not None and not name.startswith(startswith): continue
            if endswith is not None and not name.endswith(endswith): continue

            # Add property
            properties[name] = all_properties[name]

        # Get sections
        sections = get_analysis_section_names_and_descriptions()

        # Add section properties
        for section_name in sections:

            # Get properties
            all_section_properties = get_analysis_property_names_and_descriptions_for_section(section_name)

            # Make subset
            section_properties = OrderedDict()

            # Loop over the properties
            for name in all_section_properties:

                # Checks
                if contains is not None and contains not in name: continue
                if not_contains is not None and not_contains in name: continue
                if exact_name is not None and name != exact_name: continue
                if exact_not_name is not None and name == exact_not_name: continue
                if startswith is not None and not name.startswith(startswith): continue
                if endswith is not None and not name.endswith(endswith): continue

                # Add the property
                section_properties[name] = all_section_properties[name]

            # Loop over the properties
            for name in section_properties:

                # Create full property name
                full_name = section_name + "/" + name

                # Get description
                description = section_properties[name]

                # Add
                properties[full_name] = description

        # Select properties in succession
        if successive:

            # Loop over the properties
            for full_name in properties:

                # Select?
                description = properties[full_name]
                if not prompt_yn(full_name, "select '" + description + "'?", default=False): continue

                # Add the property
                if "/" in full_name:
                    section_name, property_name = full_name.split("/")
                    options.append((section_name, property_name))
                else: options.append(full_name)

        # Show all properties at the same time
        else:

            # Get the choices
            choices = prompt_choices("properties", "analysis properties", properties)

            # Add the choices
            for full_name in choices:

                # Add the property
                if "/" in full_name:
                    section_name, property_name = full_name.split("/")
                    options.append((section_name, property_name))
                else: options.append(full_name)

    # Return the choosen options
    return options

# -----------------------------------------------------------------

def compare_simulations(*simulations, **kwargs):

    """
    This function ...
    :param simulations:
    :param kwargs:
    :return:
    """

    # Create shower
    shower = SimulationShower(config=kwargs.pop("config", None))

    # Run
    shower.run(simulations=simulations)

# -----------------------------------------------------------------

def compare_analysis(*simulations, **kwargs):

    """
    This function ...
    :param simulations:
    :param kwargs:
    :return:
    """

    # Create shower
    shower = AnalysisShower(config=kwargs.pop("config", None))

    # Run
    shower.run(simulations=simulations)

# -----------------------------------------------------------------

def show_simulation(simulation, config=None):

    """
    This function ...
    :param simulation:
    :param config:
    :return:
    """

    # Create shower
    shower = SimulationShower(config=config)

    # Show
    shower.run(simulation=simulation)

# -----------------------------------------------------------------

def show_analysis(simulation, config=None):

    """
    This function ...
    :param simulation:
    :param config:
    :return:
    """

    # Create shower
    shower = AnalysisShower(config=config)

    # Show
    shower.run(simulation=simulation)

# -----------------------------------------------------------------

class SimulationShower(Configurable):

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
        super(SimulationShower, self).__init__(*args, **kwargs)

        # The simulations
        self.simulations = OrderedDict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Show
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationShower, self).setup(**kwargs)

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
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute=["host_id", "id"])
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")
        elif "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.simulations[(simulation.host_id, simulation.id)] = simulation

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
    def simulation_keys(self):

        """
        This function ...
        :return:
        """

        return self.simulations.keys()

    # -----------------------------------------------------------------

    @property
    def simulation_host_ids(self):

        """
        This function ...
        :return:
        """

        return [key[0] for key in self.simulations.keys()]

    # -----------------------------------------------------------------

    @property
    def simulation_ids(self):

        """
        This function ...
        :return:
        """

        return [key[1] for key in self.simulations.keys()]

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return [simulation.name for simulation in self.simulations.values()]

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
        return self.simulations[self.simulation_keys[0]]

    # -----------------------------------------------------------------

    @property
    def single_simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.single_simulation.name

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
                host_id = self.all_simulations[name].host_id
                simulation_id = self.all_simulations[name].id
                self.simulations[(host_id, simulation_id)] = self.all_simulations[name]

        # From IDS
        elif self.config.ids is not None:

            # Load the simulations and put them in the dictionary
            for simulation_id in self.config.ids: self.simulations[(self.config.remote, simulation_id)] = get_simulation_for_host(self.config.remote, simulation_id)

        # Nothing
        else: raise ValueError("Names or IDs must be specified")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Single simulation?
        if self.has_single_simulation: self.show_single()

        # Multiple simulations
        else: self.show_multi()

    # -----------------------------------------------------------------

    def show_single(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing single simulation ...")

        # Show the simulation name
        print("")
        print(fmt.blue + fmt.underlined + self.single_simulation_name + fmt.reset)
        print("")

        # Loop over the properties
        for name in self.property_names:

            # Get the description
            description = properties[name]

            # Get the current value
            default = get_value_for_simulation(self.single_simulation, name)

            # There is a current value
            if default is not None:

                ptype, pstring = stringify(default)

                if pstring == "":
                    type_string = "[empty " + ptype + "]" if ptype is not None else ""
                    print(" - " + fmt.bold + name + fmt.reset + " (" + description + "): " + type_string)
                else:
                    type_string = " [" + ptype + "]" if ptype is not None else ""
                    print(" - " + fmt.bold + name + fmt.reset + " (" + description + "): " + pstring + type_string)

            else: print(" - " + fmt.bold + name + fmt.reset + "(" + description +  "): None")

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

    def show_multi(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing multiple simulations ...")

        # Show simulation names
        print("")
        print("Simulations:")
        print("")
        for index, name in enumerate(self.simulation_names): print("(" + str(index) + ") " + fmt.blue + fmt.underlined + name + fmt.reset)
        print("")

        # Loop over the properties
        for name in self.property_names:

            # Get the description
            description = properties[name]

            # Get the values for all the simulations
            values_for_names = get_values_for_simulations(self.simulations, name, keys="name")
            names_for_values = containers.invert_dictionary(values_for_names)

            # Get unique values
            unique_values = sequences.unique_values(values_for_names.values())
            nunique_values = len(unique_values)

            # One unique value
            if nunique_values == 1:

                # Stringify
                ptype, string = stringify(unique_values[0])

                # Empty string, list, tuple, ...
                if string == "":

                    type_string = "[empty " + ptype + "]" if ptype is not None else "[empty]"
                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + type_string + fmt.reset)

                # Normal value
                else:

                    type_string = " [" + ptype + "]" if ptype is not None else ""
                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + type_string + fmt.reset)

            # Different unique values
            else:

                # Show name of the property
                print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + fmt.reset)

                # Loop over the unique values
                for value in unique_values:

                    # Get the corresponding simulations
                    simulation_names = names_for_values[value]
                    nsimulations = len(simulation_names)

                    # Stringify
                    ptype, string = stringify(value)

                    # Empty string, list, tuple, ...
                    if string == "":

                        type_string = "[empty " + ptype + "]" if ptype is not None else "[empty]"
                        print("   * " + type_string + fmt.bold + " (" + str(nsimulations) + ")" + fmt.reset_bold)

                    # Normal value
                    else:

                        type_string = " [" + ptype + "]" if ptype is not None else ""
                        print("   * " + string + type_string + fmt.bold + " (" + str(nsimulations) + ")" + fmt.reset_bold)

                    # Show simulation names
                    _, simulation_names_string = stringify_list_fancy(simulation_names, lines_prefix="     ")
                    print(fmt.yellow + simulation_names_string + fmt.reset)

# -----------------------------------------------------------------

class AnalysisShower(Configurable):

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
        super(AnalysisShower, self).__init__(*args, **kwargs)

        # The simulations
        self.simulations = OrderedDict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Show the options
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisShower, self).setup(**kwargs)

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
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute=["host_id", "id"])
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")
        elif "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.simulations[(simulation.host_id, simulation.id)] = simulation

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
    def simulation_keys(self):

        """
        This function ...
        :return:
        """

        return self.simulations.keys()

    # -----------------------------------------------------------------

    @property
    def simulation_host_ids(self):

        """
        This function ...
        :return:
        """

        return [key[0] for key in self.simulations.keys()]

    # -----------------------------------------------------------------

    @property
    def simulation_ids(self):

        """
        This function ...
        :return:
        """

        return [key[1] for key in self.simulations.keys()]

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return [simulation.name for simulation in self.simulations.values()]

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
        return self.simulations[self.simulation_keys[0]]

    # -----------------------------------------------------------------

    @property
    def single_simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.single_simulation.name

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
                host_id = self.all_simulations[name].host_id
                simulation_id = self.all_simulations[name].id
                self.simulations[(host_id, simulation_id)] = self.all_simulations[name]

        # From IDS
        elif self.config.ids is not None:

            # Load the simulations and put them in the dictionary
            for simulation_id in self.config.ids: self.simulations[(self.config.remote, simulation_id)] = get_simulation_for_host(self.config.remote, simulation_id)

        # Nothing
        else: raise ValueError("Names or IDs must be specified")

        # Update analysis options in each simulation
        if self.config.update:
            for key in self.simulation_keys:
                self.simulations[key].update_analysis_options()

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

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Single simulation
        if self.has_single_simulation: self.show_single()

        # Multiple simulations
        else: self.show_multi()

    # -----------------------------------------------------------------

    def show_single(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing single simulation ...")

        # Show simulation name
        print("")
        print(fmt.blue + fmt.underlined + self.single_simulation_name + fmt.reset)
        print("")

        # Show
        self.single_simulation.analysis.show_properties(contains=self.config.contains, not_contains=self.config.not_contains,
                                              exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                              startswith=self.config.startswith, endswith=self.config.endswith)

    # -----------------------------------------------------------------

    def show_multi(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing multiple simulations ...")

        # Show simulation names
        print("")
        print("Simulations:")
        print("")
        for index, name in enumerate(self.simulation_names): print("(" + str(index) + ") " + fmt.underlined + fmt.blue + name + fmt.reset)
        print("")

        # Properties
        self.show_properties()

        # Sections
        self.show_sections()

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

    def show_properties(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing properties ...")

        # Loop over the properties
        for name in self.property_names:

            # Get description
            description = self.properties[name]

            # Get the analysis options for all the simulations
            #values = get_analysis_values_for_simulations(self.simulations, name)
            values_for_names = get_analysis_values_for_simulations(self.simulations, name, keys="name")
            names_for_values = containers.invert_dictionary(values_for_names)

            # Get unique values
            unique_values = sequences.unique_values(values_for_names.values())
            nunique_values = len(unique_values)

            # Get the ptype
            #ptype = self.get_ptype(name)

            # Only one unique value
            if nunique_values == 1:

                # Stringify the value
                default = unique_values[0]
                ptype, string = stringify(default)

                # Empty string, list, tuple, ...
                if string == "":

                    type_string = "[empty " + ptype + "]" if ptype is not None else "[empty]"
                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + type_string + fmt.reset)

                # Normal value
                else:

                    type_string = " [" + ptype + "]" if ptype is not None else ""
                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + type_string + fmt.reset)

            # Multiple unique values
            else:

                # Show the name of the property
                print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "):" + fmt.reset)

                # Show the different unique values
                for value in unique_values:

                    # Get the simulations for this value
                    simulation_names = names_for_values[value]
                    nsimulations = len(simulation_names)

                    # Stringify the property
                    ptype, string = stringify(value)

                    # Empty string, list, tuple, ...
                    if string == "":

                        type_string = "[empty " + ptype + "]" if ptype is not None else "[empty]"
                        print("   * " + type_string + fmt.bold + " (" + str(nsimulations) + ")" + fmt.reset_bold)

                    # Normal value
                    else:

                        type_string = " [" + ptype + "]" if ptype is not None else ""
                        print("   * " + string + type_string + fmt.bold + " (" + str(nsimulations) + ")" + fmt.reset_bold)

                    # Show simulation names
                    _, simulation_names_string = stringify_list_fancy(simulation_names, lines_prefix="     ")
                    print(fmt.yellow + simulation_names_string + fmt.reset)

    # -----------------------------------------------------------------

    def show_sections(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Showing sections ...")

        # Loop over the sections
        for section_name in self.section_names:

            # Get section description
            section_description = self.sections[section_name]

            # Show section
            print(fmt.bold + section_name + fmt.reset + ": " + section_description)

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
                nunique_values = len(unique_values)

                # Get the ptype
                #ptype = self.get_ptype(name, section_name)
                #print(name, ptype, unique_values)

                # Only one unique value
                if nunique_values == 1:

                    default = unique_values[0]
                    ptype, string = stringify(default)
                    type_string = " [" + ptype + "]" if ptype is not None else ""
                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + type_string + fmt.reset)

                # Multiple unique values
                else:

                    print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "):" + fmt.reset)

                    for value in unique_values:

                        ptype, string = stringify(value)
                        type_string = " [" + ptype + "]" if ptype is not None else ""
                        print("    * " + string + type_string)

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

def get_analysis_values_for_simulations(simulations, name, section=None, keys="id"):

    """
    This function ...
    :param simulations:
    :param name:
    :param section:
    :param keys:
    :return:
    """

    values = dict()

    # Loop over the simulations
    for key in simulations:

        # Get simulation
        simulation = simulations[key]

        # Get value
        value = get_analysis_value_for_simulation(simulation, name, section=section)

        # Set key
        if keys == "id": key = simulation.id
        elif keys == "name": key = simulation.name
        else: raise ValueError("Invalid value for 'key'")

        # Set value
        values[key] = value

    # Return the values
    return values

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

def get_values_for_simulations(simulations, name, keys="id"):

    """
    This function ...
    :param simulations:
    :param name:
    :param keys:
    :return:
    """

    values = dict()

    # Loop over the simulations
    for key in simulations:

        # Get the simulation
        simulation = simulations[key]

        # Get value
        value = get_value_for_simulation(simulation, name)

        # Set the key
        if keys == "id": key = simulation.id
        elif keys == "name": key = simulation.name
        else: raise ValueError("Invalid value for 'key'")

        # Set value
        values[key] = value

    # Return the values
    return values

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
