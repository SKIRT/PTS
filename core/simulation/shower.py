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
from ..tools.stringify import stringify, tostr
from ..tools import sequences
from ..basics.configuration import parent_type
from ..tools import filesystem as fs
from ..simulation.remote import get_simulations_for_host
from ..tools import introspection
from ..tools import types
from ..basics import containers
from ..tools import formatting as fmt

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

def compare_simulations(*simulations):

    """
    This function ...
    :param simulations:
    :return:
    """

    # Create shower
    shower = SimulationShower()

    # Run
    shower.run(simulations=simulations)

# -----------------------------------------------------------------

def compare_analysis(*simulations):

    """
    This function ...
    :param simulations:
    :return:
    """

    # Create shower
    shower = AnalysisShower()

    # Run
    shower.run(simulations=simulations)

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

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute="id")
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")

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
        return self.simulations[self.simulation_ids[0]]

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
                simulation_id = self.all_simulations[name].id
                self.simulations[simulation_id] = self.all_simulations[name]

        # From IDS
        elif self.config.ids is not None:

            # Load the simulations and put them in the dictionary
            for simulation_id in self.config.ids: self.simulations[simulation_id] = get_simulation_for_host(self.config.remote, simulation_id)

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
                print(" - " + fmt.bold + name + fmt.reset + "(" + description + "): " + pstring + " [" + ptype + "]")

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
            values = get_values_for_simulations(self.simulations, name)

            # Get unique values
            unique_values = sequences.unique_values(values.values())
            nunique_values = len(unique_values)

            if nunique_values == 1:

                ptype, string = stringify(unique_values[0])
                print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + " [" + ptype + "]" + fmt.reset)

            else:

                print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + fmt.reset)

                for value in unique_values:
                    ptype, string = stringify(value)
                    print("   * " + string + " [" + ptype + "]")

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

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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
            if types.is_sequence(simulations) or types.is_tuple(simulations): self.simulations = containers.dict_from_sequence(simulations, attribute="id")
            elif types.is_dictionary(simulations): self.simulations = simulations
            else: raise ValueError("Simulations must be specified as sequence or dictionary")

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
        return self.simulations[self.simulation_ids[0]]

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
            values = get_analysis_values_for_simulations(self.simulations, name)

            # Get unique values
            unique_values = sequences.unique_values(values.values())
            nunique_values = len(unique_values)

            # Get the ptype
            ptype = self.get_ptype(name)

            # Only one unique value
            if len(unique_values) == 1:

                default = unique_values[0]
                ptype, string = stringify(default)

                print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + " [" + ptype + "]" + fmt.reset)

            # Multiple unique values
            else:

                # Prompt to change this property
                #change = prompt_proceed("Change the analysis option '" + name + "' for all simulations? Values are:\n - " + "\n - ".join(tostr(value) for value in unique_values))
                #if not change: continue
                #default = None
                #ptype = get_common_ptype(unique_values)
                #choices = None
                #suggestions = unique_values

                print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "):" + fmt.reset)

                for value in unique_values:

                    ptype, string = stringify(value)
                    print("   * " + string + " [" + ptype + "]")

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
                ptype = self.get_ptype(name, section_name)

                #print(name, ptype, unique_values)

                # Only one unique value
                if nunique_values == 1:

                    default = unique_values[0]
                    ptype, string = stringify(default)

                    print(fmt.green + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "): " + string + " [" + ptype + "]" + fmt.reset)

                # Multiple unique values
                else:

                    print(fmt.red + " - " + fmt.bold + name + fmt.reset_bold + " (" + description + "):" + fmt.reset)

                    for value in unique_values:

                        ptype, string = stringify(value)
                        print("    * " + string + " [" + ptype + "]")

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
