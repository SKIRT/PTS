#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.discover Contains the SimulationDiscoverer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import math
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..simulation.simulation import SkirtSimulation
from ..simulation.skifile import SkiFile
from ..tools import formatting as fmt
from ..basics.log import log
from ..simulation.logfile import LogFile
from ..basics.map import Map
from ..tools import types, numbers

# -----------------------------------------------------------------

def find_one_simulation_in_path(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    # Check paths
    ski_path = fs.find_file_in_path(path, extension="ski")
    in_path = fs.join(path, "in")
    out_path = fs.join(path, "out")

    # Determine ski prefix
    prefix = fs.strip_extension(fs.name(ski_path))

    if not fs.is_directory(in_path): in_path = None  # no input required
    if not fs.is_directory(out_path):
        log_path = fs.join(path, prefix + "_log.txt")
        out_path = path
    else: log_path = fs.join(out_path, prefix + "_log.txt")

    if not fs.is_file(log_path): raise IOError("Log file is not found for simulation '" + prefix + "'")

    # Return relevant stuff
    return prefix, ski_path, in_path, out_path

# -----------------------------------------------------------------

class SimulationDiscoverer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SimulationDiscoverer, self).__init__(*args, **kwargs)

        # The search paths
        self.search_paths = []

        # The paths
        self.ski_paths = defaultdict(list)
        self.parameter_paths = dict()
        self.log_paths = defaultdict(list)

        # The simulations
        self.simulations_ski = defaultdict(list)
        self.simulations_no_ski = defaultdict(list)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Find the simulations
        self.find()

        # List the simulations
        if self.config.list: self.list()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    @property
    def simulations_single_ski(self):

        """
        This function ...
        :return:
        """

        if len(self.simulations_ski) == 0 and len(self.simulations_no_ski) == 0: raise ValueError("No simulations found")
        elif len(self.simulations_ski) == 1 and len(self.simulations_no_ski) == 0:
            return self.simulations_ski[self.simulations_ski.keys()[0]]
        elif len(self.simulations_ski) == 0 and len(self.simulations_no_ski) == 1:
            return self.simulations_no_ski[self.simulations_no_ski.keys()[0]]
        else: raise ValueError("The found simulations correspond to multiple ski files")

    # -----------------------------------------------------------------

    @property
    def simulations(self):

        """
        This function ...
        :return:
        """

        simulations = []

        for key in self.simulations_ski: simulations += self.simulations_ski[key]
        for key in self.simulations_no_ski: simulations += self.simulations_no_ski[key]

        return simulations

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationDiscoverer, self).setup(**kwargs)

        # Set the search path
        if self.config.directories is None: self.search_paths = [self.config.path]
        else:

            # Loop over the directories
            for name in self.config.directories:

                # Determine the full path and add it to the list of search paths
                path = fs.absolute_or_in(name, self.config.path)
                self.search_paths.append(path)

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding simulations ...")

        # Find ski files
        self.find_ski_files()

        # Find parameter files
        self.find_parameter_files()

        # Find log files
        self.find_log_files()

        # Find simulations
        self.find_simulations()

    # -----------------------------------------------------------------

    def find_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding ski files ...")

        # Loop over the search paths
        for search_path in self.search_paths:

            # Search for ski files
            for path, name in fs.files_in_path(search_path, extension="ski", recursive=self.config.recursive, returns=["path", "name"]):

                # Determine the prefix
                prefix = name

                dirpath = fs.directory_of(path)

                # Add the path
                self.ski_paths[prefix].append(dirpath)

    # -----------------------------------------------------------------

    def find_parameter_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding parameter files ...")

        parameter_files = defaultdict(list)

        # Loop over the search paths
        for search_path in self.search_paths:

            # Search for parameter files
            for path, name in fs.files_in_path(search_path, extension="xml", endswith="_parameters", recursive=self.config.recursive, returns=["path", "name"]):

                # Determine the prefix
                prefix = name.split("_parameters")[0]

                dirpath = fs.directory_of(path)

                # Add the path
                parameter_files[prefix].append(dirpath)

        #parpaths = dict()

        # Compare the parameter files
        for prefix in parameter_files:

            distinct_parpaths = defaultdict(list)

            for dirpath in parameter_files[prefix]:

                parpath = fs.join(dirpath, prefix + "_parameters.xml")

                if len(distinct_parpaths) == 0: distinct_parpaths[parpath].append(parpath)
                else:

                    for path in distinct_parpaths:

                        if equaltextfiles(parpath, path, 1):
                            distinct_parpaths[path].append(parpath)
                            break

                    else: distinct_parpaths[parpath].append(parpath)

            # Set ...
            self.parameter_paths[prefix] = distinct_parpaths

    # -----------------------------------------------------------------

    def find_log_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding log files ...")

        # Loop over the search paths
        for search_path in self.search_paths:

            # Search for log files
            for path, name in fs.files_in_path(search_path, extension="txt", endswith="log", recursive=self.config.recursive, returns=["path", "name"]):

                # Determine the prefix
                prefix = name.split("_log")[0]

                # The simulation output path
                #output_path = fs.directory_of(path)

                dirpath = fs.directory_of(path)

                self.log_paths[prefix].append(dirpath)

    # -----------------------------------------------------------------

    def find_simulations(self):

        """
        This function ...
        :return:
        """

        input_paths_ski_files = dict()

        # Loop over the ...
        for prefix in self.parameter_paths:

            # Loop over the unique parameter files, find equal ski file
            for parpath in self.parameter_paths[prefix]:

                #rel_parpath = parpath.split(self.config.path)[1]
                ski_found = False

                # Loop over the ski files with this prefix
                for dirpath in self.ski_paths[prefix]:

                    skipath = fs.join(dirpath, prefix + ".ski")
                    rel_skipath = skipath.split(self.config.path)[1]

                    #print("DIFFERENCES", parpath, skipath)
                    #for line1, line2 in textfiledifferences(parpath, skipath):
                    #    print(line1, line2)

                    #output = textfiledifferences(parpath, skipath)
                    #print(output)

                    par = SkiFile(parpath)
                    ski = SkiFile(skipath)

                    #print(set(par.tree.getroot().text))
                    #print(set(ski.tree.getroot().text))

                    #par_string = str(par)
                    #ski_string = str(ski)

                    #par_string = simplify_xml(par_string)
                    #ski_string = simplify_xml(ski_string)

                    #print(par_string)
                    #print(ski_string)

                    #par_string = str(par).split("MonteCarloSimulation")[1]
                    #ski_string = str(ski).split("MonteCarloSimulation")[1]

                    #exit()

                    #par_dict = par.to_dict()
                    #ski_dict = ski.to_dict()

                    #print(par_dict)
                    #print(ski_dict)
                    #exit()

                    # Compare ski file with parameter file
                    if equaltextfiles(parpath, skipath, 1):
                        key = skipath
                        ski_found = True
                        break

                # String of relative parameter file paths concated together
                #else: key = " + ".join(map(itemgetter(1), map(lambda x: str.split(x, self.config.config_path), parpaths[prefix][parpath]))) #
                #else: key = "no ski file found"
                else: key = tuple(self.parameter_paths[prefix][parpath])

                # Loop over the parameter files, find the corresponding log file
                for parameter_file_path in self.parameter_paths[prefix][parpath]:

                    # The simulation output path
                    output_path = fs.directory_of(parameter_file_path)

                    # Get the input path
                    input_path = find_input(parameter_file_path, output_path)

                    # Determine ski file path
                    ski_path = key if ski_found else parameter_file_path

                    # Determine the log file path
                    #log_path = parameter_file_path.replace("parameters.xml", "log.txt")

                    if ski_found: input_paths_ski_files[ski_path] = input_path

                    self.log_paths[prefix].remove(output_path)

                    # Create simulation and return it
                    simulation = SkirtSimulation(prefix, input_path, output_path, ski_path)

                    # Add the simulation
                    self.simulations_ski[key].append(simulation)

        # Loop over the log files without parameter file
        for prefix in self.log_paths:

            parameters_ski_files = dict()

            # Open the ski files with this prefix
            for dirpath in self.ski_paths[prefix]:

                skipath = fs.join(dirpath, prefix + ".ski")

                # Get the input path
                #if prefix in input_paths_ski_files:
                input_path = input_paths_ski_files[skipath] if skipath in input_paths_ski_files else None
                #else: input_path = None
                parameters = comparison_parameters_from_ski(skipath, input_path=input_path)

                # Set ...
                parameters_ski_files[skipath] = parameters

            # Loop over the log files
            for output_path in self.log_paths[prefix]:

                # Determine the log file path
                logfile_path = fs.join(output_path, prefix + "_log.txt")

                # Get the parameters
                log_parameters = comparison_parameters_from_log(logfile_path)

                # Determine the number of processes
                processes = LogFile(logfile_path).processes

                # Loop over the ski files
                for skipath in parameters_ski_files:

                    #print(parameters_ski_files[skipath])
                    #print(log_parameters)

                    if log_parameters.ncells is None: print(output_path)

                    npackages_reference = parameters_ski_files[skipath].npackages
                    npackages = log_parameters.npackages

                    log_parameters.npackages = npackages_reference
                    #parameters_ski_files[skipath].npackages = None

                    if parameters_ski_files[skipath] == log_parameters and matching_npackages(npackages_reference, npackages, processes):

                        input_path = input_paths_ski_files[skipath] if skipath in input_paths_ski_files else None

                        # Create simulation and return it
                        simulation = SkirtSimulation(prefix, input_path, output_path, skipath)

                        # Add the simulation
                        key = skipath

                        # Add the simulation
                        self.simulations_ski[key].append(simulation)

                        break

                # No matching ski file found
                else:

                    # Create tuple for the parameters
                    properties = (log_parameters.npackages, log_parameters.nwavelengths, log_parameters.treegrid,
                                  log_parameters.ncells, log_parameters.npopulations, log_parameters.selfabsorption,
                                  log_parameters.transient_heating)

                    # Get the input path
                    #input_path = find_input(parameter_file_path, output_path)

                    # Create simulation and return it
                    simulation = SkirtSimulation(prefix, None, output_path, None, parameters=log_parameters)

                    # Match the npackages with properties of previous simulations without parameters
                    for key in self.simulations_no_ski:

                        npackages_reference = key[0]
                        npackages = log_parameters.npackages

                        #log_parameters.npackages = npackages_reference
                        properties = list(properties)
                        properties[0] = npackages_reference
                        properties = tuple(properties)

                        if key == properties and matching_npackages(npackages_reference, npackages, processes):
                            properties = key
                            break

                    # Add the simulation
                    #self.simulations_no_parameters.append(simulation)
                    self.simulations_no_ski[properties].append(simulation)

    # -----------------------------------------------------------------

    def list(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Listing the found simulations ...")

        # Loop over the ski files
        for ski_path in self.simulations_ski:

            if types.is_string_type(ski_path):

                rel_ski_path = ski_path.split(self.config.path)[1]
                print(fmt.green + rel_ski_path + fmt.reset + ":")

                input_path = self.simulations_ski[ski_path][0].input_path
                parameters = comparison_parameters_from_ski(ski_path, input_path)

            else:

                nsimulations = len(self.simulations_ski[ski_path])
                if nsimulations == 1: print(fmt.yellow + "ski file not found" + fmt.reset + ":")
                else: print(fmt.yellow + "ski file not found (but identical parameters)" + fmt.reset + ":")

                input_path = self.simulations_ski[ski_path][0].input_path
                parameters = comparison_parameters_from_ski(ski_path[0], input_path)

            print("")
            print(" npackages:", parameters.npackages)
            print(" nwavelengths:", parameters.nwavelengths)
            print(" treegrid:", parameters.treegrid)
            print(" ncells:", parameters.ncells)
            print(" npopulations:", parameters.npopulations)
            print(" selfabsorption:", parameters.selfabsorption)
            print(" transient heating:", parameters.transient_heating)
            print("")

            counter = 1
            for simulation in self.simulations_ski[ski_path]:

                rel_output_path = simulation.output_path.split(self.config.path)[1]
                rel_input_path = simulation.input_path.split(self.config.path)[1] if simulation.input_path is not None else None

                if rel_input_path == "": rel_input_path = "this directory"
                if rel_output_path == "": rel_output_path = "this directory"

                print(" * " + fmt.bold + "(" + str(counter) + ") '" + simulation.prefix() + "' in " + fs.name(rel_output_path) + fmt.reset + ":")
                print("")
                print("    output: " + rel_output_path)
                if rel_input_path is not None: print("    input: " + rel_input_path)
                print("    processes: " + str(simulation.processes()))
                print("    threads: " + str(simulation.threads()))
                print("    data-parallel: " + str(simulation.log_file.data_parallel))
                print("    host: " + simulation.log_file.host)
                print("")

                if self.config.output:
                    print("    output:")
                    simulation.output.show(line_prefix="      ")

                counter += 1

        #if len(self.simulations_no_parameters) > 0:
        #    print(fmt.red + "no parameters found" + fmt.reset + ":")
        #    print("")

        for properties in self.simulations_no_ski:

            print(fmt.red + "no parameter file found" + fmt.reset + ":")
            print("")

            print(" npackages:", properties[0])
            print(" nwavelengths:", properties[1])
            print(" treegrid:", properties[2])
            print(" ncells:", properties[3])
            print(" npopulations:", properties[4])
            print(" selfabsorption:", properties[5])
            print(" transient heating:", properties[6])

            print("")

            counter = 1
            for simulation in self.simulations_no_ski[properties]:

                rel_output_path = simulation.output_path.split(self.config.path)[1]

                if rel_output_path == "": rel_output_path = "this directory"

                print(" * " + fmt.bold + "(" + str(counter) + ") '" + simulation.prefix() + "' in " + fs.name(rel_output_path) + fmt.reset + ":")
                print("")
                print("    output: " + rel_output_path)
                #if rel_input_path is not None: print("    input: " + rel_input_path)
                print("    processes: " + str(simulation.processes()))
                print("    threads: " + str(simulation.threads()))
                print("    data-parallel: " + str(simulation.log_file.data_parallel))
                print("    host: " + simulation.log_file.host)
                print("")

                if self.config.output:
                    print("    output:")
                    simulation.output.show(line_prefix="      ")

                counter += 1

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def get_simulation(prefix, output_path, ski_files, parameter_files):

    """
    This function returns a simulation object based on the prefix, the output directory, and information about ski files and parameter files
    that are found on the filesystem
    """

    if prefix in ski_files:

        # Check if the directory is present
        if output_path in ski_files[prefix]:

            # Determine the ski path
            ski_path = fs.join(output_path, prefix + ".ski")

            # Get the input path
            input_path = find_input(ski_path, output_path)

            # Create simulation and return it
            return SkirtSimulation(prefix, input_path, output_path, ski_path)

    if prefix in ski_files and len(ski_files[prefix]) == 1:

        # Determine the ski path
        ski_path = fs.join(ski_files[prefix][0], prefix + ".ski")

        # Get the input path
        input_path = find_input(ski_path, output_path)

        # Create the simulation and return it
        return SkirtSimulation(prefix, input_path, output_path, ski_path)

    if prefix in parameter_files:

        # Check if the directory is present
        if output_path in parameter_files[prefix]:

            # Determine the parameters.xml path
            parameters_path = fs.join(output_path, prefix + "_parameters.xml")

            # Get the input path
            input_path = find_input(parameters_path, output_path)

            # Create simulation and return it
            return SkirtSimulation(prefix, input_path, output_path, parameters_path)

    # Error
    raise RuntimeError("Information needed for simulation '" + prefix + "' in '" + output_path + "' not found")

# -----------------------------------------------------------------

def find_input(ski_path, output_path):

    """
    This function tries to determined the input directory based on a ski path and an output path
    :param ski_path:
    :param output_path:
    :return:
    """

    # Load the ski file
    ski = SkiFile(ski_path)

    # If the ski file needs input
    if ski.needs_input:

        # Loop over the input files
        for filename in ski.input_files:
            filepath = fs.join(output_path, filename)
            if not fs.is_file(filepath): break
        else: return output_path  # input_path = output_path

        # Look in subdirectories of the output path if this path is also the directory of the ski file
        if output_path == fs.directory_of(ski_path):

            for path in fs.directories_in_path(output_path):

                for filename in ski.input_files:
                    if not fs.contains_file(path, filename): break
                else: return path # input_path = path

        # Look in subdirectories of the ski file directory if this is only one directory 'up' with respect to the output directory
        elif fs.directory_of(output_path) == fs.directory_of(ski_path):

            for path in fs.directories_in_path(fs.directory_of(ski_path), exact_not_name=fs.name(output_path)):

                for filename in ski.input_files:
                    if not fs.contains_file(path, filename): break
                else: return path

        # Input not found
        else: raise RuntimeError("Input directory could not be found")

    # Needs no input
    else: return None

# -----------------------------------------------------------------

def equaltextfiles(filepath1, filepath2, max_ndiffs, ignore_empty=True):

    """
    This function ...
    :param filepath1:
    :param filepath2:
    :param max_ndiffs:
    :param ignore_empty:
    :return:
    """

    # Get lines
    lines1 = fs.get_lines(filepath1)
    lines2 = fs.get_lines(filepath2)

    if ignore_empty:

        lines1 = [line for line in lines1 if line]
        lines2 = [line for line in lines2 if line]

    return equallists(lines1, lines2, max_ndiffs)

# -----------------------------------------------------------------

def textfiledifferences(filepath1, filepath2):

    """
    This function ...
    :param filepath1:
    :param filepath2:
    :return:
    """

    # Get lines
    lines1 = fs.get_lines(filepath1)
    lines2 = fs.get_lines(filepath2)

    indices = listdifferences(lines1, lines2)

    # Return the different lines
    return [(lines1[index], lines2[index]) for index in indices]

# -----------------------------------------------------------------

def equallists(list1, list2, max_ndiffs):

    """
    This function ...
    :param list1:
    :param list2:
    :param max_ndiffs:
    :return:
    """

    # The lists must have the same length, which must be at least 2 (to avoid everything being read into 1 line)
    length = len(list1)
    if length < 2 or length != len(list2): return False

    # Compare the lists item by item
    diffs = 0
    for index in range(length):

        if list1[index] != list2[index]:

            # verify against allowed number of differences
            diffs += 1

            if diffs > max_ndiffs: return False

            # verify that the differing items are identical up to numerics and months
            pattern = re.compile(r"(\d{1,4}-[0-9a-f]{7,7}(-dirty){0,1})|\d{1,4}|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec")
            item1 = re.sub(pattern, "*", list1[index])
            item2 = re.sub(pattern, "*", list2[index])

            item1 = item1.replace(" k", "000 ")
            item2 = item2.replace(" k", "000 ")

            if item1.strip() != item2.strip(): return False

    # no relevant differences
    return True

# -----------------------------------------------------------------

def listdifferences(list1, list2):

    """
    This function ...
    :param list1:
    :param list2:
    :return:
    """

    # The lists must have the same length, which must be at least 2 (to avoid everything being read into 1 line)
    length = min(len(list1), len(list2))
    #if length < 2 or length != len(list2): return [-1]

    indices = []

    # compare the lists item by item
    for index in range(length):

        if list1[index] != list2[index]:

            # verify that the differing items are identical up to numerics and months
            pattern = re.compile(r"(\d{1,4}-[0-9a-f]{7,7}(-dirty){0,1})|\d{1,4}|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec")
            item1 = re.sub(pattern, "*", list1[index])
            item2 = re.sub(pattern, "*", list2[index])

            item1 = item1.replace(" k", "000 ")
            item2 = item2.replace(" k", "000 ")

            if item1.strip() != item2.strip(): indices.append(index)

    # Return the indices
    return indices

# -----------------------------------------------------------------

def comparison_parameters_from_ski(ski_path, input_path=None):

    """
    This function ...
    :param ski_path:
    :param input_path:
    """

    # Open the ski file
    ski = SkiFile(ski_path)

    # Get the number of wavelengths
    if ski.wavelengthsfile():
        if ski.find_wavelengthsfile(input_path): nwavelengths = ski.nwavelengthsfile(input_path)
        else: nwavelengths = None
    else: nwavelengths = ski.nwavelengths()

    # Get the dust grid type
    #grid_type = ski.gridtype()

    # If the grid is a tree grid, get additional properties
    if ski.treegrid():

        ncells = None

        if ski.filetreegrid(): pass
        else:

            min_level = ski.tree_min_level()
            max_level = ski.tree_max_level()
            search_method = ski.tree_search_method()
            sample_count = ski.tree_sample_count()
            max_optical_depth = ski.tree_max_optical_depth()
            max_mass_fraction = ski.tree_max_mass_fraction()
            max_dens_disp = ski.tree_max_dens_disp()

    # Else, set all properties to None
    else:

        ncells = ski.ncells()
        min_level = max_level = search_method = sample_count = max_optical_depth = max_mass_fraction = max_dens_disp = None

    # Check whether dust self-absorption was enabled for the simulation
    selfabsorption = ski.dustselfabsorption()

    # Check whether transient heating was enabled for the simulation
    transient_heating = ski.transientheating()

    # Determine the total number of pixels from all the instruments defined in the ski file
    npixels = ski.nspatialpixels()

    # Get the number of dust populations
    npopulations = ski.npopulations()

    # Set the parameters
    parameters = Map()
    parameters.npackages = ski.packages()
    parameters.nwavelengths = nwavelengths
    parameters.treegrid = ski.treegrid()
    parameters.ncells = ncells
    parameters.npopulations = npopulations
    parameters.selfabsorption = selfabsorption
    parameters.transient_heating = transient_heating

    # Return the parameters
    return parameters

# -----------------------------------------------------------------

def comparison_parameters_from_log(log_path):

    """
    This function ...
    :param log_path:
    :return:
    """

    # Open the log file
    log_file = LogFile(log_path)

    # Set the parameters
    log_parameters = Map()

    # Get the properties
    log_parameters.npackages = log_file.stellar_packages
    log_parameters.nwavelengths = log_file.wavelengths
    log_parameters.treegrid = log_file.uses_tree
    log_parameters.ncells = log_file.dust_cells if not log_file.uses_tree else None
    log_parameters.npopulations = log_file.npopulations
    log_parameters.selfabsorption = log_file.selfabsorption
    log_parameters.transient_heating = log_file.uses_transient_heating

    # Return the parameters
    return log_parameters

# -----------------------------------------------------------------

def simplify_xml(xml_string):

    new_string = xml_string
    index = 0
    last_name = None
    left_or_right = None
    last_end_index = None
    while index < len(xml_string):
        if xml_string[index] == "<" and xml_string[index + 1] != "/":
            last_name = xml_string[index:].split(" ")[0]
            left_or_right = "<"
        elif xml_string[index] == ">" and xml_string[index + 1] != "/":
            last_end_index = index
        elif xml_string[index] == "<" and xml_string[index + 1] == "/":
            name = xml_string[index + 1:].split(">")[0]
            if name == last_name:
                new_string = new_string[:last_end_index - 1] + "/>" + new_string[index + len(name) + 1:]
    return new_string

# -----------------------------------------------------------------

def matching_npackages(npackages1, npackages2, nprocesses, nprocesses2=None):

    """
    This function ...
    :param npackages1:
    :param npackages2:
    :param nprocesses:
    :param nprocesses2:
    :return:
    """

    if nprocesses2 is None: nprocesses2 = nprocesses

    # Checks
    if not types.is_integer_type(nprocesses): raise ValueError("Number of processes must be integer")
    if not types.is_integer_type(nprocesses2): raise ValueError("Number of processes must be integer")
    if not numbers.is_integer(npackages1): raise ValueError("Number of photon packages must be integer") # doesn't need to be a literal integer (e.g. 1e6 is fine = float)
    if not numbers.is_integer(npackages2): raise ValueError("Number of photon packages must be integer") # doesn't need to be a literal integer

    # Return
    return math.ceil(npackages1 / float(nprocesses)) * nprocesses == math.ceil(npackages2 / float(nprocesses2)) * nprocesses2

# -----------------------------------------------------------------
