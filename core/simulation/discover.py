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
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..simulation.simulation import SkirtSimulation
from ..simulation.skifile import SkiFile
from ..tools import formatting as fmt
from ..tools.logging import log

# -----------------------------------------------------------------

class SimulationDiscoverer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(SimulationDiscoverer, self).__init__(config)

        # The simulations
        self.simulations = defaultdict(list)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

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

        if len(self.simulations) > 1: raise ValueError("The found simulations correspond to multiple ski files")
        return self.simulations[self.simulations.keys()[0]]

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        search_path = self.config.path

        ski_files = defaultdict(list)
        parameter_files = defaultdict(list)

        # Search for ski files
        for path, name in fs.files_in_path(search_path, extension="ski", recursive=self.config.recursive, returns=["path", "name"]):

            # Determine the prefix
            prefix = name

            dirpath = fs.directory_of(path)

            # Add the path
            ski_files[prefix].append(dirpath)

        # Search for parameter files
        for path, name in fs.files_in_path(search_path, extension="xml", endswith="_parameters", recursive=self.config.recursive, returns=["path", "name"]):

            # Determine the prefix
            prefix = name.split("_parameters")[0]

            dirpath = fs.directory_of(path)

            # Add the path
            parameter_files[prefix].append(dirpath)

        # Search for log files
        for path, name in fs.files_in_path(search_path, extension="txt", endswith="log", recursive=self.config.recursive, returns=["path", "name"]):

            # The simulation output path
            output_path = fs.directory_of(path)

            # Determine the prefix
            prefix = name.split("_log")[0]

            # Load simulation
            simulation = get_simulation(prefix, output_path, ski_files, parameter_files)

            # Add the simulation
            #simulations.append(simulation)
            self.simulations[simulation.ski_path].append(simulation)

    # -----------------------------------------------------------------

    def list(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Listing the found simulations ...")

        # Loop over the ski files
        for ski_path in self.simulations:

            rel_ski_path = ski_path.split(self.config.path)[1]

            print(fmt.green + rel_ski_path + fmt.reset + ":")

            print("")

            counter = 1
            for simulation in self.simulations[ski_path]:

                rel_output_path = simulation.output_path.split(self.config.path)[1]
                rel_input_path = simulation.input_path.split(self.config.path)[1] if simulation.input_path is not None else None

                print(" * " + fmt.bold + "(" + str(counter) + ") " + fs.name(rel_output_path) + fmt.reset + ":")
                print("")
                print("    output: " + rel_output_path)
                if rel_input_path is not None: print("    input: " + rel_input_path)
                print("    processes: " + str(simulation.processes()))
                print("    threads: " + str(simulation.threads()))
                print("")

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
