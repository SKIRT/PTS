#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.screen Contains the ScreenScript class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import strings
from .arguments import SkirtArguments
from ..tools.utils import lazyproperty
from ..remote.remote import load_remote
from ..tools import sequences

# -----------------------------------------------------------------

class ScreenScript(object):

    """
    An instance of the ...
    """

    def __init__(self, name, host_id, output_path, process_binding=True, skirt_path=None, mpirun_path=None,
                 remote=None, check_ski=True):

        """
        The constructor ...
        :param name:
        :param host_id:
        :param output_path:
        :param process_binding:
        :param skirt_path:
        :param mpirun_path:
        :param remote:
        :param check_ski:
        :return:
        """

        # Set the screen name
        self.name = name

        # Set paths
        self._skirt_path = skirt_path
        self._mpirun_path = mpirun_path

        # Output path
        self.output_path = output_path

        # Set host ID
        self.host_id = host_id

        # Flags
        self.process_binding = process_binding

        # The arguments
        self.arguments = OrderedDict()

        # Set flag
        self.check_ski = check_ski

        # The path
        self.path = None

        # Reference to the remote
        self._remote = remote

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        script = cls.from_lines(fs.get_lines(path))
        script.path = path
        return script

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        return cls.from_lines(remote.get_lines(path))

    # -----------------------------------------------------------------

    @classmethod
    def from_lines(cls, lines):

        """
        This function ...
        :param lines:
        :return:
        """

        # Get info
        # filename = fs.strip_extension(fs.name(path))
        # queue_name = strings.split_at_last(filename, "_")[0]

        # Initiliaze variables
        host_id = None
        screen_output_path = None
        screen_name = None

        # Get launch info
        header = get_header_lines(lines)
        next_output_path = False
        for line in header:
            if next_output_path:
                screen_output_path = strings.unquote(line.strip())
                next_output_path = False
            elif "running SKIRT on remote host" in line: host_id = strings.unquote(line.split("on remote host ")[1])
            elif "upload this file to the remote filesystem in the following directory" in line: next_output_path = True
            elif line.startswith("screen -S"): screen_name = line.split("screen -S ")[1].split()[0]

        # Create screen script
        screen = cls(screen_name, host_id, screen_output_path, check_ski=False) # check ski will be set to True if conditions are encountered

        # Loop over the simulation lines
        for line in lines:

            # Skip comments and empty lines
            if line.startswith("#"): continue
            if not line: continue
            if line.endswith("--version"): continue

            # Conditional statements
            if line.startswith("if [ -e"):
                screen.check_ski = True
                continue
            elif line.startswith("then"): continue
            elif line.startswith("else"): continue
            elif line.startswith("fi"): continue
            line = line.strip()
            if line.startswith("echo"): continue

            # Get SKIRT arguments
            arguments = SkirtArguments.from_command(line)

            # Check whether the simulation name is defined
            if arguments.simulation_name is None: raise ValueError("Simulation name is not defined")

            # Add
            screen.add_simulation(arguments.simulation_name, arguments)

        # Return the screen script
        return screen

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        if self._remote is not None: return self._remote
        else: return load_remote(self.host_id)

    # -----------------------------------------------------------------

    @property
    def skirt_paths(self):

        """
        This function ...
        :return:
        """

        return [arguments.skirt_path for arguments in self.arguments.values()]

    # -----------------------------------------------------------------

    @property
    def mpirun_paths(self):

        """
        This function ...
        :return:
        """

        return [arguments.mpirun_path for arguments in self.arguments.values()]

    # -----------------------------------------------------------------

    @property
    def ski_paths(self):

        """
        This function ...
        :return:
        """

        return [arguments.ski_pattern for arguments in self.arguments.values()]

    # -----------------------------------------------------------------

    @property
    def skirt_path(self):

        """
        This function ...
        :return:
        """

        if self._skirt_path is not None: return self._skirt_path
        else: return sequences.get_all_equal_value(self.skirt_paths, return_none=True)

    # -----------------------------------------------------------------

    @property
    def mpirun_path(self):

        """
        This function ...
        :return:
        """

        if self._mpirun_path is not None: return self._mpirun_path
        else: return sequences.get_all_equal_value(self.mpirun_paths, return_none=True, ignore_none=True)

    # -----------------------------------------------------------------

    def add_simulation(self, name, arguments):

        """
        This function ...
        :param name:
        :param arguments:
        :return:
        """

        # Check name
        if name in self.simulation_names: raise ValueError("Already a simulation with the name '" + name + "'")

        # Add the arguments
        self.arguments[name] = arguments

    # -----------------------------------------------------------------

    def remove_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.arguments.pop(name)

    # -----------------------------------------------------------------

    def get_simulations_before(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return sequences.before(self.simulation_names, simulation_name)

    # -----------------------------------------------------------------

    def get_simulations_before_ski_path(self, ski_path):

        """
        This function ...
        :param ski_path:
        :return:
        """

        # Get indices before specified ski path
        indices = sequences.indices_before(self.ski_paths, ski_path)

        # Return the simulation names
        names = self.simulation_names
        return [names[index] for index in indices]

    # -----------------------------------------------------------------

    def get_ski_paths_before_ski_path(self, ski_path):

        """
        This function ...
        :param ski_path:
        :return:
        """

        # Get indices before specified ski path
        indices = sequences.indices_before(self.ski_paths, ski_path)

        # Return
        ski_paths = self.ski_paths
        return [ski_paths[index] for index in indices]

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.arguments.keys()

    # -----------------------------------------------------------------

    def has_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.arguments

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_names)

    # -----------------------------------------------------------------

    @property
    def has_simulations(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations > 0

    # -----------------------------------------------------------------

    def get_arguments(self, simulation_name):

        """
        This functio n...
        :param simulation_name:
        :return:
        """

        return self.arguments[simulation_name]

    # -----------------------------------------------------------------

    def get_logging_options(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.arguments[simulation_name].logging

    # -----------------------------------------------------------------

    def get_skirt_path(self, simulation_name):

        """
        Thisf unction ...
        :param simulation_name:
        :return:
        """

        return self.arguments[simulation_name].skirt_path

    # -----------------------------------------------------------------

    def get_mpirun_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.arguments[simulation_name].mpirun_path

    # -----------------------------------------------------------------

    def get_parallelization(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.arguments[simulation_name].parallelization

    # -----------------------------------------------------------------

    def set_parallelization(self, simulation_name, parallelization):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :return:
        """

        self.arguments[simulation_name].parallelization = parallelization

    # -----------------------------------------------------------------

    @property
    def remote_script_file_name(self):

        """
        This function ...
        :return:
        """

        return self.name + ".sh"

    # -----------------------------------------------------------------

    def to_lines(self):

        """
        This function ...
        :return:
        """

        # Initialize a list for the lines
        lines = []

        # Add header
        lines.append("#!/bin/sh")
        lines.append("# Batch script for running SKIRT on remote host " + self.host_id)
        lines.append("# To execute manualy, upload this file to the remote filesystem in the following directory:")
        lines.append("# " + self.output_path)
        lines.append("# under the name '" + self.remote_script_file_name + "' and enter the following commmands:")
        lines.append("# cd '" + self.output_path + "' # navigate to the screen output directory")
        lines.append("# screen -S " + self.name + " -L -d -m " + self.remote_script_file_name + "'")
        lines.append("")

        # Show version of MPI
        if self.mpirun_path is not None:
            lines.append(self.mpirun_path + " --version")
            lines.append("")

        # Loop over the arguments
        for simulation_name in self.simulation_names:

            # Get arguments
            arguments = self.arguments[simulation_name]

            # Get command
            command = arguments.to_command(scheduler=False, skirt_path=self.skirt_path, mpirun_path=self.mpirun_path,
                                           bind_to_cores=self.process_binding, to_string=True, remote=self.remote)

            # Add simulation name to command
            command_line = command + " # " + simulation_name

            # Check presence of ski file?
            if self.check_ski:

                # Add check
                line = "if [ -e '" + arguments.ski_pattern + "' ]"
                lines.append(line)

                lines.append("then")

                line = "    " + command_line
                lines.append(line)

                lines.append("else")

                line = '    echo "Simulation ' + "'" + simulation_name + "' cannot be executed: ski file is missing (" + arguments.ski_pattern + ")" + '"'
                lines.append(line)

                lines.append("fi")

            # Don't check the presence of the ski file
            else: lines.append(command_line)

            # Add empty line
            lines.append("")

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Get the lines
        lines = self.to_lines()

        # Write the lines
        fs.write_lines(path, lines)

        # Update path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether a path is defined for the script file
        if self.path is None: raise RuntimeError("The screen script file does not exist yet")

        # Save to the original path
        self.saveto(self.path)

# -----------------------------------------------------------------

def get_header_lines(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    header_lines = []

    # Loop over the lines
    for line in lines:

        # We are no longer at the header
        if not line.startswith("#"): break

        # Clean line
        line = line.split("#", 1)[1].strip()
        header_lines.append(line)

    # Return the lines
    return header_lines

# -----------------------------------------------------------------