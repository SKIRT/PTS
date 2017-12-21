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

# -----------------------------------------------------------------

class ScreenScript(object):

    """
    An instance of the ...
    """

    def __init__(self, name, host_id, output_path, process_binding=True, skirt_path=None, mpirun_path=None, remote=None):

        """
        The constructor ...
        :param name:
        :param host_id:
        :param output_path:
        :param process_binding:
        :param skirt_path:
        :param mpirun_path:
        :param remote:
        :return:
        """

        # Set the screen name
        self.name = name

        # Set paths
        self.skirt_path = skirt_path
        self.mpirun_path = mpirun_path

        # Output path
        self.output_path = output_path

        # Set host ID
        self.host_id = host_id

        # Flags
        self.process_binding = process_binding

        # The arguments
        self.arguments = OrderedDict()

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

        # Get info
        # filename = fs.strip_extension(fs.name(path))
        # queue_name = strings.split_at_last(filename, "_")[0]

        # Initiliaze variables
        host_id = None
        screen_output_path = None
        screen_name = None

        # Get launch info
        header = fs.get_header_lines(path)
        next_output_path = False
        for line in header:
            if next_output_path:
                screen_output_path = strings.unquote(line.strip())
                next_output_path = False
            elif "running SKIRT on remote host" in line: host_id = strings.unquote(line.split("on remote host ")[1])
            elif "upload this file to the remote filesystem in the following directory" in line: next_output_path = True
            elif line.startswith("screen -S"): screen_name = line.split("screen -S ")[1].split()[0]

        # Create screen script
        screen = cls(screen_name, host_id, screen_output_path)

        # Loop over the simulation lines
        for line in fs.read_lines(path):

            # Skip comments and empty lines
            if line.startswith("#"): continue
            if not line: continue

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

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.arguments.keys()

    # -----------------------------------------------------------------

    @property
    def remote_script_file_name(self):

        """
        This function ...
        :return:
        """

        return self.name + ".sh"

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
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
            command = arguments.to_command(scheduler=False, skirt_path=self.skirt_path, mpirun_path=self.mpirun_path, bind_to_cores=self.process_binding, to_string=True, remote=self.remote)

            # Add simulation name to command
            line = command + " # " + simulation_name

            # Add line
            lines.append(line)

        # Write the lines
        fs.write_lines(path, lines)

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
