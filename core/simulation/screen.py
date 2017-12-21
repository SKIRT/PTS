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

# -----------------------------------------------------------------

class ScreenScript(object):

    """
    An instance of the ...
    """

    def __init__(self, name, skirt_path, mpirun_path=None):

        """
        The constructor ...
        :param name:
        :param skirt_path:
        :param mpirun_path:
        :return:
        """

        # Set the screen name
        self.name = name

        # Set paths
        self.skirt_path = skirt_path
        self.mpirun_path = mpirun_path

        # The arguments
        self.arguments = OrderedDict()

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

    # -----------------------------------------------------------------

    def add_simulation(self, name, arguments):

        """
        This function ...
        :param name:
        :param arguments:
        :return:
        """

        # Write the command string to the job script
        threads_per_core = self.threads_per_core if self.use_hyperthreading else 1
        command = arguments.to_command(scheduler=False, skirt_path=self.skirt_path, mpirun_path=self.host.mpi_command,
                                       bind_to_cores=self.host.force_process_binding,
                                       threads_per_core=threads_per_core, to_string=True, remote=self)

        script_file.write(command + "\n")
        # Write to disk
        script_file.flush()

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # If a path is given, create a script file at the specified location
        script_file = open(local_script_path, 'w')

        # If no screen output path is set, create a directory
        if screen_output_path is None: screen_output_path = self.create_directory_in(self.pts_temp_path, screen_name, recursive=True)

        # Write a general header to the batch script
        remote_script_file_name = screen_name + ".sh"
        script_file.write("#!/bin/sh\n")
        script_file.write("# Batch script for running SKIRT on remote host " + self.host_id + "\n")
        script_file.write("# To execute manualy, upload this file to the remote filesystem in the following directory:\n")
        script_file.write("# " + screen_output_path + "\n")
        script_file.write("# under the name '" + remote_script_file_name + "' and enter the following commmands:\n")
        script_file.write("# cd '" + screen_output_path + "' # navigate to the screen output directory\n")
        script_file.write("# screen -S " + screen_name + " -L -d -m " + remote_script_file_name + "'\n")
        script_file.write("\n")

        # Close the script file (if it is temporary it will automatically be removed)
        script_file.close()

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
