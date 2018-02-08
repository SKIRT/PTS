#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.jobscript Contains the JobScript class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class JobScript(object):

    """
    This class ...
    """

    def __init__(self, name, walltime, nodes, ppn, output_path=None, error_path=None, mail=False, extra_header_lines=None):

        """
        This function ...
        :param name:
        :param output_path:
        :param error_path:
        :param mail:
        :param extra_header_lines:
        :return:
        """

        # Attributes
        self.header = []
        self.commands = []
        self.modules = []
        self.pbs_options = dict()

        # Determine the walltime in "hours, minutes and seconds" format
        minutes, seconds = divmod(walltime, 60)
        hours, minutes = divmod(minutes, 60)

        # Set PBS options
        self.pbs_options["N"] = name
        if output_path is not None: self.pbs_options["o"] = output_path
        if error_path is not None: self.pbs_options["e"] = error_path
        if mail: self.pbs_options["m"] = "bae"
        self.pbs_options["l"] = []
        self.pbs_options["l"].append("walltime=%d:%02d:%02d" % (hours, minutes, seconds))
        self.pbs_options["l"].append("nodes=" + str(nodes) + ":ppn=" + str(ppn))

        # Set header
        self.header.append("#!/bin/sh")
        self.header.append("# Batch script created with PTS")
        if extra_header_lines is not None:
            for line in extra_header_lines: self.header.append("# " + line)

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

        name = None
        walltime = None
        output_path = None
        error_path = None
        mail = False

        with open(path, 'r') as fh:
            for line in fh:
                line = line[:-1]


        jobscript = cls()

    # -----------------------------------------------------------------

    def import_module(self, module):

        """
        This function ...
        :param module:
        :return:
        """

        self.modules.append(module)

    # -----------------------------------------------------------------

    def add_command(self, command, comment=None):

        """
        This function ...
        :param command:
        :param comment:
        :return:
        """

        self.commands.append((command, comment))

    # -----------------------------------------------------------------

    def to_lines(self):

        """
        This function ...
        :return:
        """

        lines = []

        # Add the header
        for line in self.header: lines.append(line)

        # Empty line
        lines.append("")

        # Add PBS options
        for label in self.pbs_options:
            value = self.pbs_options[label]
            if isinstance(value, list):
                for line in value: lines.append("#PBS -" + label + " " + line)
            else: lines.append("#PBS -" + label + " " + self.pbs_options[label])

        # Empty line
        lines.append("")

        # Import modules
        lines.append("# Load the necessary modules")
        for module in self.modules:
            lines.append("module load " + module)

        # Empty line
        #lines.append("")

        # Add the commands
        for command, comment in self.commands:

            lines.append("")
            if comment is not None: lines.append("# " + comment)
            lines.append(command)

        # Return the list of lines
        return lines

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise RuntimeError("Does not have a path: use the saveto method")
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get lines
        lines = self.to_lines()

        # Write lines
        with open(path, 'w') as fh:
            for line in lines:
                fh.write(line + "\n")

        # Set the new path
        self.path = path

# -----------------------------------------------------------------
