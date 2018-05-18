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

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self.pbs_options["N"] if "N" in self.pbs_options else None

    # -----------------------------------------------------------------

    @property
    def output_path(self):

        """
        This function ...
        :return:
        """

        return self.pbs_options["o"] if "o" in self.pbs_options else None

    # -----------------------------------------------------------------

    @property
    def error_path(self):

        """
        This property ...
        :return:
        """

        return self.pbs_options["e"] if "e" in self.pbs_options else None

    # -----------------------------------------------------------------

    @property
    def nnodes(self):

        """
        This property ...
        :return:
        """

        if "l" not in self.pbs_options: return None
        for line in self.pbs_options["l"]:
            if line.startswith("nodes="):
                return int(line.split("nodes=")[1].split(":")[0])
        return None

    # -----------------------------------------------------------------

    @property
    def ppn(self):

        """
        This function ...
        :return:
        """

        if "l" not in self.pbs_options: return None
        for line in self.pbs_options["l"]:
            if line.startswith("nodes="):
                return int(line.split("ppn=")[1])
        return None

    # -----------------------------------------------------------------

    @property
    def mail(self):

        """
        This property ...
        :return:
        """

        return "m" in self.pbs_options

    # -----------------------------------------------------------------

    @property
    def walltime(self):

        """
        This property ...
        :return:
        """

        if "l" not in self.pbs_options: return None
        for line in self.pbs_options["l"]:
            if line.startswith("walltime="):
                hours, minutes, seconds = line.split("walltime=")[1].split(":")
                return float(hours) * 3600 + float(minutes) * 60 + float(seconds)
        return None

    # -----------------------------------------------------------------

    @property
    def walltime_seconds(self):

        """
        This function ...
        :return:
        """

        walltime = self.walltime
        if walltime is None: return walltime
        else: return self.walltime * 3600

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, **kwargs):

        """
        This function ...
        :param path:
        :param kwargs:
        :return:
        """

        from ..tools import filesystem as fs

        # Get the lines
        lines = fs.get_lines(path)

        # Create from lines
        script = cls.from_lines(lines, **kwargs)
        script.path = path

        # Return
        return script

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote, **kwargs):

        """
        This function ...
        :param path:
        :param remote:
        :param kwargs:
        :return:
        """

        # Get the lines
        lines = remote.get_lines(path)

        # Create and return
        return cls.from_lines(lines, remote=remote, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_lines(cls, lines, **kwargs):

        """
        This function ...
        :param lines:
        :param kwargs:
        :return:
        """

        # Get the properties
        name, walltime, output_path, error_path, mail, nodes, ppn, header_lines, pbs_lines, modules, commands, extra_header_lines = get_jobscript_properties(lines)

        # Create
        jobscript = cls(name, walltime, nodes, ppn, output_path=output_path, error_path=error_path, mail=mail, extra_header_lines=extra_header_lines)

        # Set the commands and modules
        jobscript.commands = commands
        jobscript.modules = modules

        # Create the jobscript object
        # TODO: FIX API: CONSTRUCTOR ARGUMENTS OF BASE CLASS JOBSCRIPT IS DIFFERENT FROM THAT OF DERIVED CLASS SKIRTJOBSCRIPT:
        # SO THE LINE BELOW FAILS WHEN CREATING SKIRTJOBSCRIPT.FROM_FILE
        # ARGUMENTS FOR SKIRTJOBSCRIPT CONSTRUCTOR: name, arguments, host_id, cluster, skirt_path, mpi_command, walltime, modules, mail=False, bind_to_cores=False, extra_header_lines=None, remote=None
        #jobscript = cls(name, walltime, nodes, ppn, output_path=output_path, error_path=error_path, mail=mail, extra_header_lines=extra_header_lines)
        #jobscript = JobScript(name, walltime, nodes, ppn, output_path=output_path, error_path=error_path, mail=mail, extra_header_lines=extra_header_lines)
        # NOW THE FROM_LINES FUNCTION IS OVERRIDDEN IN THE DERIVED SKIRTJOBSCRIPT CLASS

        # Set the commands and modules
        #jobscript.commands = commands
        #jobscript.modules = modules

        # Return the object
        return jobscript

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

def get_jobscript_properties(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    name = None
    walltime = None
    output_path = None
    error_path = None
    mail = False
    nodes = None
    ppn = None

    header_lines = []
    pbs_lines = []
    modules = []
    commands = []

    # Loop over the lines
    body = False
    last_comment = None
    for line in lines:

        # PBS options
        if line.startswith("#PBS"):
            pbs_lines.append(line.split("#PBS ")[1])
            last_comment = None
            body = True

        # Comment or header
        elif line.startswith("#"):

            if body:
                last_comment = line[2:]
            else:
                header_lines.append(line[2:])
                last_comment = None

        # Module
        elif line.startswith("module load"):
            modules.append(line.split("module load ")[1])
            last_comment = None

        # Other not empty lines
        elif line.strip():
            commands.append((line, last_comment))
            last_comment = None

    # print(header_lines)
    # print(pbs_lines)
    # print(modules)
    # print(commands)

    extra_header_lines = []
    for line in header_lines:
        if line.startswith("/bin/sh"):
            continue
        elif line.startswith("Batch script created with"):
            continue
        else:
            extra_header_lines.append(line)

    for line in pbs_lines:
        if line.startswith("-m"):
            mail = True
        elif line.startswith("-o"):
            output_path = line.split("-o ")[1]
        elif line.startswith("-e"):
            error_path = line.split("-e ")[1]
        elif line.startswith("-N"):
            name = line.split("-N ")[1]
        elif line.startswith("-l"):
            if "walltime=" in line:
                timeformat = line.split("walltime=")[1]
                hours, minutes, seconds = timeformat.split(":")
                walltime = float(hours) + float(minutes) / 60 + float(seconds) / 3600
            elif "nodes=" in line and "ppn=" in line:
                nodes = int(line.split("nodes=")[1].split(":")[0])
                ppn = int(line.split("ppn=")[1])
            else:
                raise ValueError("Unrecognized option: '" + line + "'")
        else: raise ValueError("Unrecognized option: '" + line + "'")

    # Return the properties
    return name, walltime, output_path, error_path, mail, nodes, ppn, header_lines, pbs_lines, modules, commands, extra_header_lines

# -----------------------------------------------------------------
