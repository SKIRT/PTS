#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.python Contains the RemotePythonSession class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules

# Import the relevant PTS classes and modules
from pts.core.tools import time
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log

# -----------------------------------------------------------------

class RemotePythonSession(object):

    """
    This class ...
    """

    def __init__(self, remote, assume_pts=False):

        """
        This function ...
        :param remote:
        :param assume_pts:
        :return:
        """

        # The remote connection
        self.remote = remote

        # The last number of lines in the session output
        self.previous_length = 0

        # Generate session ID and screen name
        self.session_id = timestamp = time.unique_name()
        self.screen_name = "pts_remotepython_" + self.session_id

        # Create pipe file
        out_pipe_filename = "out_pipe_" + self.session_id + ".txt"
        self.out_pipe_filepath = fs.join(self.remote.home_directory, out_pipe_filename)

        # Create the pipe file
        if self.remote.is_file(self.out_pipe_filepath): self.remote.remove_file(self.out_pipe_filepath)
        self.remote.touch(self.out_pipe_filepath)

        # Start the screen
        start_screen_command = "screen -dmS " + self.screen_name
        self.remote.execute(start_screen_command)
        #print("START SCREEN COMMAND: ", start_screen_command)

        # Start python in the screen session
        start_python_command = "screen -r " + self.screen_name + " -X stuff $'python > " + self.out_pipe_filepath + "\n'"
        #print("START PYTHON COMMAND: ", start_python_command)
        # Start python
        self.remote.execute(start_python_command)

        # Import PTS stuff
        if assume_pts:

            # Import standard PTS tools
            self.import_package("filesystem", as_name="fs", from_name="pts.core.tools")

            # Set logging level to match that of local PTS
            if log.is_debug():

                self.import_package("setup_log", from_name="pts.core.tools.logging")
                self.send_line("log = setup_log('DEBUG')")

            else: self.import_package("log", from_name="pts.core.tools.logging")

    # -----------------------------------------------------------------

    @classmethod
    def from_host_id(cls, host_id, assume_pts=False):

        """
        This function ...
        :param host_id:
        :param assume_pts:
        :return:
        """

        from .remote import Remote
        remote = Remote()
        remote.setup(host_id)
        return cls(remote, assume_pts=assume_pts)

    # -----------------------------------------------------------------

    def import_package(self, name, as_name=None, from_name=None):

        """
        This function ...
        :param name:
        :param as_name:
        :param from_name:
        :return:
        """

        command = ""

        if from_name is not None:
            command += "from " + from_name + " "

        command += "import " + name

        if as_name is not None:
            command += " as " + as_name

        # Execute the import command
        #output = self.send_python_line(command, output=True, show_output=log.is_debug())
        output = self.send_line(command, output=True)

        # If output is given, this is normally not so good
        if len(output) > 0:

            # Check output
            last_line = output[-1]
            if "cannot import" in last_line: log.warning(last_line)
            if "ImportError" in last_line: log.warning(last_line)

            return False

        return True

    # -----------------------------------------------------------------

    def remove_variable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.send_line("del " + name)

    # -----------------------------------------------------------------

    def define_simple_variable(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        command = name + " = " + str(value)
        self.send_line(command)

    # -----------------------------------------------------------------

    def get_attributes(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simple_variable("vars(" + name + ")")

    # -----------------------------------------------------------------

    def get_simple_variable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output = self.send_line(name, output=True)

        #if len(output) < 1: raise NameError("No such variable: '" + name + "'")
        if len(output) == 0: return None
        elif len(output) > 1: raise RuntimeError("Unexpected output: " + str(output))

        return eval(output[0])

    # -----------------------------------------------------------------

    def get_simple_property(self, variable, name):

        """
        This function ...
        :param variable:
        :param name:
        :return:
        """

        try: return self.get_simple_variable(variable + "." + name)
        except NameError: raise AttributeError("Variable '" + variable + "' has no attribute '" + name + "'")

    # -----------------------------------------------------------------

    def get_string(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output = self.send_line(name, output=True)
        return output[0][1:-1]

    # -----------------------------------------------------------------

    def variables(self):

        """
        This function ...
        :return:
        """

        # Get the complete dir() list of variables
        variables = self.get_simple_variable("dir()")

        # '__builtins__', '__doc__', '__name__', '__package__'
        variables.remove("__builtins__")
        variables.remove("__doc__")
        variables.remove("__name__")
        variables.remove("__package__")

        # Return the list of user-defined variables
        return variables

    # -----------------------------------------------------------------

    def send_line(self, line, output=True):

        """
        This function ...
        :param line:
        :return:
        """

        return self.execute(line, output=output)

    # -----------------------------------------------------------------

    def execute(self, line, output=True):

        """
        This function ...
        :param line:
        :param output:
        :return:
        """

        send_command = "screen -r " + self.screen_name + " -X stuff $'" + line + "\n'"

        #print("COMMAND: ", send_command)

        # Send the line
        self.remote.execute(send_command, output=False)

        # Check the output
        if output:

            # Get output
            output = self.remote.read_text_file(self.out_pipe_filepath)
            lines = output[self.previous_length:]
            self.previous_length = len(output)
            return lines

    # -----------------------------------------------------------------

    def send_lines(self, lines, show_output=False):

        """
        This function ...
        :param lines:
        :param show_output:
        :return:
        """

        output_lines = []

        # Send the lines consecutively
        for line in lines:

            # Send the line
            output = self.send_line(line, output=True)

            output_lines.append(output)

            # Show output
            if show_output:
                for o in output: print(o)

        return output_lines

    # -----------------------------------------------------------------

    def do_loop(self, top_statement, in_loop_lines, show_output=False):

        """
        This function ...
        :param top_statement:
        :param in_loop_lines:
        :param show_output:
        :return:
        """

        # Top statement
        self.send_line(top_statement)

        # In-loop lines
        for line in in_loop_lines:

            # Send the line
            output = self.send_line("    " + line, output=True)

            # Show output
            if show_output:
                for o in output: print(o)

        # Finish the loop
        self.send_line("")

    # -----------------------------------------------------------------

    def rename_file(self, directory, old_name, new_name):

        """
        This function ...
        :param directory:
        :param old_name:
        :param new_name:
        :return:
        """

        self.send_line("fs.rename_file('" + directory + "', '" + old_name + "', '" + new_name + "')")

    # -----------------------------------------------------------------

    def remove_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.send_line("fs.remove_directory('" + path + "')")

    # -----------------------------------------------------------------

    def remove_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.send_line("fs.remove_file('" + path + "')")

    # -----------------------------------------------------------------

    def change_cwd(self, path):

        """
        THis function ...
        :param path:
        :return:
        """

        self.send_line("fs.change_cwd('" + path + "')")

    # -----------------------------------------------------------------

    def directories_in_path(self, path, startswith=None, recursive=False):

        """
        This function ...
        :param path:
        :param startswith:
        :param recursive:
        :return:
        """

        return self.get_simple_variable("fs.directories_in_path('" + path + "', startswith=" + repr(startswith) + ", recursive=" + repr(recursive) + ")")

    # -----------------------------------------------------------------

    def files_in_path(self, path, recursive=False):

        """
        This function ...
        :return:
        """

        return self.get_simple_variable("fs.files_in_path('" + path + "', recursive=" + repr(recursive) + ")")

    # -----------------------------------------------------------------

    def to_home_directory(self):

        """
        This function ...
        :return:
        """

        self.send_line("fs.to_home_directory()")

    # -----------------------------------------------------------------

    def create_directory(self, path):

        """
        This function ...
        :return:
        """

        self.send_line("fs.create_directory('" + path + "')")

    # -----------------------------------------------------------------

    def create_directories(self, *paths):

        """
        This function ...
        :param paths:
        :return:
        """

        self.send_line("fs.create_directories(*" + str(paths) + ")")

    # -----------------------------------------------------------------

    def file_or_directory(self, path):

        """
        This fucntion ...
        :param path:
        :return:
        """

        return self.get_simple_variable("fs.file_or_directory('" + path + "')")

    # -----------------------------------------------------------------

    def is_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the corresponding PTS function
        return self.get_simple_variable("fs.is_directory('" + path + "')")

    # -----------------------------------------------------------------

    def is_file(self, path):

        """
        THis function ...
        :param path:
        :return:
        """

        # Call the corresponding PTS function
        return self.get_simple_variable("fs.is_file('" + path + "')")

    # -----------------------------------------------------------------

    def is_subdirectory(self, path, parent_path):

        """
        This function ...
        :return:
        """

        return self.get_simple_variable("fs.is_subdirectory('" + path + "', '" + parent_path + "')")

    # -----------------------------------------------------------------

    @property
    def installed_packages(self):

        """
        This function ...
        :return:
        """

        # Import PTS introspection tools
        self.import_package("introspection", from_name="pts.core.tools")

        # Get packages
        packages = self.get_simple_variable("introspection.installed_python_packages()")

        # Return the dictionary of python packages with their version numbers
        return packages

    # -----------------------------------------------------------------

    def is_present_package(self, package):

        """"
        This function ...
        :param package:
        """

        # Import PTS introspection tools
        self.import_package("introspection", from_name="pts.core.tools")

        # Check if present
        present = self.get_simple_variable("introspection.is_present_package('" + package + "')")

        # Return
        return present

    # -----------------------------------------------------------------

    @property
    def home_directory(self):

        """
        This function ...
        :return:
        """

        # If we are in a python session
        self.import_package("expanduser", from_name="os.path")
        return self.get_string("expanduser('~')")

    # -----------------------------------------------------------------

    @property
    def session_temp_directory(self):

        """
        This function ...
        :return:
        """

        return self.remote.session_temp_directory

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Stop the python session
        end_python_command = "screen -r " + self.screen_name + " -X stuff $'exit()\n'"
        self.remote.execute(end_python_command)

        # Stop screen session
        self.remote.kill_screen(self.screen_name)

        # Remove the pipe file
        self.remote.remove_file(self.out_pipe_filepath)

# -----------------------------------------------------------------
