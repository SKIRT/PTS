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
from abc import ABCMeta, abstractmethod, abstractproperty
import pexpect

# Import the relevant PTS classes and modules
from pts.core.tools import time
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log

# -----------------------------------------------------------------

class RemotePythonSession(object):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, output_path=None):

        """
        This function ...
        :param output_path:
        """

        # Set the output directory path
        self.output_path = output_path

        # The remote instance
        self.remote = None

    # -----------------------------------------------------------------

    def import_pts(self):

        """
        This function ...
        :return:
        """

        # Import standard PTS tools
        self.import_package("filesystem", as_name="fs", from_name="pts.core.tools")

        # Set logging level to match that of local PTS
        if log.is_debug():

            self.import_package("setup_log", from_name="pts.core.tools.logging")
            self.send_line("log = setup_log('DEBUG')")

        else: self.import_package("log", from_name="pts.core.tools.logging")

    # -----------------------------------------------------------------

    def import_package(self, name, as_name=None, from_name=None, show_output=False):

        """
        This function ...
        :param name:
        :param as_name:
        :param from_name:
        :param show_output:
        :return:
        """

        command = ""

        if from_name is not None:
            command += "from " + from_name + " "

        command += "import " + name

        if as_name is not None:
            command += " as " + as_name

        # Execute the import command
        #output = self.send_line(command, output=True, show_output=log.is_debug())
        self.send_line(command, show_output=show_output)

        # If output is given, this is normally not so good
        #if len(output) > 0:

            # Check output
            #last_line = output[-1]
            #if "cannot import" in last_line: log.warning(last_line)
            #if "ImportError" in last_line: log.warning(last_line)

            #return False

        #return True

    # -----------------------------------------------------------------

    @abstractmethod
    def send_line(self, command, show_output=False, timeout=None):

        """
        This function ...
        :param command:
        :param show_output:
        :param timeout:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def send_lines(self, lines):

        """
        This function ...
        :param lines:
        :return:
        """

        # Send the lines consecutively
        for line in lines:

            # Send the line
            self.send_line(line)

    # -----------------------------------------------------------------

    def with_statement_and_loop(self, with_statement, for_statement, in_loop_lines):

        """
        This function ...
        :param with_statement:
        :param for_statement:
        :param in_loop_lines:
        :param show_output:
        :return:
        """

        # With statement
        self.send_line(with_statement)

        # For loop
        self.send_line("    " + for_statement)

        # In-loop lines
        for line in in_loop_lines:
            # Send the line
            self.send_line("        " + line)

            # Show output
            # if show_output:
            #    for o in output: print(o)

        # Finish the loop
        self.send_line("    ")

    # -----------------------------------------------------------------

    def do_loop(self, top_statement, in_loop_lines):

        """
        This function ...
        :param top_statement:
        :param in_loop_lines:
        :return:
        """

        # Top statement
        self.send_line(top_statement)

        # In-loop lines
        for line in in_loop_lines:
            # Send the line
            # output = self.send_line("    " + line, output=True)

            self.send_line("    " + line)

            # Show output
            # if show_output:
            #    for o in output: print(o)

        # Finish the loop
        self.send_line("")

    # -----------------------------------------------------------------

    @abstractmethod
    def get_simple_variable(self, name, show_output=False):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_simple_property(self, variable, name, show_output=False):

        """
        This function ...
        :param variable:
        :param name:
        :param show_output:
        :return:
        """

        try: return self.get_simple_variable(variable + "." + name, show_output=show_output)
        except NameError: raise AttributeError("Variable '" + variable + "' has no attribute '" + name + "'")

    # -----------------------------------------------------------------

    def get_string(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simple_variable(name)[1:-1]

        #return output[0][1:-1]

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

    def files_in_path(self, path, recursive=False, extension=None, returns="path"):

        """
        This function ...
        :param path:
        :param recursive:
        :param extension:
        :param returns:
        :return:
        """

        return self.get_simple_variable("fs.files_in_path('" + path + "', recursive=" + repr(recursive) + ", extension=" + repr(extension) + ", returns=" + repr(returns) + ")")

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

    def read_lines(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        for line in self.get_simple_variable("list(fs.read_lines('" + path + "'))"): yield line

    # -----------------------------------------------------------------

    def read_lines_reversed(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        for line in self.get_simple_variable("list(fs.read_lines_reversed('" + path + "'))"): yield line

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

class AttachedPythonSession(RemotePythonSession):

    """
    This function ...
    """

    def __init__(self, remote, python_command="python", assume_pts=False, output_path=None):

        """
        This function ...
        :param remote:
        :param assume_pts:
        :param output_path:
        """

        # Call the constructor of the base class
        super(AttachedPythonSession, self).__init__(output_path)

        from .remote import Remote

        # The remote connection
        if isinstance(remote, basestring):

            # Set the original remote
            self.remote = Remote()
            if not self.remote.setup(remote): raise RuntimeError("Could not connect to the remote")

            # Set the parent remote
            self.parent = Remote()
            if not self.parent.setup(remote): raise RuntimeError("Could not create a seperate connection to the remote for the attached python session")

        else:

            # Set the original remote
            self.remote = remote

            # Create a new remote
            self.parent = Remote()
            connected = self.parent.setup(self.remote.host_id)
            if not connected: raise RuntimeError("Could not create a seperate connection to the remote for the attached python session")

        # Start python
        self.parent.execute(python_command, expect=">>>")

        # Import PTS stuff
        if assume_pts: self.import_pts()

    # -----------------------------------------------------------------

    def send_line(self, line, show_output=False, timeout=None):

        """
        This function ...
        :param line:
        :param show_output:
        :param timeout:
        :return:
        """

        # Execute the line
        self.execute(line, timeout=timeout, show_output=show_output)

    # -----------------------------------------------------------------

    def execute(self, line, timeout=None, show_output=False):

        """
        This function ...
        :param line:
        :param timeout:
        :param show_output:
        :return:
        """

        # Debugging
        #log.debug("The command to execute the line is:")
        #log.debug(send_command)

        # Send the line
        self.parent.execute(line, show_output=show_output, timeout=timeout, expect=">>>")

        # Sleep for a while so that we are sure that the actual python stuff has reached the interactive python session within the screen
        #time.wait(5) # in seconds

    # -----------------------------------------------------------------

    def get_simple_variable(self, name, show_output=False):

        """
        This function ...
        :param name:
        :param show_output:
        :return:
        """

        self.send_line(name, show_output=show_output)

        # Return the value
        return eval(self.parent.ssh.before)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Stop the python session
        end_python_command = "exit()"
        self.parent.execute(end_python_command)

# -----------------------------------------------------------------

class DetachedPythonSession(RemotePythonSession):

    """
    This class ...
    """

    def __init__(self, remote, python_command="python", assume_pts=False, tmux=False, output_path=None):

        """
        This function ...
        :param remote:
        :param assume_pts:
        :param tmux: use tmux instead of screen
        :return:
        """

        # Set tmux flag
        self.tmux = tmux

        # Call the constructor of the base class
        super(DetachedPythonSession, self).__init__(output_path)

        # The remote connection
        if isinstance(remote, basestring):
            from .remote import Remote
            self.remote = Remote()
            self.remote.setup(remote)
        else: self.remote = remote

        # The path to the out pipe file
        self.out_pipe_filepath = None

        # The last number of lines in the session output
        self.previous_length = 0

        # Generate session ID and screen name
        self.session_id = time.unique_name()
        self.screen_name = "pts_remotepython_" + self.session_id

        # Flag that indicates whether the session is attached
        self.attached = False

        # Create the pipe
        self.create_pipe()

        # Start the session
        self.start_session()

        # Set the output pipe
        self.set_out_pipe()

        # Import PTS stuff
        if assume_pts: self.import_pts()

    # -----------------------------------------------------------------

    def create_pipe(self):

        """
        This function ...
        :return:
        """

        # Create pipe file
        out_pipe_filename = "out_pipe_" + self.session_id + ".txt"
        self.out_pipe_filepath = fs.join(self.remote.pts_temp_path, out_pipe_filename)

        # Debugging
        log.debug("Creating pipe file '" + self.out_pipe_filepath + "' on remote ...")

        # Create the pipe file for output
        if self.remote.is_file(self.out_pipe_filepath): self.remote.remove_file(self.out_pipe_filepath)
        #self.remote.touch(self.out_pipe_filepath)
        #self.remote.write_line(self.out_pipe_filepath, "")
        self.remote.touch_alternative(self.out_pipe_filepath)

    # -----------------------------------------------------------------

    @property
    def screenlog_path(self):

        """
        This function ...
        :return:
        """

        if self.output_path is None: return None

        the_index = 0
        the_path = None

        # Find screenlog file
        for path in self.remote.files_in_path(self.output_path):

            name = fs.name(path)
            if name.startswith("screenlog"):

                index = int(name.split(".")[1])
                if the_path is None:
                    the_path = path
                    the_index = index
                elif index > the_index:
                    the_path = path
                    the_index = index

        # Return the path
        return the_path

    # -----------------------------------------------------------------

    def start_session(self, python_command="python"):

        """
        This function ...
        :param python_command:
        :return:
        """

        # Inform the user
        log.info("Starting session for python ...")

        # Start the screen
        #start_screen_command = "screen -dmS " + self.screen_name

        # Debugging
        #log.debug("Starting screen session with command:")
        #log.debug(start_screen_command)

        # Start the screen
        #self.remote.execute(start_screen_command, show_output=True)

        # Determine command
        #command = "tmux new-session -d -n " + self.screen_name + " python > " + self.out_pipe_filepath

        if self.tmux: command = "tmux new -d -n " + self.screen_name + " " + python_command
        elif self.output_path is not None: command = "screen -dmS " + self.screen_name + " -L"
        else: command = "screen -dmS " + self.screen_name

        #command = "tmux new -d -n " + self.screen_name + " python > " + self.out_pipe_filepath

        #command = "tmux new -d -n " + self.screen_name + " python \; pipe-pane 'cat > " + self.out_pipe_filepath + "'"

        # Debugging
        log.debug("Starting session with the command:")
        log.debug(command)

        # Execute the command
        self.remote.execute(command, show_output=True, cwd=self.output_path)

        # Tmux ls
        #self.remote.execute("tmux ls", show_output=True)

        #to_out_pipe_command = "tmux pipe-pane -o -t " + self.screen_name + " 'cat >> " + self.out_pipe_filepath + "'"

        # Debugging
        #log.debug("Piping output with the command:")
        #log.debug(to_out_pipe_command)

        #self.remote.execute(to_out_pipe_command, show_output=True)

        if not self.tmux:

            start_python_command = "screen -S " + self.screen_name + " -p 0 -X stuff '" + python_command + "\n'"
            self.remote.execute(start_python_command)

    # -----------------------------------------------------------------

    def set_out_pipe(self):

        """
        This function ...
        :return:
        """

        # Start python in the screen session
        # start_python_command = "screen -r " + self.screen_name + " -X stuff $'python > " + self.out_pipe_filepath + " < " + self.in_pipe_filepath + "\n'"
        #start_python_command = "screen -r -S " + self.screen_name + " -X stuff $'python > " + self.out_pipe_filepath + "\n'"

        #start_python_command = "screen -r -S " + self.screen_name + " -X stuff $'python\n'"

        #start_python_command = "screen -r -X stuff $'python\n'"
        #start_python_command = "screen -r -S " + self.screen_name + " -p 0 -X stuff 'python > " + self.out_pipe_filepath + "\n'"
        #start_python_command = 'screen -r ' + self.screen_name + ' -X stuff "python > ' + self.out_pipe_filepath + '"'

        #start_python_command = "screen -r -S " + self.screen_name

        # Debugging
        #log.debug("Starting python session with command:")
        #log.debug(start_python_command)

        # Start python
        #self.remote.execute(start_python_command, expect="$")
        #self.remote.execute(start_python_command)

        #start_python_command = "python > " + self.out_pipe_filepath

        # Start python again
        #self.remote.execute(start_python_command, expect=">>>")

        # Detach
        #self.remote.execute("$'\x01'")
        #self.remote.execute("d")

        #self.remote.execute("\n")

        self.send_line("import sys")
        #self.send_line(r"newline = '\\\n'")
        self.send_line("import os")
        self.send_line("newline = os.linesep")
        #if self.tmux: self.do_loop("def my_display(x):", ["with open('" + self.out_pipe_filepath + r"', 'a') as fh: fh.write(str(x) + '\n')"])
        #else: self.do_loop("def my_display(x):", ["with open('" + self.out_pipe_filepath + "', 'a') as fh: fh.write(str(x) + newline)"])

        # Define function
        self.send_line("def my_display(x):")
        self.send_line("    with open('" + self.out_pipe_filepath + "', 'a') as fh:")
        self.send_line("        fh.write(str(x) + newline)")
        self.send_line("        fh.flush()")
        self.send_line("")

        self.send_line("sys.displayhook = my_display")

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

    def attach(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Attaching to the screen session " + self.screen_name + " ...")

        # Tmux or screen
        if self.tmux: self.remote.execute("tmux a -t " + self.screen_name, expect=">>>")
        else:

            #self.remote.execute("screen -r " + self.screen_name, expect=">>>")
            self.remote.ssh.sendline("screen -r " + self.screen_name)

            while True:

                index = self.remote.ssh.expect([">>>", pexpect.TIMEOUT])
                if index == 1: break

        # Set attached flag
        self.attached = True

        # Set displayhook to default
        #self.send_line("sys.displayhook = sys.__displayhook__")
        self.remote.execute("sys.displayhook = sys.__displayhook__", expect=">>>")

    # -----------------------------------------------------------------

    def detach(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Detaching from the screen session " + self.screen_name + " ...")

        # Check
        if not self.attached: log.warning("Not attached")

        # Set displayhook back to custom function
        #self.send_line("sys.displayhook = my_display")
        self.remote.execute("sys.displayhook = my_display", expect=">>>")

        # Tmux
        if self.tmux:

            # Send Ctrl+b d to detach
            self.remote.execute("^B")
            self.remote.execute("d")

            # Match the prompt
            self.remote.ssh.prompt()

        # Screen
        else:

            # Send Ctrl+A d to detach
            self.remote.ssh.send("^A")
            self.remote.ssh.send("d")

            # Match the prompt
            self.remote.ssh.prompt()

    # -----------------------------------------------------------------

    def get_simple_variable(self, name, show_output=False, max_attempts=4):

        """
        This function ...
        :param name:
        :param show_output:
        :param max_attempts:
        :return:
        """

        if show_output:

            output = self.send_line(name)
            print(output)

        else:

            # Marker
            self.send_line("'[PTS]'")

            # Spit out the value of the variable with name 'name'
            self.send_line(name)

            # Wait for some time
            time.wait(10)

            # Get output
            lines = []
            for line in self.remote.read_lines_reversed(self.out_pipe_filepath):
                #print("LINE", line)
                if line == "[PTS]": break
                lines.append(line)
            lines.reverse()

            #print("1", lines)

            attempts = 1
            while (len(lines) == 0 or lines[0] == "") and attempts <= max_attempts:

                time.wait(15)

                # Get output
                lines = []
                for line in self.remote.read_lines_reversed(self.out_pipe_filepath):
                    #print("LINE", line)
                    if line == "[PTS]": break
                    lines.append(line)
                lines.reverse()

                # Increment the number of attempts
                attempts += 1

            #print("2", lines)

            # LOOK FOR AN ERROR
            if len(lines) == 0 or lines[0] == "":

                # LOOK FOR ERRORS
                #print("NAME", name)
                after = []
                # Get output from screen output file
                screenlog_lines_reversed = list(self.remote.read_lines_reversed(self.screenlog_path))
                for index in range(len(screenlog_lines_reversed)):
                    #print("LINE", line)
                    #print(line, ">>> " + name)
                    line = screenlog_lines_reversed[index]
                    #print(line)
                    if ">>> " + name in line:
                        after = screenlog_lines_reversed[1:index]
                        break

                after.reverse()

                for line in after:
                    if "Traceback (most recent call last)" in line: raise RuntimeError(after[-1])

                else: raise RuntimeError("Could not determine the value of '" + name)

                #lines = after

            # Return the value
            return eval(lines[0])

    # -----------------------------------------------------------------

    def send_line(self, line, show_output=False, timeout=None):

        """
        This function ...
        :param line:
        :param show_output:
        :param timeout:
        :return:
        """

        # Attached or detached
        if show_output: return self.execute_attached(line, timeout=timeout)
        else: return self.execute_detached(line, timeout=timeout)

    # -----------------------------------------------------------------

    def execute_attached(self, line, timeout=None):

        """
        This function ...
        :param line:
        :param timeout:
        :return:
        """

        # First attach
        self.attach()

        # Send the line, show output since that is the reason we executed in attached mode
        output = self.remote.execute(line, timeout=timeout, show_output=True, expect=">>>")

        # Detach again
        self.detach()

        # Return the output
        return output

    # -----------------------------------------------------------------

    def execute_detached(self, line, timeout=None):

        """
        This function ...
        :param line:
        :param timeout:
        :return:
        """

        if self.tmux: send_command = 'tmux send-keys -t ' + self.screen_name + ' "' + line + '" Enter'
        else:

            if "'" not in line: send_command = "screen -S " + self.screen_name + " -p 0 -X stuff '" + line + "\n'"
            elif '"' not in line: send_command = 'screen -S ' + self.screen_name + ' -p 0 -X stuff "' + line + '\n"'
            else: raise ValueError("Line cannot contain both single quotes and double quotes")

        # Debugging
        log.debug("The command to execute the line is:")
        log.debug(send_command)

        # Send the line
        self.remote.execute(send_command, show_output=log.is_debug(), timeout=timeout)

        # Sleep for a while so that we are sure that the actual python stuff has reached the interactive python session within the screen
        time.wait(5) # in seconds

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

        # Debugging
        if self.output_path is not None: log.debug("Screen output has been placed in '" + self.output_path + "'")

        # Detach if we are still attached
        if self.attached: self.detach()

        # Remove things if we are not in debug mode
        if not log.is_debug():

            # Using tmux
            if self.tmux:

                end_session_command = "tmux kill-session -t " + self.screen_name
                self.remote.execute(end_session_command)

            # Using screen
            else:

                # Stop the python session
                end_python_command = "screen -r -S " + self.screen_name + " -X stuff 'exit()\n'"
                self.remote.execute(end_python_command)

                # Stop screen session
                self.remote.kill_screen(self.screen_name)

            # Remove the pipe file
            self.remote.remove_file(self.out_pipe_filepath)

        else:

            # Debugging info
            log.debug("Not closing session and removing pipe in debug mode")
            if self.tmux: log.debug(" - Tmux session name: " + self.screen_name)
            else: log.debug(" - Screen session name: " + self.screen_name)
            log.debug(" - pipe file path: " + self.out_pipe_filepath)
            if self.output_path is not None: log.debug(" - screen log file: " + self.screenlog_path)

# -----------------------------------------------------------------
