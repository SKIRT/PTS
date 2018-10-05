#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.terminal Provides functions for interacting with the terminal (launching commands and getting output).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from cStringIO import StringIO
import sys
import psutil
import subprocess

# Import the relevant PTS classes and modules
from . import filesystem as fs
from . import introspection
from . import strings
from . import types

# -----------------------------------------------------------------

class Capturing(list):

    """
    This class ...
    """

    def __init__(self, pipe):

        """
        This function ...
        :param pipe:
        """

        self._pipe = pipe
        super(Capturing, self).__init__()

    # -----------------------------------------------------------------

    def __enter__(self):

        self._stdout = self._pipe
        sys.stdout = self._stringio = StringIO()
        return self

    # -----------------------------------------------------------------

    def __exit__(self, *args):

        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

# -----------------------------------------------------------------

def is_existing_executable(name):

    """
    This function ...
    :param name:
    :return:
    """

    import os

    try:
        dvnll = open(os.devnull)
        subprocess.Popen(name, stdout=dvnll, stderr=dvnll).communicate()
        return True
    except: return False

# -----------------------------------------------------------------

def executable_path(name, no_pexpect=False):

    """
    This function ...
    :param name:
    :param no_pexpect:
    :return:
    """

    if no_pexpect: output = execute_no_pexpect("which " + name)
    else: output = execute("which " + name)
    return output[0]

# -----------------------------------------------------------------

def make_executable(filepath):

    """
    This fucntion ...
    :param filepath:
    :return:
    """

    # Make executable
    subprocess.call("chmod +rx " + filepath, shell=True)

# -----------------------------------------------------------------

def run_script(filepath, options="", output=True, show_output=False, timeout=None, expect=None, no_pexpect=False):

    """
    This function ...
    :param filepath:
    :param options:
    :param output:
    :param show_output:
    :param timeout:
    :param expect:
    :param cwd:
    :return:
    """

    # Determine the path of the directory where the script is
    dir_path = fs.directory_of(filepath)
    filename = fs.name(filepath)

    make_executable(filepath)
    if no_pexpect: return execute_no_pexpect("./" + filename + " " + options, output=output, show_output=show_output, cwd=dir_path)
    else: return execute("./" + filename + " " + options, output=output, show_output=show_output, timeout=timeout, expect=expect, cwd=dir_path)

# -----------------------------------------------------------------

def execute_no_pexpect(command, output=True, show_output=False, cwd=None):

    """
    This function ...
    :param command:
    :param output:
    :param show_output:
    :param cwd:
    :return:
    """

    the_output = subprocess.check_output(command, shell=True, stderr=sys.stderr, cwd=cwd)
    if output: return the_output.split("\n")[:-1]

    #import os

    #if show_output: pipe = sys.stdout
    #else: pipe = subprocess.PIPE

    # Capture the output
    #with Capturing(pipe) as lines:

        #subprocess.call(command, shell=True, stdout=pipe, stderr=sys.stderr, cwd=cwd)
        #if output: return lines

# -----------------------------------------------------------------

def execute(command, output=True, show_output=False, timeout=None, expect=None, cwd=None, return_first=False):

    """
    This function ...
    :param output:
    :param show_output:
    :param timeout:
    :param expect:
    :param cwd:
    :param return_first:
    :return:
    """

    # Import here to accomodate fresh python installations
    import pexpect

    # Create the process
    child = pexpect.spawn(command, timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Expect
    if expect is not None: child.expect(expect)
    else: child.expect(pexpect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

    # Ignore the first and the last line (the first is the command itself, the last is always empty)
    if output:
        lines = child.before.replace('\x1b[K', '').split("\r\n")
        if return_first: return lines[:-1]
        else: return lines[1:-1]

# -----------------------------------------------------------------

def launch_fetch_lines(command, show_output=False, timeout=None, cwd=None):

    """
    This function ...
    :param command:
    :param show_output:
    :param timeout:
    :param cwd:
    :return:
    """

    # Import here to accomodate fresh python installations
    import pexpect

    # Create the process
    child = pexpect.spawn(command, timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Expect end-of-line over and over again
    end_of_line = "\r\n"
    end_of_file = pexpect.EOF
    while True:

        index = child.expect([end_of_line, end_of_file])

        if index == 0:
            lines = child.before.replace('\x1b[K', '').split("\r\n")
            if len(lines) != 1: raise ValueError("Encountered multiple lines: " + str(lines))
            yield lines[0]
        elif index == 1:
            child.logfile = None
            raise StopIteration
        else: raise RuntimeError("Something went wrong")

# -----------------------------------------------------------------

def fetch_lines(child):

    """
    This function ...
    :param child:
    :return:
    """

    import re
    import pexpect

    ansi_escape = re.compile(r'\x1b[^m]*m')

    # Expect end-of-line over and over again
    end_of_line = "\r\n"
    end_of_file = pexpect.EOF
    while True:

        index = child.expect([end_of_line, end_of_file])

        if index == 0:
            lines = child.before.replace('\x1b[K', '').split("\r\n")
            if len(lines) != 1: raise ValueError("Encountered multiple lines: " + str(lines))
            line = lines[0]
            line = ansi_escape.sub('', line).replace('\x1b[K', '').replace("\x08", "").replace("\xe2\x80\x98", "'").replace("\xe2\x80\x99", "'")
            yield line
        elif index == 1:
            child.logfile = None
            raise StopIteration
        else: raise RuntimeError("Something went wrong")

# -----------------------------------------------------------------

def launch_return(command, show_output=False, timeout=None, cwd=None):

    """
    This function ...
    :param command:
    :param show_output:
    :param timeout:
    :param cwd:
    :return:
    """

    # Import here to accomodate fresh python installations
    import pexpect

    # Create the process
    child = pexpect.spawn(command, timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None
    #print(child.logfile)
    # Return the child
    return child

# -----------------------------------------------------------------

def execute_lines_no_pexpect(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    output = kwargs.pop("output", True)
    show_output = kwargs.pop("show_output", False)
    cwd = kwargs.pop("cwd", None)

    # Remember the output lines
    lines = []

    # Check if all strings
    for line in args: assert types.is_string_type(line)

    # Execute the lines
    for line in args: lines += execute_no_pexpect(line, show_output=show_output, cwd=cwd)

    # Return output
    if output: return lines

# -----------------------------------------------------------------

def execute_lines_expect_clone(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    # Import our own copy of pexpect because we cannot interact with the process otherwise
    from . import expect

    # Get arguments
    output = kwargs.pop("output", True)
    show_output = kwargs.pop("show_output", False)
    timeout = kwargs.pop("timeout", None)
    cwd = kwargs.pop("cwd", None)

    # Execute first line
    assert types.is_string_type(args[0])
    child = expect.spawn(args[0], timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Loop over the lines
    for line in args[1:]:

        # Just a command where completion is expected
        if types.is_string_type(line):

            # Send the command
            child = child.sendline(line)
            # child.expect()
            # child.expect("$", timeout=timeout)

        # Tuple: something is expected and must be filled in
        elif isinstance(line, tuple):

            # Expect
            if len(line) == 3 and line[2]:

                # index = self.ssh.expect([self.ssh.PROMPT, line[0]]) # this is not working, why?
                index = child.expect(["$", line[0]], timeout=timeout)
                if index == 0: pass
                elif index == 1: child.sendline(line[1])
                # eof = self.ssh.prompt()
            else:
                # self.ssh.expect(line[0])
                # self.ssh.sendline(line[1])
                child.expect(line[0], timeout=timeout)
                child.sendline(line[1])

        # Invalid
        else: raise ValueError("Lines must be strings or tuples")

    # Expect
    child.expect(expect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

    # Return the output
    if output: return child.before.split("\r\n")[1:-1]

# -----------------------------------------------------------------

def execute_lines(*args, **kwargs):

    """
    This function ...
    :return:
    """

    # Import here to accomodate fresh python installations
    import pexpect

    # Get arguments
    output = kwargs.pop("output", True)
    show_output = kwargs.pop("show_output", False)
    timeout = kwargs.pop("timeout", None)
    cwd = kwargs.pop("cwd", None)

    # Execute first line
    assert types.is_string_type(args[0])
    child = pexpect.spawn(args[0], timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Loop over the lines
    for line in args[1:]:

        # Just a command where completion is expected
        if types.is_string_type(line):

            # Send the command
            child = child.sendline(line)

        # Tuple: something is expected and must be filled in
        elif types.is_tuple(line):

            # Expect
            if len(line) == 3 and line[2]:

                index = child.expect([pexpect.EOF, "$", line[0]], timeout=timeout)
                if index == 0: pass
                elif index == 0: pass
                elif index == 2: child.sendline(line[1])
            else:
                child.expect(line[0], timeout=timeout)
                child.sendline(line[1])

        # Invalid
        else: raise ValueError("Lines must be strings or tuples")

    # Expect
    child.expect(pexpect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

    # Return the output
    if output: return child.before.split("\r\n")[1:-1]

# -----------------------------------------------------------------

def remove_alias_from(alias, configuration_path):

    """
    This function ...
    :param alias:
    :param configuration_path:
    :return:
    """

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(configuration_path):

        if line.startswith("alias " + alias):
            remove_next = True
            if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
        elif remove_next and line.strip() == "": pass
        elif remove_next:
            lines.append(line)
            remove_next = False
        else: lines.append(line)

    # Write lines
    fs.write_lines(configuration_path, lines)

# -----------------------------------------------------------------

def remove_aliases(*args):

    """
    This function ...
    :param args:
    :return:
    """

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(introspection.shell_configuration_path()):

        if strings.startswith_any(line, ["alias " + alias for alias in args]):
            remove_next = True
            if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
        elif remove_next and line.strip() == "": pass
        elif remove_next:
            lines.append(line)
            remove_next = False
        else: lines.append(line)

    # Write lines
    fs.write_lines(introspection.shell_configuration_path(), lines)

    # Check if aliases don't exist anymore
    current_aliases = aliases()

    # Loop over the aliases that had to be removed
    for alias in args:
        if alias in current_aliases:
            for configuration_path in introspection.other_configuration_paths(): remove_alias_from(alias, configuration_path)

# -----------------------------------------------------------------

def remove_aliases_and_variables_with_comment(comment):

    """
    This function ...
    :param comment:
    :return:
    """

    # First make backup
    #fs.copy_file(introspection.shell_configuration_path(), fs.home, new_name="backup_profile")

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(introspection.shell_configuration_path()):

        if comment in line:
            remove_next = True
        elif remove_next:
            if line.strip() == "": remove_next = False
            else: pass
        else: lines.append(line)

    # Write lines
    fs.write_lines(introspection.shell_configuration_path(), lines)

# -----------------------------------------------------------------

def add_to_path_variable(value, comment=None, in_shell=False):

    """
    This function ...
    :param value:
    :param comment:
    :param in_shell:
    :return:
    """

    add_to_environment_variable("PATH", value, comment=comment, in_shell=in_shell)

# -----------------------------------------------------------------

def add_to_python_path_variable(value, comment=None, in_shell=False):

    """
    This function ...
    :param value:
    :param comment:
    :param in_shell:
    :return:
    """

    add_to_environment_variable("PYTHONPATH", value, comment=comment, in_shell=in_shell)

# -----------------------------------------------------------------

def add_to_environment_variable(variable_name, value, comment=None, in_shell=False):

    """
    This function ...
    :param variable_name:
    :param value:
    :param comment:
    :param in_shell:
    :return:
    """

    # Determine command
    export_command = "export " + variable_name + "=" + value + ":$PATH"

    # Define lines
    lines = []
    lines.append("")
    if comment is not None: lines.append("# " + comment)
    lines.append(export_command)
    lines.append("")

    # Add lines
    fs.append_lines(introspection.shell_configuration_path(), lines)

    # Run export path in the current shell to make variable visible
    if in_shell:
        #subprocess.call(export_command, shell=True) #execute(export_command) # not always working?
        subprocess.call("source " + introspection.shell_configuration_path(), shell=True)

# -----------------------------------------------------------------

def define_alias(name, alias_to, comment=None, in_shell=False):

    """
    This function ...
    :param name:
    :param alias_to:
    :param comment:
    :param in_shell:
    :return:
    """

    # Generate the command
    alias_command = 'alias ' + name + '="' + alias_to + '"'

    # Define lines
    lines = []
    lines.append("")
    if comment is not None: lines.append("# " + comment)
    lines.append(alias_command)
    lines.append("")

    # Add lines
    fs.append_lines(introspection.shell_configuration_path(), lines)

    # Execute in shell
    if in_shell:
        #execute(alias_command) # not always working?
        subprocess.call("source " + introspection.shell_configuration_path(), shell=True)

# -----------------------------------------------------------------

def aliases():

    """
    This function ...
    :return:
    """

    # DOESN'T WORK ANYMORE??
    # alias_dict = dict()
    #
    # output = execute_no_pexpect("alias")
    #
    # for line in output:
    #
    #     first, second = line.split("=")
    #
    #     alias = first.split("alias ")[1].strip()
    #
    #     command = second
    #     if command.startswith("'") and command.endswith("'"): command = command[1:-1]
    #     elif command.startswith('"') and command.endswith('"'): command = command[1:-1]
    #
    #     alias_dict[alias] = command
    #
    # # Return the dictionary
    # return alias_dict

    alias_dict = dict()

    for line in fs.read_lines(introspection.shell_configuration_path()):

        if not line.startswith("alias"): continue

        key = line.split("alias ")[1].split("=")[0].strip()
        value = line.split("=", 1)[1].strip()
        if value.startswith("'") and value.endswith("'"): value = value[1:-1]
        elif value.startswith('"') and value.endswith('"'): value = value[1:-1]

        # Set
        alias_dict[key] = value

    # Return the aliases
    return alias_dict

# -----------------------------------------------------------------

def resolve_alias(alias):

    """
    This function ...
    :param alias:
    :return:
    """

    return aliases()[alias]

# -----------------------------------------------------------------

def paths_in_path_variable():

    """
    This function ...
    :return:
    """

    output = execute_no_pexpect("echo $PATH")
    assert len(output) == 1
    return output[0].split(":")

# -----------------------------------------------------------------

def paths_in_python_path_variable():

    """
    This function ...
    :return:
    """

    output = execute_no_pexpect("echo $PYTHONPATH")
    assert len(output) == 1
    return output[0].split(":")

# -----------------------------------------------------------------

def remove_from_path_variable(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(introspection.shell_configuration_path()):

        if path + ":$PATH" in line and line.startswith("export PATH"):
            if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
            remove_next = True
        elif remove_next:
            if line.strip() == "": pass
            else:
                lines.append(line)
                remove_next = False
        else: lines.append(line)

    # Write lines
    fs.write_lines(introspection.shell_configuration_path(), lines)

# -----------------------------------------------------------------

def remove_from_path_variable_containing(string):

    """
    This function ...
    :param string:
    :return:
    """

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(introspection.shell_configuration_path()):

        if not line.startswith("export PATH="): continue

        what = line.split("PATH=")[1].split(":$PATH")[0]

        if string in what:
            if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
            remove_next = True
        elif remove_next:
            if line.strip() == "": pass
            else:
                lines.append(line)
                remove_next = False
        else: lines.append(line)

    # Write lines
    fs.write_lines(introspection.shell_configuration_path(), lines)

# -----------------------------------------------------------------

def remove_from_python_path_variable(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Lines to keep
    lines = []

    remove_next = False
    for line in fs.read_lines(introspection.shell_configuration_path()):

        if path + ":$PYTHONPATH" in line and line.startswith("export PYTHONPATH"):
            if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
            remove_next = True
        elif remove_next:
            if line.strip() == "": pass
            else:
                lines.append(line)
                remove_next = False
        else: lines.append(line)

    # Write lines
    fs.write_lines(introspection.shell_configuration_path(), lines)

# -----------------------------------------------------------------

def kill(proc_pid):

    """
    This function ...
    :param proc_pid:
    :return:
    """

    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()

# -----------------------------------------------------------------

def size():

    """
    Thisf unction ...
    :return:
    """

    import os
    rows, columns = os.popen('stty size', 'r').read().split()
    return int(rows), int(columns)

# -----------------------------------------------------------------

def nrows():

    """
    This function ...
    :return:
    """

    rows, columns = size()
    return rows

# -----------------------------------------------------------------

def ncolumns():

    """
    This function ...
    :return:
    """

    rows, columns = size()
    return columns

# -----------------------------------------------------------------

# NO, NO: this is all not possible! Python is running in the current shell so no other commands can be fired into this very shell
# I THINK
# def set_size(rows, columns):
#
#     """
#     Thisf unction ...
#     :param rows:
#     :param columns:
#     :return:
#     """
#
#     if introspection.is_macos():
#
#         if not types.is_integer_type(rows): raise ValueError("Number of rows must be integer")
#         if not types.is_integer_type(columns): raise ValueError("Number of columns must be integer")
#         command = "printf '\e[8;" + str(rows) + ";" + str(columns) + "t'"
#         print(command)
#         execute_no_pexpect(command)
#         import os
#         os.popen(command, 'r').read()
#
#     else:
#
#         set_nrows(rows)
#         set_ncolumns(columns)

# -----------------------------------------------------------------

# NO, NO: this is all not possible! Python is running in the current shell so no other commands can be fired into this very shell
# I THINK
# def set_nrows(rows):
#
#     """
#     This function ...
#     :param rows:
#     :return:
#     """
#
#     # Compose the command
#     if not types.is_integer_type(rows): raise ValueError("Number of rows must be integer")
#     command = "stty rows " + str(rows)
#     execute_no_pexpect(command)

# -----------------------------------------------------------------

# NO, NO: this is all not possible! Python is running in the current shell so no other commands can be fired into this very shell
# I THINK
# def set_ncolumns(columns):
#
#     """
#     Thisf unctino ...
#     :param columns:
#     :return:
#     """
#
#     # Compose the command
#     if not types.is_integer_type(columns): raise ValueError("Number of columns must be integer")
#     command = "stty columns " + str(columns)
#     execute_no_pexpect(command)

# -----------------------------------------------------------------
