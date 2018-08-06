#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.remote Contains the Remote class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import sys
import pexpect
from pexpect import pxssh, ExceptionPexpect
import tempfile
import StringIO
import subprocess
from lxml import etree

# Import the relevant PTS classes and modules
from .host import Host, load_host
from .utils import HostDownException
from .vpn import VPN
from ..basics.log import log
from ..tools import parsing
from ..tools import filesystem as fs
from ..tools import time
from ..basics.task import Task
from ..tools import introspection
from ..tools.introspection import possible_cpp_compilers, possible_mpi_compilers, possible_mpirun_names
from .python import AttachedPythonSession, DetachedPythonSession
from ..units.parsing import parse_unit as u
from ..basics.map import Map
from ..tools import strings, types
from ..tools.utils import lazyproperty
from ..tools.utils import memoize_method

# -----------------------------------------------------------------

class TimeOutReached(RuntimeError):

    """
    This class ...
    """

    def __init__(self, message, timeout, command=None):

        """
        Thisf unction ...
        :param message:
        :param timeout: the max timeout that was used
        :param command:
        """

        # Call the base class constructor with the parameters it needs
        super(TimeOutReached, self).__init__(message)

        # The timeout and command
        self.timeout = timeout
        self.command = command

# -----------------------------------------------------------------

def is_available(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    remote = Remote()
    success = remote.setup(host_id)
    del remote
    return success

# -----------------------------------------------------------------

def active_keys():

    """
    This function ...
    :return:
    """

    try:
        output = subprocess.check_output(["ssh-add", "-L"])
    except subprocess.CalledProcessError as e:
        #log.error(e.output)
        #exit()
        return []

    names = []
    for line in output.split("\n"):
        if not line: continue
        _, key, path = line.split(" ")
        name = fs.name(path)
        names.append(name)

    return names

# -----------------------------------------------------------------

def add_key(name, password=None, show_output=False):

    """
    This function ...
    :param name:
    :param password:
    :param show_output:
    :return:
    """

    # Check if the key exists
    key_path = fs.absolute_path(fs.join("~/.ssh", name))
    if not fs.is_file(key_path): raise ValueError("The key '" + name + "' does not exist in the .ssh directory")

    command = ["ssh-add", key_path]
    command_string =  " ".join(command)

    # Password can be asked or not
    child = pexpect.spawn(command_string)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    index = child.expect([pexpect.EOF, name + ":"])

    if index == 1:

        if password is None: raise RuntimeError("Password is asked but none is passed to this function")

        # Send the password
        child.sendline(password)

        # Expect EOF
        child.expect(pexpect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

# -----------------------------------------------------------------

def load_remote(remote, silent=False):

    """
    This function ...
    :param remote:
    :param silent:
    :return:
    """

    # Make remote
    if types.is_string_type(remote): remote = Remote(host_id=remote, silent=silent)
    elif isinstance(remote, Host): remote = Remote(host_id=remote, silent=silent)
    elif not isinstance(remote, Remote): raise ValueError("Remote must be string, Host or Remote object")

    # Return the remote instance
    return remote

# -----------------------------------------------------------------

def get_host_id(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    # Get host ID
    if isinstance(host_id, Host): the_host_id = host_id.id
    elif types.is_string_type(host_id): the_host_id = host_id
    elif isinstance(host_id, Remote): the_host_id = host_id.host_id
    else: raise ValueError("Invalid value for 'host_id'")

    # Return the host iD
    return the_host_id

# -----------------------------------------------------------------

def get_home_path(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    return Remote(host_id=host_id).home_directory

# -----------------------------------------------------------------

class Remote(object):

    """
    This function ...
    """

    def __init__(self, host_id=None, silent=False, log_conda=False):

        """
        The constructor ...
        :param host_id:
        :param silent:
        :param log_conda:
        :return:
        """

        # Call the constructor of the base class
        super(Remote, self).__init__()

        # -- Attributes --

        # The SSH interface, an instance of the pxssh class
        self.ssh = pxssh.pxssh()

        # The host instance
        self.host = None

        # The VPN service
        self.vpn = None

        # A flag indicating whether the connection with the remote has been established
        self.connected = False

        # A regular expression object that strips away special unicode characters, used on the remote console output
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

        # Flags
        self.log_conda = log_conda
        self.conda_activated = False

        # Remember the commands that were executed on the remote host
        self.commands = []

        # Set silent flag
        self.silent = silent

        # If host ID is given, setup
        if host_id is not None:
            if not self.setup(host_id, silent=self.silent): log.warning("The connection could not be made. Run setup().")

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster_name=None, login_timeout=30, nrows=None, ncols=200, one_attempt=False,
              retry_factor=4, silent=False):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :param login_timeout:
        :param nrows:
        :param ncols:
        :param one_attempt:
        :param retry_factor:
        :param silent:
        :return:
        """

        # Create the host object
        if isinstance(host_id, Host):
            self.host = host_id
            if cluster_name is not None: raise ValueError("Cannot specify cluster name if host is given")
        elif types.is_string_type(host_id): self.host = load_host(host_id, cluster_name)
        else: raise ValueError("Invalid value for 'host_id'")

        # Set the host ID
        host_id = self.host.id

        # If a VPN connection is required for the remote host
        if self.host.requires_vpn: self.connect_to_vpn()

        # Check if key is active
        if self.host.key is not None:
            if self.host.key not in active_keys(): add_key(self.host.key, self.host.key_password)

        # Set the silent flag
        self.silent = silent

        # Make the connection
        try: self.login(login_timeout, silent=silent)
        except HostDownException:

            if one_attempt:
                self.warning("Could not connect to the remote host")
                self.ssh = pxssh.pxssh()
                return False

            # Warning
            log.warning("Connection to host '" + host_id + "' failed, trying again ...")
            self.ssh = pxssh.pxssh()
            try: self.login(login_timeout * retry_factor, silent=silent) # try now with a timeout that is four times as long
            except HostDownException:
                self.warning("Could not connect to the remote host")
                self.ssh = pxssh.pxssh()
                return False

        # Swap cluster
        if self.host.cluster_name is not None: self.swap_cluster(self.host.cluster_name)

        # LMOD_DISABLE_SAME_NAME_AUTOSWAP
        #self.define_environment_variable("LMOD_DISABLE_SAME_NAME_AUTOSWAP", "yes")

        # Set the nrows and ncols
        if nrows is not None: self.nrows = nrows
        if ncols is not None: self.ncols = ncols

        # Return whether the connection was made
        return self.connected

    # -----------------------------------------------------------------

    def __del__(self):

        """
        The destructor ...
        :return:
        """

        # Clear the temporary data
        #self.clear_pts_temp()

        # Disconnect from the remote host
        if self.connected: self.logout()

    # -----------------------------------------------------------------

    @property
    def ncols(self):

        """
        This function ...
        :return:
        """

        command = "tput cols"
        return int(self.execute(command)[0])

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        command = "tput lines"
        return int(self.execute(command)[0])

    # -----------------------------------------------------------------

    @nrows.setter
    def nrows(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Compose the command
        assert int(value) == value
        command = "stty rows " + str(value)
        self.execute(command, output=False)

    # -----------------------------------------------------------------

    @ncols.setter
    def ncols(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Compose the command
        assert int(value) == value
        command = "stty columns " + str(value)
        self.execute(command, output=False)

    # -----------------------------------------------------------------

    @property
    def size(self):

        """
        This function ...
        :return:
        """

        return self.nrows, self.ncols

    # -----------------------------------------------------------------

    @size.setter
    def size(self, size_tuple):

        """
        This function ...
        :param size_tuple:
        :return:
        """

        self.nrows = size_tuple[0]
        self.ncols = size_tuple[1]

    # -----------------------------------------------------------------

    def info(self, message):

        """
        This function ...
        :param message:
        :return:
        """

        # Generate prefix
        prefix = self.host_id
        if self.log_conda and self.connected:
            name = self.conda_active_environment_fast(assert_activated=self.conda_activated)
            if name is not None: prefix += ", " + name

        # Inform the user
        log.info("[" + prefix + "] " + message)

    # -----------------------------------------------------------------

    def debug(self, message):

        """
        This function ...
        :param message:
        :return:
        """

        # Generate prefix
        prefix = self.host_id
        if self.log_conda:
            name = self.conda_active_environment_fast(assert_activated=self.conda_activated)
            if name is not None: prefix += ", " + name

        # Debugging
        log.debug("[" + prefix + "] " + message)

    # -----------------------------------------------------------------

    def warning(self, message):

        """
        This function ...
        :param message:
        :return:
        """

        # Generate prefix
        prefix = self.host_id
        if self.log_conda:
            name = self.conda_active_environment_fast(assert_activated=self.conda_activated)
            if name is not None: prefix += ", " + name

        # Give warning
        log.warning("[" + prefix + "] " + message)

    # -----------------------------------------------------------------

    def error(self, message):

        """
        This function ...
        :param message:
        :return:
        """

        # Generate prefix
        prefix = self.host_id
        if self.log_conda:
            name = self.conda_active_environment_fast(assert_activated=self.conda_activated)
            if name is not None: prefix += ", " + name

        # Give error message
        log.error("[" + prefix + "] " + message)

    # -----------------------------------------------------------------

    def success(self, message):

        """
        This function ...
        :param message:
        :return:
        """

        # Generate prefix
        prefix = self.host_id
        if self.log_conda:
            name = self.conda_active_environment_fast(assert_activated=self.conda_activated)
            if name is not None: prefix += ", " + name

        # Show success message
        log.success("[" + prefix + "] " + message)

    # -----------------------------------------------------------------

    def define_environment_variable(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        command = "export " + name + "=" + value
        self.execute(command)

    # -----------------------------------------------------------------

    def purge(self):

        """
        This function ...
        :return:
        """

        self.execute("module purge")

    # -----------------------------------------------------------------

    def load_modules(self, *args):

        """
        This function ...
        :return:
        """

        for name in args: self.load_module(name)

    # -----------------------------------------------------------------

    def load_module(self, module_name, show_output=False, return_success=False):

        """
        This function ...
        :param module_name:
        :param show_output:
        :param return_success:
        :return:
        """

        output = self.execute("module load " + module_name, show_output=show_output)

        # Look for errors
        for line in output:

            if "Lmod has detected the following error" in line:

                # Create error message
                error_message = line.split("error:")[1].strip()
                if "see output of 'ml'" in error_message:
                    error_message.replace("(see output of 'ml')", "")
                    error_message += " Loaded modules: " + ", ".join(self.loaded_modules)

                if return_success:
                    log.warning(error_message)
                    return False
                else: raise RuntimeError(error_message)

        # Success: no error encountered
        return True

    # -----------------------------------------------------------------

    def unload_module(self, module_name, show_output=False):

        """
        This function ...
        :param module_name:
        :param show_output:
        :return:
        """

        self.execute("module del " + module_name, show_output=show_output)

    # -----------------------------------------------------------------

    def reload_all_modules(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        self.execute("module update", show_output=show_output)

    # -----------------------------------------------------------------

    def unload_all_modules(self):

        """
        This function ...
        :return:
        """

        # This remote does not use module system lmod
        if not self.has_lmod: return

        # Inform the user
        self.info("Unloading all modules ...")

        # Module purge
        self.execute("module purge", output=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def available_modules(self):

        """
        This function ...
        :return:
        """

        output = self.execute("module spider")

        modules = dict()

        name = None
        description = None
        was_empty = True

        for line in output:

            if not line:

                modules[name] = description
                was_empty = True
                name = None
                description = None
                continue

            if not line.startswith("  "): continue

            if was_empty:
                name = line.split("  ")[1].split(":")[0]
                was_empty = False
            elif description is None:
                description = line.split("    ")[1]
                was_empty = False
            else:
                #print(line)
                description += line.split("    ")[1]
                was_empty = False

        return modules

    # -----------------------------------------------------------------

    def get_module_versions(self, module_name, exact=True, show_output=False):

        """
        This function ...
        :param module_name:
        :param exact:
        :param show_output:
        :return:
        """

        if exact: command = "module spider " + module_name
        else: command = "module spider -r " + module_name

        #print(command)
        #print(self.ssh.before)
        output = self.execute(command, show_output=show_output)

        versions = []

        triggered = False
        for line in output:

            if triggered:

                if not line.strip(): break
                if "Other possible modules matches" in line: break

                version = line.strip()
                versions.append(version)

            elif "Versions:" in line: triggered = True

            elif module_name + ":" in line:

                if line.split(":")[1].strip() == "": continue
                else:

                    # only one version
                    return [line.split(":")[1].strip()]

        return versions

    # -----------------------------------------------------------------

    @property
    def has_lmod(self):

        """
        This function ...
        :return:
        """

        output = self.execute("module -v")

        if "command not found" in output[0]: return False

        for line in output:
            if "Modules based on Lua" in line: return True

        return False # other executable called 'module'

    # -----------------------------------------------------------------

    @property
    def loaded_modules(self):

        """
        This function ...
        :return:
        """

        #output = self.execute("module list")
        output = self.execute("ml")

        modules = []

        for line in output:

            if not line: continue

            if not line.startswith("  "): continue

            line = line[2:]

            if ") " not in line: continue

            splitted = line.split("  ")

            for part in splitted:

                if not part.strip(): continue
                if ") " not in part: continue
                name = part.split(") ")[1].strip()
                if not name: continue
                modules.append(name)

        return modules

    # -----------------------------------------------------------------

    @property
    def bashrc_path(self):

        """
        This function ...
        :return:
        """

        # Check if bashrc exists, if not create it if we are on Linux
        bashrc_path = fs.join(self.home_directory, ".bashrc")
        if self.is_linux and not self.is_file(bashrc_path): self.touch(bashrc_path)

        # Return the path
        return bashrc_path

    # -----------------------------------------------------------------

    @property
    def profile_path(self):

        """
        This function ...
        :return:
        """

        # Check whether .profile exists, create it if we are on MacOS
        profile_path = fs.join(self.home_directory, ".profile")
        if self.is_macos and not self.is_file(profile_path): self.touch(profile_path)

        # Return the path
        return profile_path

    # -----------------------------------------------------------------

    @property
    def bash_profile_path(self):

        """
        This function ...
        :return:
        """

        bash_profile_path = fs.join(self.home_directory, ".bash_profile")
        return bash_profile_path

    # -----------------------------------------------------------------

    @property
    def shell_configuration_path(self):

        """
        This function ...
        :return:
        """

        if self.is_macos: return self.profile_path
        elif self.is_linux: return self.bashrc_path
        else: raise NotImplemented("System must be running MacOS or Linux")

    # -----------------------------------------------------------------

    @property
    def other_configuration_paths(self):

        """
        This function ...
        :return:
        """

        paths = []

        if self.is_macos:
            if self.is_file(self.bashrc_path): paths.append(self.bashrc_path)
            if self.is_file(self.bash_profile_path): paths.append(self.bash_profile_path)
        elif self.is_linux:
            if self.is_file(self.bash_profile_path): paths.append(self.bash_profile_path)
            if self.is_file(self.profile_path): paths.append(self.profile_path)
        else: raise NotImplementedError("System must be running MacOS or Linux")

        return paths

    # -----------------------------------------------------------------

    def fix_configuration_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.info("Checking and fixing shell configuration files ...")

        # On MacOS, don't need fixing
        if self.is_macos: return

        # If bash_profile file exists
        if self.is_file(self.bash_profile_path):

            # Check if it points to the bashrc file
            links_bashrc = ".bashrc" in ";".join(self.read_lines(self.bash_profile_path))

            # Add link to bashrc
            if not links_bashrc: self._add_bashrc_link(self.bash_profile_path)

        # If profile file exists
        elif self.is_file(self.profile_path):

            # Check if it points to the bashrc file
            links_bashrc = ".bashrc" in ";".join(self.read_lines(self.profile_path))

            # Add link to bashrc
            if not links_bashrc: self._add_bashrc_link(self.profile_path)

        # Profile and bash_profile files do not exist
        else:

            # Make bash profile
            self.touch(self.bash_profile_path)

            # Make it link to bashrc
            self._add_bashrc_link(self.bash_profile_path)

    # -----------------------------------------------------------------

    def _add_bashrc_link(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        lines = []
        lines.append("")
        lines.append("# When logging in from a console, .bashrc will be called.")
        lines.append("if [ -f ~/.bashrc ]; then")
        lines.append("   source ~/.bashrc")
        lines.append("fi")
        lines.append("")

        # Add the lines
        self.append_lines(path, lines)

    # -----------------------------------------------------------------

    @property
    def paths_in_path_variable(self):

        """
        This function ...
        :return:
        """

        output = self.execute("echo $PATH")
        assert len(output) == 1
        return output[0].split(":")

    # -----------------------------------------------------------------

    @property
    def paths_in_python_path_variable(self):

        """
        This function ...
        :return:
        """

        output = self.execute("echo $PYTHONPATH")
        assert len(output) == 1
        return output[0].split(":")

    # -----------------------------------------------------------------

    def add_to_path_variable(self, value, comment=None, in_shell=False):

        """
        This function ...
        :param value:
        :param comment:
        :param in_shell:
        :return:
        """

        # add to PATH
        self.add_to_environment_variable("PATH", value, comment=comment, in_shell=in_shell)

    # -----------------------------------------------------------------

    def add_to_python_path_variable(self, value, comment=None, in_shell=False):

        """
        This function ...
        :param value:
        :param comment:
        :param in_shell:
        :return:
        """

        self.add_to_environment_variable("PYTHONPATH", value, comment=comment, in_shell=in_shell)

    # -----------------------------------------------------------------

    def add_to_environment_variable(self, variable_name, value, comment=None, in_shell=False):

        """
        This function ...
        :param variable_name:
        :param value:
        :param comment:
        :param in_shell:
        :return:
        """

        # Set PYTHONPATH
        export_command = "export " + variable_name + "=" + value + ":$" + variable_name

        # Define lines
        lines = []
        lines.append("")
        if comment is not None: lines.append("# " + comment)
        lines.append(export_command)
        lines.append("")

        # Add lines
        self.append_lines(self.shell_configuration_path, lines)

        # Execute in shell
        if in_shell: self.execute(export_command)

    # -----------------------------------------------------------------

    def define_alias(self, name, alias_to, comment=None, in_shell=False):

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
        self.append_lines(self.shell_configuration_path, lines)

        # Execute in shell
        if in_shell: self.execute(alias_command)

    # -----------------------------------------------------------------

    def remove_aliases_and_variables_with_comment(self, comment):

        """
        This function ...
        :param comment:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(self.shell_configuration_path):

            if comment in line: remove_next = True
            elif remove_next:
                if line.strip() == "": remove_next = False
                else: pass
            else: lines.append(line)

        # Write lines
        self.write_lines(self.shell_configuration_path, lines)

    # -----------------------------------------------------------------

    def remove_alias_from(self, alias, configuration_path):

        """
        This function ...
        :param alias:
        :param configuration_path:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(configuration_path):

            if line.startswith("alias " + alias):
                remove_next = True
                if len(lines) > 0 and lines[-1].startswith("#"): del lines[-1]
            elif remove_next and line.strip() == "":
                pass
            elif remove_next:
                lines.append(line)
                remove_next = False
            else: lines.append(line)

        # Write lines
        self.write_lines(configuration_path, lines)

    # -----------------------------------------------------------------

    def remove_aliases(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(self.shell_configuration_path):

            if strings.startswith_any(line, ["alias " + alias for alias in args]): remove_next = True
            elif remove_next and line.strip() == "": pass
            elif remove_next:
                lines.append(line)
                remove_next = False
            else: lines.append(line)

        # Write lines
        self.write_lines(self.shell_configuration_path, lines)

        # Check if aliases don't exist anymore
        current_aliases = self.aliases

        # Loop over the aliases that had to be removed
        for alias in args:
            if alias in current_aliases:
                for configuration_path in self.other_configuration_paths: self.remove_alias_from(alias, configuration_path)

    # -----------------------------------------------------------------

    def remove_from_path_variable(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(self.shell_configuration_path):

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
        self.write_lines(self.shell_configuration_path, lines)

    # -----------------------------------------------------------------

    def remove_from_path_variable_containing(self, string):

        """
        This function ...
        :param string:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(self.shell_configuration_path):

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
        self.write_lines(self.shell_configuration_path, lines)

    # -----------------------------------------------------------------

    def remove_from_python_path_variable(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Lines to keep
        lines = []

        remove_next = False
        for line in self.read_lines(self.shell_configuration_path):

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
        self.write_lines(self.shell_configuration_path, lines)

    # -----------------------------------------------------------------

    @property
    def aliases(self):

        """
        This function ...
        :return:
        """

        alias_dict = dict()

        output = self.execute("alias")

        for line in output:

            splitted = line.split("=")

            if len(splitted) == 2: first, second = splitted[0], splitted[1]
            else:
                first = splitted[0]
                second = "".join(splitted[1:])

            alias = first.split("alias ")[1].strip()

            command = second
            if command.startswith("'") and command.endswith("'"): command = command[1:-1]
            elif command.startswith('"') and command.endswith('"'): command = command[1:-1]

            # Set the alias
            alias_dict[alias] = command

        # Return the dictionary
        return alias_dict

    # -----------------------------------------------------------------

    def resolve_alias(self, alias):

        """
        This function ...
        :param alias:
        :return:
        """

        return self.aliases[alias]

    # -----------------------------------------------------------------

    @property
    def python_path(self):

        """
        This function ...
        :return:
        """

        return self.find_executable("python")

    # -----------------------------------------------------------------

    @property
    def has_cpp_compiler(self):

        """
        This function ...
        :return:
        """

        for name in possible_cpp_compilers:
            if self.is_executable(name): return True
        return False

    # -----------------------------------------------------------------

    @property
    def cpp_compiler_path(self):

        """
        This function ...
        :return:
        """

        for name in possible_cpp_compilers:
            path = self.find_executable(name)
            if path is not None: return path
        return None

    # -----------------------------------------------------------------

    @property
    def has_mpi(self):

        """
        This function ...
        :return:
        """

        for name in possible_mpirun_names:
            if self.is_executable(name): return True
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def mpirun_path(self):

        """
        This function ...
        :return:
        """

        for name in possible_mpirun_names:
            path = self.find_executable(name)
            if path is not None: return path
        return None # No mpirun

    # -----------------------------------------------------------------

    @property
    def has_mpi_compiler(self):

        """
        This function ...
        :return:
        """

        for name in possible_mpi_compilers:
            if self.is_executable(name): return True

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_compiler_path(self):

        """
        This function ...
        :return:
        """

        for name in possible_mpi_compilers:
            path = self.find_executable(name)
            if path is not None: return path
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_version(self):

        """
        This function ...
        :return:
        """

        if not self.has_mpi: return None
        else: return self.version_of(self.mpirun_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def is_openmpi(self):

        """
        This function ...
        :return:
        """

        version_string = self.mpi_version
        return "Open MPI" in version_string or "OpenRTE" in version_string

    # -----------------------------------------------------------------

    @lazyproperty
    def openmpi_version(self):

        """
        This function ...
        :return:
        """

        if not self.is_openmpi: raise RuntimeError("Not using OpenMPI")
        return self.mpi_version.split(")")[1].strip()

    # -----------------------------------------------------------------

    def mpi_has_option(self, string):

        """
        This function ...
        :param string:
        :return:
        """

        output = self.execute(self.mpirun_path + " --" + string)

        # Check output
        for line in output:

            if "unrecognized argument " + string in line: return False # certainly not available
            if 'option "' + string + '" did not have enough parameters' in line: return True # certainly available

        # Probably available, found nothing suspect
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_has_cpus_per_proc_option(self):

        """
        This function ...
        :return:
        """

        return self.mpi_has_option("cpus-per-proc")

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_has_report_bindings_option(self):

        """
        This function ...
        :return:
        """

        return self.mpi_has_option("report-bindings")

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_has_map_by_option(self):

        """
        This function ...
        :return:
        """

        return self.mpi_has_option("map-by")

    # -----------------------------------------------------------------

    @lazyproperty
    def mpi_has_bind_to_option(self):

        """
        This function ...
        :return:
        """

        return self.mpi_has_option("bind-to")

    # -----------------------------------------------------------------

    def swap_cluster(self, cluster_name, check_modules=False, timeout=4):

        """
        This function ...
        :param cluster_name:
        :param check_modules:
        :param timeout:
        :return:
        """

        # Check if not already the loaded cluster
        if check_modules and "cluster/" + cluster_name in self.loaded_modules: return True

        # Swap to requested cluster
        try:
            self.execute("module swap cluster/" + cluster_name, timeout=timeout)
            return True
        except TimeOutReached:
            log.warning("Swapping to the '" + cluster_name + "' cluster failed: timeout reached")
            return False

    # -----------------------------------------------------------------

    def find_and_load_python(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        # Load module
        if self.has_lmod: module_name = self.load_python_module(show_output=show_output)
        else: module_name = None

        # Return
        if return_module: return self.python_path, module_name
        else: return self.python_path

    # -----------------------------------------------------------------

    @property
    def python_version_long(self):

        """
        This function ...
        :return:
        """

        return self.custom_python_version_long("python")

    # -----------------------------------------------------------------

    def custom_python_version_long(self, path, show_output=False):

        """
        This function ...
        :param path:
        :param show_output:
        :return:
        """

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Launch interactive python session
        self.ssh.sendline(path)

        #self.ssh.expect(">>>")
        index = self.ssh.expect([">>>", self.ssh.PROMPT, pexpect.EOF, "[PEXPECT]$"])
        #print(index)
        if index != 0: return None
        output = self.ssh.before.split("\r\n")[1:-1]

        # Close python again
        self.ssh.sendline("exit()")
        self.ssh.prompt()

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Get distribution and version, and architecture
        distribution_and_version = output[0].split("|")[0].split("(")[0].strip()
        if "|" in output[0]:
            # Compose version string with architecture
            architecture = output[0].split("|")[1].strip()
            return distribution_and_version + " " + architecture
        else: return distribution_and_version

    # -----------------------------------------------------------------

    @property
    def python_version_short(self):

        """
        This function ...
        :return:
        """

        return self.version_of(self.python_path)

    # -----------------------------------------------------------------

    def load_python_module(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        # Find Python module versions
        versions = self.get_module_versions("Python", show_output=show_output)

        latest_version = None
        latest_floating_python_version = None
        latest_intel_version = None
        latest_intel_year = None

        #print(versions)

        # Loop over the versions
        for version in versions:

            #print(version)

            # Get python version
            python_version = version.split("/")[1].split("-")[0]
            if int(python_version[0]) >= 3: continue # no python3

            compiler_version = version.split(python_version + "-")[1]

            if not ("intel" in compiler_version or "ictce" in compiler_version): continue

            splitted = python_version.split(".")
            floating_python_version = float(splitted[0]) + float(splitted[1]) / 100. + float(splitted[2]) / 10000.

            if latest_version is None or floating_python_version > latest_floating_python_version:
                latest_version = version
                latest_floating_python_version = floating_python_version

            elif floating_python_version == latest_floating_python_version:

                # Parse to get intel version
                intel_version = version.split("-")[2]

                if "." not in intel_version:

                    intel_year = int(intel_version) if len(intel_version) == 4 else int(intel_version[0:4])
                    if latest_intel_year is None or intel_year > latest_intel_year:
                        latest_intel_year = intel_year
                        latest_version = version

                elif latest_intel_year is not None:  # geef voorrang aan jaartallen
                    continue

                else:

                    if latest_intel_version is None or intel_version > latest_intel_version:
                        latest_version = version
                        latest_intel_version = intel_version

            else: continue

        # Return the latest version
        #return latest_version

        self.load_module(latest_version, show_output=show_output)  # Load the python module

        # Return the version
        return latest_version

    # -----------------------------------------------------------------

    def find_and_load_git(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        # Load intel compiler toolkit if possible
        if self.has_lmod:

            versions = self.get_module_versions("git", show_output=show_output)

            latest_version = None
            latest_module_name = None
            for version in versions:

                # Get version
                the_version = version.split("/")[1].split("-")[0]

                if latest_version is None or the_version > latest_version:
                    latest_version = the_version
                    latest_module_name = version

            if latest_module_name is None: raise RuntimeError("Git not available from the modules")

            # load the module
            loaded = self.load_module(latest_module_name, return_success=True)
            if not loaded: log.warning("Could not load the git module: finding system git installation ...")

        else: latest_module_name = None

        # Get path and version
        git_path = self.find_executable("git", show_output=show_output)
        git_version = self.version_of("git", show_output=show_output)

        # Return
        if return_module: return git_path, git_version, latest_module_name
        else: return git_path, git_version

    # -----------------------------------------------------------------

    def find_and_load_cpp_compiler(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        # Load intel compiler toolkit if possible
        if self.has_lmod: module_name = self.load_intel_compiler_toolkit(show_output=show_output)
        else: module_name = None

        # Search for the compiler and return its path
        if return_module: return self.cpp_compiler_path, module_name
        else: return self.cpp_compiler_path

    # -----------------------------------------------------------------

    def load_intel_compiler_toolkit(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        # Find latest Intel Compiler Toolkit version and load it
        intel_version = self._find_latest_iimpi_version_module(show_output=show_output)
        if intel_version is None: self.warning("Intel Cluster Toolkit Compiler Edition could not be found")
        elif intel_version not in self.loaded_modules:
            loaded = self.load_module(intel_version, show_output=show_output, return_success=True) # Load the module
            if not loaded: log.warning("The intel compiler module '" + intel_version + "' could not be loaded")

        # Return the module name
        return intel_version

    # -----------------------------------------------------------------

    def _find_latest_iimpi_version_module(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        # Load the Intel Cluster Toolkit Compiler Edition (inludes Intel MPI)
        versions = self.get_module_versions("iimpi", show_output=show_output)

        latest_version = None
        latest_iimpi_version = None
        latest_iimpi_year = None
        for version in versions:

            # Parse to get simple iimpi version
            iimpi_version = version.split("/")[1].split("-")[0]

            if "." not in iimpi_version:

                iimpi_year = int(iimpi_version) if len(iimpi_version) == 4 else int(iimpi_version[0:4])
                if latest_iimpi_year is None or iimpi_year > latest_iimpi_year:
                    latest_iimpi_year = iimpi_year
                    latest_version = version

            elif latest_iimpi_year is not None: # geef voorrang aan jaartallen
                continue

            else:

                if latest_iimpi_version is None or iimpi_version > latest_iimpi_version:
                    latest_version = version
                    latest_iimpi_version = iimpi_version

        # Return the latest version
        return latest_version

    # -----------------------------------------------------------------

    def find_and_load_mpi_compiler(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        if self.has_lmod: module_name = self.load_intel_compiler_toolkit(show_output=show_output)
        else: module_name = None

        # Return
        if return_module: return self.mpi_compiler_path, module_name
        else: return self.mpi_compiler_path

    # -----------------------------------------------------------------

    def find_and_load_qmake(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        # Use modules or not
        if self.has_lmod: return self._check_qt_remote_lmod(return_module=return_module, show_output=show_output)
        else:

            if return_module: return self._check_qt_remote_no_lmod(show_output=show_output), None
            else: return self._check_qt_remote_no_lmod(show_output=show_output)

    # -----------------------------------------------------------------

    def _check_qt_remote_lmod(self, return_module=False, show_output=False):

        """
        This function ...
        :param return_module:
        :param show_output:
        :return:
        """

        # Find latest Qt5 version
        qt5_version = self._find_latest_qt_version_module(show_output=show_output)

        # If no version could be found, give an error
        if qt5_version is None: raise RuntimeError("No Intel-compiled version of Qt 5 could be found between the available modules")

        # Load the module
        self.load_module(qt5_version, show_output=show_output)

        # Get the qmake path
        qmake_path = self.find_executable("qmake", show_output=show_output)

        # Get the Intel compiler version
        #if "intel" in qt5_version: intel_version = qt5_version.split("intel-")[1].split("-")[0]
        #else: intel_version = qt5_version.split("ictce-")[1].split("-")[0]

        # Return the qmake path
        if return_module: return qmake_path, qt5_version
        else: return qmake_path

    # -----------------------------------------------------------------

    def _check_qt_remote_no_lmod(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        # Keep the qmake paths in a list to decide later which one we can use
        qmake_paths = []

        # Search for qmake in the home directory
        command = "find " + self.home_directory + "/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.execute(command, show_output=show_output)

        # Search for qmake in the /usr/local directory
        command = "find /usr/local/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.execute(command, show_output=show_output)

        ## TOO SLOW?
        # Search for Qt directories in the home directory
        #for directory_path in self.remote.directories_in_path(self.remote.home_directory, startswith="Qt"):
            # Search for 'qmake' executables
            #for path in self.remote.files_in_path(directory_path, recursive=True):
                # Add the path
                #qmake_paths.append(path)

        ## TOO SLOW?
        # Search for Qt directories in /usr/local
        #for directory_path in self.remote.directories_in_path("/usr/local", startswith="Qt"):
            # Search for 'qmake' executables
            #for path in self.remote.files_in_path(directory_path, recursive=True):
                # Add the path
                #qmake_paths.append(path)

        # Check if qmake can be found by running 'which'
        qmake_path = self.find_executable("qmake", show_output=show_output)
        if qmake_path is not None: qmake_paths.append(qmake_path)

        latest_version = None

        latest_qmake_path = None

        # Loop over the qmake paths
        for qmake_path in qmake_paths:

            # Get the version
            output = self.execute(qmake_path + " -v", show_output=show_output)
            for line in output:
                if "Qt version" in line:
                    qt_version = line.split("Qt version ")[1].split(" in")[0]
                    break
            else: raise RuntimeError("Qt version could not be determined")

            if qt_version < "5.2.0": continue # oldest supported version
            if "conda" in qmake_path: continue
            if "canopy" in qmake_path: continue
            if "epd" in qmake_path: continue
            if "enthought" in qmake_path: continue

            if latest_version is None or qt_version > latest_version:

                latest_version = qt_version
                latest_qmake_path = qmake_path

        # Return the path
        return latest_qmake_path

    # -----------------------------------------------------------------

    def _find_latest_qt_version_module(self, show_output=False):

        """
        This function ...
        :param show_output:
        :return:
        """

        versions = []

        # Check 'Qt5' module versions
        versions += self.get_module_versions("Qt5", show_output=show_output)

        # Check 'Qt' module versions
        versions += self.get_module_versions("Qt", show_output=show_output)

        # Get the latest Qt version
        latest_version = None
        latest_qt_version = None
        for version in versions:

            # Only use Intel-compiled stuff
            if not ("ictce" in version or "intel" in version): continue

            # Parse to get simple Qt version
            qt_version = version.split("/")[1].split("-")[0]

            # Skip version below Qt 5
            if int(qt_version[0]) < 5: continue

            if latest_qt_version is None or qt_version > latest_qt_version:
                latest_version = version
                latest_qt_version = qt_version

        # Return the latest version
        return latest_version

    # -----------------------------------------------------------------

    def connect_to_vpn(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.info("Connecting to vpn service '" + self.host.vpn.service + "' ...")

        # Connect to the VPN service
        self.vpn = VPN(self.host.vpn.service)
        self.vpn.connect(self.host.vpn.user, self.host.vpn.password, self.host.vpn.secret, self.host.vpn.prompt_time_delay)

    # -----------------------------------------------------------------

    def login(self, login_timeout=20, silent=False):

        """
        This function ...
        :param login_timeout:
        :param silent:
        :return:
        """

        # Inform the user
        if not silent: log.info("Logging in to the remote environment on host '" + self.host.id + "' ...")

        # Try connecting to the remote host
        try: self.connected = self.ssh.login(self.host.name, self.host.user, self.host.password, port=self.host.port, login_timeout=login_timeout)
        except ExceptionPexpect: raise HostDownException(self.host.id)

        # Check whether connection was succesful
        if not self.connected: raise RuntimeError("Connection failed")

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        if log is not None and not self.silent: log.info("Logging out from the remote environment on host '" + self.host_id + "' ...")
        # the conditional statement is because of this error message during destruction at the end of a script:
        # Exception AttributeError: "'NoneType' object has no attribute 'info'" in <bound method SKIRTRemote.__del__ of <pts.core.simulation.remote.SKIRTRemote object at 0x118628d50>> ignored

        # Disconnect
        if self.connected:

            # TODO: Check if there are python sessions open?

            self.ssh.logout()
            self.connected = False

        # Disconnect from the VPN service if necessary
        #if self.vpn is not None: self.vpn.disconnect()

    # -----------------------------------------------------------------

    def start_screen(self, name, local_script_path, script_destination, screen_output_path=None,
                     keep_remote_script=True, attached=False):

        """
        This function ...
        :param name:
        :param local_script_path:
        :param script_destination:
        :param screen_output_path:
        :param keep_remote_script:
        :param attached:
        :return:
        """

        # Copy the script to the remote host
        self.upload(local_script_path, script_destination)

        # Rename the remote script
        local_script_name = fs.name(local_script_path)
        remote_script_name = name + ".sh"
        remote_script_path = fs.join(script_destination, remote_script_name)
        self.rename_file(script_destination, local_script_name, remote_script_name)

        # Make the shell script executable
        self.execute("chmod +x " + remote_script_path, output=False)

        # Record the screen output: 'script' command
        #if screen_output_path is not None: self.execute("script " + screen_output_path)

        # Create the screen session and execute the batch script
        if attached:
            #self.execute("screen -S " + name + " -L -m " + remote_script_path, output=False, show_output=True)
            self.execute("sh " + remote_script_path, output=False, show_output=True)
        else: self.execute("screen -S " + name + " -L -d -m " + remote_script_path, output=False, timeout=None, cwd=screen_output_path)

        # Remove the remote shell script
        if not keep_remote_script: self.execute("rm " + remote_script_path, output=False)

        # Return the remote screen script file path
        return remote_script_path

    # -----------------------------------------------------------------

    def kill_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Quit the specified screen session
        self.execute("screen -S " + name + " -X quit", output=False)

    # -----------------------------------------------------------------

    def create_screen_session(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        command = "screen -dmS " + name
        self.execute(command)

    # -----------------------------------------------------------------

    def close_screen_session(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.kill_screen(name)

    # -----------------------------------------------------------------

    def close_all_screen_sessions(self):

        """
        This function ...
        :return:
        """

        for name in self.screen_names(): self.close_screen_session(name)

    # -----------------------------------------------------------------

    def kill_job(self, job_id, cluster_name=None, show_output=False):

        """
        This function ...
        :param job_id:
        :param cluster_name:
        :param show_output:
        :return:
        """

        # If cluster name is given, switch
        original_cluster_name = None
        if cluster_name is not None:
            active_cluster = self.active_cluster_name
            if active_cluster != cluster_name:
                original_cluster_name = active_cluster
                self.swap_cluster(cluster_name)

        # Stop the job with the specified ID
        self.execute("qdel " + str(job_id), output=False, show_output=show_output)

        # Swap back to the original cluster?
        if original_cluster_name is not None: self.swap_cluster(original_cluster_name)

    # -----------------------------------------------------------------

    def tmux_names(self):

        """
        This function ...
        :return:
        """

        if not self.is_executable("tmux"): return []
        output = self.execute("tmux ls")
        if len(output) == 1 and output[0] == "failed to connect to server": return []
        else: return [line.split(":")[0] for line in output]

    # -----------------------------------------------------------------

    def create_tmux_session(self, name, in_session_command=None):

        """
        This function ...
        :param name:
        :param in_session_command:
        :return:
        """

        command = "tmux new -d -s " + name
        if in_session_command is not None: command += " " + in_session_command
        self.execute(command)

    # -----------------------------------------------------------------

    def close_tmux_session(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        end_session_command = "tmux kill-session -t " + name
        output = self.execute(end_session_command)

        for line in output:
            if "session not found" in line: raise ValueError("Session does not exist")

    # -----------------------------------------------------------------

    def close_all_tmux_sessions(self):

        """
        This function ...
        :return:
        """

        for name in self.tmux_names(): self.close_tmux_session(name)

    # -----------------------------------------------------------------

    def screen_names(self):

        """
        This function ...
        :return:
        """

        output = self.execute("screen -ls")

        # Get the names
        names = []
        for line in output:

            # Skip lines that do not state a session
            if "(Attached)" not in line and "(Detached)" not in line: continue

            # Get the name
            name = line.split(".")[1].split("\t")[0]
            names.append(name)

        # Return the screen names
        return names

    # -----------------------------------------------------------------

    def screen_numbers(self):

        """
        This function ...
        :return:
        """

        output = self.execute("screen -ls")

        # Get the numbers
        numbers = []
        for line in output:

            # Skip lines that do not state a session
            if "(Attached)" not in line and "(Detached)" not in line: continue

            # Get the number
            number = int(line.split(".")[0].split("\t")[1])
            numbers.append(number)

        # Return the screen numbers
        return numbers

    # -----------------------------------------------------------------

    def screen_sessions(self):

        """
        This function ...
        :return:
        """

        output = self.execute("screen -ls")

        # Initialize dictionary
        screens = dict()

        # Loop over the lines
        for line in output:

            # Skip lines that do not state a session
            if "(Attached)" not in line and "(Detached)" not in line: continue

            # Get the name
            name = line.split(".")[1].split("\t")[0]

            # Get the number
            number = int(line.split(".")[0].split("\t")[1])

            # Get the timestamp
            timestamp = line.split("\t")[2][1:-1]

            # Get the state
            state = line.split("\t")[3][1:-1]

            # Create map and add to the dictionary
            screens[name] = Map(number=number, timestamp=timestamp, state=state)

        # Return the dictionary
        return screens

    # -----------------------------------------------------------------

    def screen_states(self):

        """
        This function ...
        :return:
        """

        states = dict()
        for name in self.screen_names(): states[name] = self.screen_state(name)
        return states

    # -----------------------------------------------------------------

    def screen_state(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Execute the 'screen -ls' command
        output = self.execute("screen -ls")

        # Loop over the different active screen sessions
        for entry in output:

            # Check if the specified screen name corresponds to the current entry
            if name in entry:

                # Check the state of this screen session
                if "Detached" in entry: return "detached"
                elif "Attached" in entry: return "attached"
                else: raise ValueError("Screen " + name + " in unkown state")

        # If the screen name was not encountered, return None (the screen session has finished or has been aborted)
        return None

    # -----------------------------------------------------------------

    def is_active_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        state = self.screen_state(name)
        return state == "detached" or state == "attached"

    # -----------------------------------------------------------------

    def is_attached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "attached"

    # -----------------------------------------------------------------

    def is_detached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "detached"

    # -----------------------------------------------------------------

    def screen_number(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Execute the 'screen -ls' command
        output = self.execute("screen -ls | grep " + name)

        if len(output) > 1: raise ValueError("Multiple screens with name '" + name + "'")

        line = output[0]
        return int(line.split(".")[0])

    # -----------------------------------------------------------------

    def get_requirements(self, processors):

        """
        This function calculates the required amount of nodes and processors per node, given a certain number of
        processors.
        :param processors:
        :return:
        """

        # Calculate the necessary amount of nodes
        nodes = processors // self.cores_per_node + (processors % self.cores_per_node > 0)

        # Determine the number of processors per node
        ppn = processors if nodes == 1 else self.cores_per_node

        # Return the number of nodes and processors per node
        return nodes, ppn

    # -----------------------------------------------------------------

    def launch_pts_command(self, command, arguments, show_output=True):

        """
        This function ...
        :param command:
        :param arguments:
        :param show_output:
        :return:
        """

        # Convert the list of commands back to a single string
        portions = []
        for argument in arguments:
            if " " in argument: portions.append("'" + argument + "'")
            else: portions.append(argument)
        argument_string = " ".join(portions)

        # Execute the command
        self.execute("pts " + command + " " + argument_string, show_output=show_output)

    # -----------------------------------------------------------------

    def download_from_url_to(self, url, path, show_output=True, overwrite=False):

        """
        This function ...
        :param url:
        :param path:
        :param show_output:
        :param overwrite:
        :return:
        """

        # Download to directory of file path
        if self.is_directory(path): return self.download_from_url_to_directory(url, path, show_output, overwrite)
        else:
            self.download_from_url_to_file(url, path, show_output, overwrite)
            return path

    # -----------------------------------------------------------------

    def download_from_url_to_file(self, url, filepath, show_output=True, overwrite=False):

        """
        This function ...
        :param url:
        :param filepath:
        :param show_output:
        :param overwrite:
        :return:
        """

        # Is an existing file
        if self.is_file(filepath):
            if overwrite: self.remove_file(filepath)
            else: raise IOError("File already exists: " + filepath)

        # Execute the download command
        command = "wget " + url + " -O " + filepath
        self.execute(command, show_output=show_output)

        # Check
        if not self.is_file(filepath): raise RuntimeError("Something went wrong downloading the file '" + url + "' to '" + filepath + "'")

    # -----------------------------------------------------------------

    def download_from_url_to_directory(self, url, path, show_output=True, overwrite=False):

        """
        This function ...
        :param url:
        :param path:
        :param show_output:
        :param overwrite:
        :return:
        """

        # Determine the path of the downloaded file
        filename = fs.name(url)
        filepath = fs.join(path, filename)

        # If is already a file
        if self.is_file(filepath):
            filepath = fs.join(path, time.filename_timestamp() + "_" + filename)

        # Download to file path
        self.download_from_url_to_file(url, filepath, show_output=show_output, overwrite=overwrite)

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def run_pts(self, command, config, input_dict=None, keep_remote_output=False, remove_local_output=False):

        """
        This function ...
        :param command:
        :param config:
        :param input_dict:
        :param keep_remote_output:
        :param remove_local_output:
        :return:
        """

        # Create a name for this PTS session
        unique_session_name = time.unique_name(command.replace("/", "-"))

        # Create a remote temporary directory
        remote_temp_path = self.new_temp_directory()

        ##

        ## CHANGE THE LOG PATH TO A REMOTE PATH AND CHANGE THE CWD

        # Get the local working directory path
        local_cwd = config.path

        # Always create a log file while executing remotely
        config.report = True
        config.log_path = remote_temp_path
        config.path = remote_temp_path

        # Don't save the config remotely again (we will upload it explicitly)
        config.config_path = None

        config.write_config = False

        ##

        # DETERMINE REMOTE OUTPUT PATH

        # If the 'output' parameter is specified in the configuration, get the full absolute path to the corresponding output directory
        if "output" in config:

            # Determine the local output path that is desired by the user
            local_output_path = fs.absolute_or_in(config.output, local_cwd)

            # Create the local output path if necessary
            if not fs.is_directory(local_output_path): fs.create_directory(local_output_path)

            # Set the remote output path
            remote_output_path = fs.join(config.path, "out")

            # Change the out path for running on the remote
            config.output = remote_output_path

        # else, the output directory is the remote working directory
        # THIS IS FOR 'OLDER' CONFIGURABLE DERIVED CLASSES THAT DON'T USE THE 'OUTPUT_PATH' and 'OUTPUT_PATH_FILE' MECHANICS ...
        else:

            # Set local output path and remote output path
            local_output_path = local_cwd
            remote_output_path = config.path

        ###

        # Debugging
        self.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, unique_session_name + ".cfg")

        # Save the configuration file to the temporary directory
        config.saveto(temp_conf_path)

        # Debugging
        self.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

        # Upload the config file
        remote_conf_path = fs.join(remote_temp_path, fs.name(temp_conf_path))
        self.upload(temp_conf_path, remote_temp_path)

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        #### UPLOAD THE INPUT

        if input_dict is not None:

            # Create remote input path
            input_path = fs.join(remote_temp_path, "input")
            self.create_directory(input_dict)

            # TODO: use the write_input function in the configurable module instead

            for name in input_dict:

                filename = name + "." + input_dict[name].default_extension

                local_filename = fs.join(temp_path, filename) # local temporary path

                # TODO: save the input item here

                filepath = fs.join(input_path, filename) # remote temporary path

                # TODO: UPLOAD THE INPUT item here

                # TODO: delete the local input item

        ####

        # Debugging
        self.debug("Creating a script for remote execution ...")

        # Determine the path to the remote equivalent of this file
        remote_main_path = fs.join(self.pts_package_path, "do", "__main__.py")

        # Create a bash script
        temp_script_path = fs.join(temp_path, unique_session_name + ".sh")

        # Determine python path
        python_path = self.conda_pts_environment_python_path

        # Write the lines to the script file
        with open(temp_script_path, 'w') as script_file:

            #script_file.write("#!/usr/bin/env python\n")
            #script_file.write("# -*- coding: utf8 -*-\n")
            #script_file.write("\n")
            script_file.write(python_path + " " + remote_main_path + " --configfile " + remote_conf_path + " " + command + "\n")

        # Execute the script
        # name, local_script_path, script_destination, screen_output_path=None, keep_remote_script=False
        self.start_screen(unique_session_name, temp_script_path, remote_temp_path, screen_output_path=remote_temp_path, keep_remote_script=keep_remote_output)

        # Remove the local script file
        fs.remove_file(temp_script_path)

        # Remove the remote temporary directory
        if keep_remote_output: self.info("Remote output will be placed in '" + remote_output_path + "'")

        # Create a new Task object
        task = Task(command, config.to_string())

        # Set the host ID and cluster name (if applicable)
        task.host_id = self.host_id
        task.cluster_name = self.cluster_name

        # Generate a new task ID
        task_id = self._new_task_id()

        # Determine the path to the task file
        task_file_path = fs.join(self.local_pts_host_run_dir, str(task_id) + ".task")
        task.path = task_file_path

        # Set properties such as the task ID and name and the screen name
        task.id = task_id
        task.remote_temp_pts_path = remote_temp_path
        task.name = unique_session_name
        task.screen_name = unique_session_name
        task.remote_screen_output_path = remote_temp_path

        # Set local and remote output path
        task.local_output_path = local_output_path
        task.remote_output_path = remote_output_path

        # Other
        task.remove_remote_output = not keep_remote_output
        task.remove_local_output = remove_local_output

        # Save the task
        task.save()

        # Return the task
        return task

    # -----------------------------------------------------------------

    def clear_commands(self):

        """
        This function ...
        :return:
        """

        self.commands = []

    # -----------------------------------------------------------------

    def add_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Debugging
        self.debug('Sent command: "'  + command + '"')

        # Add the command
        self.commands.append(command)

    # -----------------------------------------------------------------

    def make_executable(self, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Make executable
        self.execute("chmod +rx '" + filepath + "'")

    # -----------------------------------------------------------------

    def run_script(self, filepath, options, output=True, show_output=False, timeout=None, expect=None):

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

        # Detemrine directory path and filename
        dir_path = fs.directory_of(filepath)
        filename = fs.name(filepath)

        self.make_executable(filepath)
        return self.execute("./" + filename + " " + options, output=output, show_output=show_output, timeout=timeout,
                            expect=expect, cwd=dir_path)

    # -----------------------------------------------------------------

    def execute(self, command, output=True, expect_eof=True, contains_extra_eof=False, show_output=False, timeout=None,
                expect=None, cwd=None, output_start=1):

        """
        This function ...
        :param command:
        :param output:
        :param expect_eof:
        :param contains_extra_eof:
        :param show_output:
        :param timeout:
        :param expect:
        :param cwd:
        :param output_start:
        :return:
        """

        # Check whether connected
        if not self.connected: raise RuntimeError("The remote is not available")

        # Change the working directory if necessary
        if cwd is not None: original_cwd = self.change_cwd(cwd)
        else: original_cwd = None

        # Send the command
        self.ssh.sendline(command)

        # Add the command to the list of commands
        self.add_command(command)

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Retrieve the output if requested
        if expect is None:

            try: matched = self.ssh.prompt(timeout=timeout)
            except pexpect.EOF: raise RuntimeError("End-of-file encountered: connection broken?")

            # If an extra EOF is used before the actual output line (don't ask me why but I encounter this on the HPC UGent infrastructure), do prompt() again
            if contains_extra_eof: matched = self.ssh.prompt()

            # If the command could not be sent, raise an error
            if not matched and expect_eof and not contains_extra_eof: raise TimeOutReached("Time out for command '" + command + "'", timeout=timeout, command=command)

        # Check for expected characters
        else:

            try: index = self.ssh.expect(expect, timeout=timeout)
            except pexpect.EOF: raise RuntimeError("End-of-file encountered: connection broken?")
            assert index == 0

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Ignore the first and the last line (the first is the command itself, the last is always empty)
        # Trial and error to get it right for HPC UGent login nodes; don't know what is happening
        if contains_extra_eof:
            splitted = self.ssh.before.replace('\x1b[K', '').split("\r\n")
            for i in range(len(splitted)): splitted[i] = splitted[i].replace(" \r", "").replace("\x08", "").replace("\xe2\x80\x98", "'").replace("\xe2\x80\x99", "'")
            if splitted[-1] == "": the_output = splitted[output_start:-1]
            else: the_output = splitted[output_start:]
        else:
            splitted = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")
            for i in range(len(splitted)): splitted[i] = splitted[i].replace(" \r", "").replace("\x08", "").replace("\xe2\x80\x98", "'").replace("\xe2\x80\x99", "'")
            if splitted[-1] == "": the_output = splitted[output_start:-1]
            else: the_output = splitted[output_start:]

        # Set the working directory back to the original
        if original_cwd is not None: self.change_cwd(original_cwd)

        # Check whether the output contains the (env_name) string from activated conda environment
        if len(the_output) > 0 and the_output[-1].strip().startswith("(") and the_output[-1].strip().endswith(")"):
            the_output = the_output[:-1]

        # Return the output
        if output: return the_output

    # -----------------------------------------------------------------

    def execute_lines(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Get arguments
        output = kwargs.pop("output", True)
        show_output = kwargs.pop("show_output", False)
        timeout = kwargs.pop("timeout", None)
        cwd = kwargs.pop("cwd", None)

        # Change the working directory if necessary
        if cwd is not None: original_cwd = self.change_cwd(cwd)
        else: original_cwd = None

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        #print(args)

        # Loop over the lines
        for line in args:

            #print("LINE:", line)

            #print("BEFORE", self.ssh.before)
            #print("AFTER", self.ssh.after)

            # If string, send the command
            if types.is_string_type(line): self.ssh.sendline(line)

            # If tuple
            elif isinstance(line, tuple):

                # Expect
                if len(line) == 3 and line[2]:

                    # index = self.ssh.expect([self.ssh.PROMPT, "$", line[0]], timeout=timeout)
                    #index = self.ssh.expect([line[0], "$", self.ssh.PROMPT], timeout=timeout)
                    #print("here", self.ssh.PROMPT)
                    index = self.ssh.expect([line[0], self.ssh.PROMPT, "[PEXPECT]$"], timeout=timeout)
                    #print(index)
                    if index == 0: self.ssh.sendline(line[1])
                    else: pass

                    #if index == 0: pass
                    #elif index == 1: pass
                    #elif index == 2: self.ssh.sendline(line[1])

                else:

                    self.ssh.expect(line[0], timeout=timeout)
                    self.ssh.sendline(line[1])

            # Invalid
            else: raise ValueError("Lines must be strings or tuples")

        # Synchronize
        eof = self.ssh.prompt(timeout=timeout)

        # If the command could not be sent, raise an error
        if not eof: raise RuntimeError("The commands could not be sent")

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Return the output
        the_output = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]

        # Set the working directory back to the original
        if original_cwd is not None: self.change_cwd(original_cwd)

        # Return the output
        if output: return the_output

    # -----------------------------------------------------------------

    def in_python_virtual_environment(self):

        """
        This function ...
        :return:
        """

        # python -c 'import sys; print (sys.real_prefix)'

        lines = []
        lines.append("import sys")
        lines.append("answer = False")
        lines.append("answer = True if hasattr(sys, 'real_prefix') else answer")
        lines.append("answer = True if 'conda' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'canopy' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'epd' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'enthought' in sys.prefix.lower() else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        if output[0] == "True": return True
        elif output[0] == "False": return False
        else: raise RuntimeError("Something went wrong: answer from remote = '" + str(output))

    # -----------------------------------------------------------------

    def start_python_session(self, assume_pts=True, output_path=None, attached=False, new_connection_for_attached=False):

        """
        This function ...
        :param assume_pts: assume PTS is present, import some basic PTS tools
        :param output_path:
        :param attached:
        :param new_connection_for_attached:
        :return:
        """

        # Check whether we find a python installation made for PTS
        if self.conda_environment_for_pts is not None:
            python_path = fs.join(self.conda_installation_path, "envs", self.conda_environment_for_pts, "bin", "python")
        else: python_path = "python"

        # Start python session and return it
        if attached:
            if new_connection_for_attached:
                new_remote = Remote(host_id=self.host_id)
                return AttachedPythonSession(new_remote, python_path, assume_pts=assume_pts, output_path=output_path)
            else: return AttachedPythonSession(self, python_path, assume_pts=assume_pts, output_path=output_path)
        else: return DetachedPythonSession(self, python_path, assume_pts=assume_pts, output_path=output_path)

    # -----------------------------------------------------------------

    def execute_python_interactive(self, lines, output_path=None, attached=False):

        """
        This function ...
        :param lines:
        :param output_path:
        :param attached:
        :return:
        """

        # Create python session
        python = self.start_python_session(output_path=output_path, attached=attached)

        # Execute the lines
        output = python.send_lines(lines)

        # Return the output
        return output

    # -----------------------------------------------------------------

    def touch(self, path):

        """
        Touch a file
        :param path:
        :return:
        """

        self.execute("touch " + path)

    # -----------------------------------------------------------------

    def touch_alternative(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.write_line(path, "")

    # -----------------------------------------------------------------

    def rename_file(self, directory, old_name, new_name):

        """
        This function ...
        :param directory:
        :param old_name:
        :param new_name:
        :return:
        """

        # Determine the old and new file path
        old_path = fs.join(directory, old_name)
        new_path = fs.join(directory, new_name)

        # Check if file exists
        if not self.is_file(old_path): raise IOError("Not a file: '" + old_path + "'")

        # Use the 'mv' command to rename the file
        self.execute("mv " + old_path + " " + new_path)

        # Return the new filepath
        return new_path

    # -----------------------------------------------------------------

    def rename_directory(self, parent, old_name, new_name):

        """
        This function ...
        :param parent:
        :param old_name:
        :param new_name:
        :return:
        """

        # Check if directory exists
        path = fs.join(parent, old_name)
        new_path = fs.join(parent, new_name)
        if not self.is_directory(path): raise IOError("Not a directory: '" + path + "'")

        # Move
        self.execute("mv " + path + " " + new_path)

    # -----------------------------------------------------------------

    def remove_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check if directory
        if not self.is_directory(path): raise ValueError("Not a directory: '" + path + "'")

        # Debugging
        self.debug("Removing directory '" + path + "' ...")

        # Execute the command
        output = self.execute("rm -rf '" + path + "'", output=True)

        # Check for error
        for line in output:
            if "cannot remove" in line: raise RuntimeError(line)

    # -----------------------------------------------------------------

    def clear_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        self.debug("Clearing directory '" + path + "' ...")

        # Remove all files
        for file_path in self.files_in_path(path):
            print(file_path)
            self.remove_file(file_path)

        # Remove all d irectories
        for directory_path in self.directories_in_path(path): self.remove_directory(directory_path)

    # -----------------------------------------------------------------

    def remove_file(self, path, show_output=False):

        """
        This function ...
        :param path:
        :param show_output:
        :return:
        """

        # Check if file
        if not self.is_file(path): raise ValueError("Not a file: '" + path + "'")

        # Debugging
        self.debug("Removing file '" + path + "' ...")

        # Execute the command
        output = self.execute("rm '" + path + "'", output=True, show_output=show_output)

        # Check for error
        for line in output:
            if "cannot remove" in line: raise RuntimeError(line)

    # -----------------------------------------------------------------

    def remove_files(self, paths, show_output=False):

        """
        This function ...
        :param paths:
        :param show_output:
        :return:
        """

        for path in paths: self.remove_file(path)

    # -----------------------------------------------------------------

    def remove_files_in_path(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        show_output = kwargs.pop("show_output", False)
        self.remove_files(self.files_in_path(*args, **kwargs), show_output=show_output)

    # -----------------------------------------------------------------

    def remove_directories(self, paths):

        """
        This function ...
        :param paths:
        :return:
        """

        for path in paths: self.remove_directory(path)

    # -----------------------------------------------------------------

    def remove_directories_in_path(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self.remove_directories(self.directories_in_path(*args, **kwargs))

    # -----------------------------------------------------------------

    def change_cwd(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get current working directory
        original_cwd = self.working_directory

        # Try to change the directory, give an error if this fails
        #output = self.execute("cd " + path)
        #if len(output) > 0 and "No such file or directory" in output[0]: raise RuntimeError("The directory does not exist")

        #return original_cwd

        # Send the command
        self.ssh.sendline("cd '" + path + "'")
        self.ssh.prompt()
        output = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]
        if len(output) > 0 and "No such file or directory" in output[0]: raise RuntimeError("The directory '" + path + "' does not exist")

        # Return the original working directory
        return original_cwd

    # -----------------------------------------------------------------

    def items_in_path(self, path, recursive=False):

        """
        This function ...
        :param path:
        :param recursive:
        :return:
        """

        if recursive:

            command = "ls -R '" + path + "' | awk '/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\n'
            command += "NF&&f{ print "
            command += 's"/"$0 }'
            command += "'"

            output = self.execute(command, cwd=path)

            paths = [line for line in output if not line.startswith(">")]

            # Return the paths
            return paths

        else:

            directories = self.execute("for i in $(ls */); do echo ${i%%/}; done")
            if "cannot access */" in directories[0]: directories = []

            files = self.execute("for f in *; do [[ -d $f ]] || echo $f; done", cwd=path)

            # Get paths
            items = directories + files
            paths = [fs.join(path, name) for name in items]

            # Return the paths
            return paths

    # -----------------------------------------------------------------

    def directories_in_path(self, path, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                        returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None):

        """
        This function ...
        :param path:
        :param recursive:
        :param ignore_hidden:
        :param contains:
        :param not_contains:
        :param returns:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :return:
        """

        if not ignore_hidden: self.warning("Not ignoring hidden folders is not yet supported")

        # List the directories in the provided path
        if recursive:
            command = "ls -R -d */ | awk '\n/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\nNF&&f{ print '
            command += 's"/"$0 }'
            command += "'"
            output = self.execute(command, cwd=path)
            paths = [line for line in output if not line.startswith(">")]
            paths = [dirpath for dirpath in paths if self.is_directory(dirpath)]
        else:
            #print(path)
            output = self.execute("for i in $(ls -d */); do echo ${i%%/}; done", cwd=path)
            #print(output)
            if "cannot access */" in output[0]: return []
            elif "cannot access '*/'" in output[0]: return []
            paths = [fs.join(path, name) for name in output]

        if returns == "dict":
            returns = ["name", "path"]
            return_dict = True
        else: return_dict = False

        # Filter
        result = []

        # List all items in the specified directory
        for item_path in paths:

            # Get the name of the directory
            item = fs.name(item_path)

            # Get the path to the containing directory
            directory = fs.directory_of(item_path)

            # Ignore hidden directories if requested
            if ignore_hidden and item.startswith("."): continue

            # Ignore names that do not contain a certain string, if specified
            if contains is not None and contains not in item: continue

            # Ignore names that do contain a certain string that it should not contain, if specified
            if not_contains is not None and not_contains in item: continue

            # If the directory name does not match the exact name, skip it
            if exact_name is not None and exact_name != item: continue

            # If the directory name matches the 'exact not name', skip it
            if exact_not_name is not None:

                if types.is_string_type(exact_not_name):
                    if exact_not_name == item: continue
                elif isinstance(exact_not_name, list):
                    if item in exact_not_name: continue
                else: raise ValueError("Invalid option for 'exact_not_name': must be string or list")

            # Ignore directory names that do not start or end with the specified strings
            if startswith is not None and not item.startswith(startswith): continue
            if endswith is not None and not item.endswith(endswith): continue

            # Create the return value
            if types.is_string_type(returns):

                if returns == "path": thing = item_path
                elif returns == "name": thing = item
                elif returns == "directory": thing = directory
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

            else:  # Assume 'returns' is a list or tuple

                thing = []

                for return_value in returns:

                    if return_value == "path": thing.append(item_path)
                    elif return_value == "name": thing.append(item)
                    elif return_value == "directory": thing.append(directory)
                    else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

            # Add an entry to the result
            result.append(thing)

        # Return the result
        if return_dict: return dict(result)
        else: return result

    # -----------------------------------------------------------------

    def find_directory_in_path(self, path, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                               exact_name=None, exact_not_name=None, startswith=None, endswith=None):

        """
        This function ...
        :param path: 
        :param recursive: 
        :param ignore_hidden: 
        :param contains: 
        :param not_contains: 
        :param exact_name: 
        :param exact_not_name: 
        :param startswith: 
        :param endswith: 
        :return: 
        """

        paths = self.directories_in_path(path, recursive=recursive, ignore_hidden=ignore_hidden, contains=contains,
                                    not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name,
                                    startswith=startswith, endswith=endswith)
        if len(paths) == 1: return paths[0]
        elif len(paths) == 0: raise ValueError("Not found")
        else: raise ValueError("Multiple directories found")

    # -----------------------------------------------------------------

    def files_in_path(self, path, recursive=False, ignore_hidden=True, extension=None, not_extension=None,
                      contains=None, not_contains=None, exact_name=None, exact_not_name=None, startswith=None,
                      endswith=None, returns="path", extensions=False):

        """
        This function ...
        :param path:
        :param recursive:
        :param ignore_hidden:
        :param extension:
        :param not_extension:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :param returns:
        :param extensions:
        :return:
        """

        if not ignore_hidden: self.warning("Not ignoring hidden files is not yet supported")

        # List the files in the provided path
        if recursive:
            command = "ls -R '" + path + "' | awk '\n/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\nNF&&f{ print '
            command += 's"/"$0 }'
            command += "'"
            output = self.execute(command, cwd=path)
            #print(output)
            paths = [line for line in output if not line.startswith(">")]
            paths = [filepath for filepath in paths if self.is_file(filepath)]
        else:
            output = self.execute("for f in *; do [[ -d $f ]] || echo $f; done", cwd=path)
            if len(output) == 1 and output[0] == "*": return []
            paths = [fs.join(path, name) for name in output]

        if returns == "dict":
            returns = ["name", "path"]
            return_dict = True
        else: return_dict = False

        #print(paths)

        # Filter
        result = []

        # Loop over the paths
        for path in paths:

            item_path = path

            # Determine the name of the item
            item = fs.name(path)

            # Get the path to the containing directory
            directory = fs.directory_of(item_path)

            # Get the file name and extension
            item_name = fs.strip_extension(item)
            item_extension = fs.get_extension(item)

            # Ignore files with extension different from the one that is specified
            if extension is not None:
                if types.is_string_type(extension):
                    if item_extension != extension: continue
                elif types.is_string_sequence(extension):
                    if item_extension not in extension: continue
                else: raise ValueError("Unknown type for 'extension': " + str(extension))

            # Ignore files with extensions that are in 'not_extension'
            if not_extension is not None:
                if types.is_string_type(not_extension):
                    if item_extension == not_extension: continue
                elif types.is_string_sequence(not_extension):
                    if item_extension in not_extension: continue
                else: raise ValueError("Unknown type for 'not_extension': " + str(not_extension))

            # Ignore filenames that do not contain a certain string, if specified
            if contains is not None and contains not in item_name: continue

            # Ignore filenames that do contain a certain string that it should not contain, if specified
            if not_contains is not None and not_contains in item_name: continue

            # Ignore filenames that do not match the exact filename, if specified
            if exact_name is not None and exact_name != item_name: continue

            # If the filename matches the 'exact not name', skip it
            if exact_not_name is not None and exact_not_name == item_name: continue

            # Ignore filenames that do not start or end with the specified strings
            if startswith is not None and not item_name.startswith(startswith): continue
            if endswith is not None and not item_name.endswith(endswith): continue

            # Create the return value
            if types.is_string_type(returns):

                if returns == "path": thing = item_path
                elif returns == "name": thing = item_name + "." + item_extension if extensions else item_name
                elif returns == "directory": thing = directory
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

            else:  # Assume 'returns' is a list

                thing = []

                for return_value in returns:

                    if return_value == "path": thing.append(item_path)
                    elif return_value == "name": thing.append(item_name + "." + item_extension if extensions else item_name)
                    elif return_value == "directory": thing.append(directory)
                    else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

            # Add entry to the result
            result.append(thing)

        # Return the result
        if return_dict: return dict(result)
        else: return result

    # -----------------------------------------------------------------

    def nfiles_in_path(self, *args, **kwargs):

        """
        Thisf unction ...
        :param args:
        :param kwargs:
        :return:
        """

        return len(self.files_in_path(*args, **kwargs))

    # -----------------------------------------------------------------

    def has_files_in_path(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self.nfiles_in_path(*args, **kwargs) > 0

    # -----------------------------------------------------------------

    def contains_files(self, directory, filenames=None):

        """
        This function ...
        :param directory:
        :param filenames:
        :return:
        """

        # Filenames are given
        if filenames is not None:

            # Loop over the filenames
            for filename in filenames:

                # Check
                filepath = fs.join(directory, filename)
                if not self.is_file(filepath): return False

            # All checks passed
            return True

        # No filenames are given
        else: return self.has_files_in_path(directory)

    # -----------------------------------------------------------------

    def find_file_in_path(self, path, recursive=False, ignore_hidden=True, extension=None, contains=None, not_contains=None,
                            exact_name=None, exact_not_name=None, startswith=None, endswith=None):

        """
        This function ...
        :param path:
        :param recursive:
        :param ignore_hidden:
        :param extension:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :return: 
        """

        # Get paths
        paths = self.files_in_path(path, recursive=recursive, ignore_hidden=ignore_hidden, extension=extension,
                              contains=contains, not_contains=not_contains, exact_name=exact_name,
                              exact_not_name=exact_not_name, startswith=startswith, endswith=endswith)
        if len(paths) == 1: return paths[0]
        elif len(paths) == 0: raise ValueError("Not found")
        else: raise ValueError("Multiple files found")

    # -----------------------------------------------------------------

    def copy_file(self, file_path, directory_path, new_name=None):

        """
        This function ...
        :param file_path:
        :param directory_path:
        :param new_name:
        :return:
        """

        if new_name is not None: destination = fs.join(directory_path, new_name)
        else: destination = directory_path

        # Copy
        self.execute("cp '" + file_path + "' '" + destination + "'")

        # Return
        if self.is_file(destination): return destination
        elif self.is_directory(destination): return fs.join(destination, fs.name(file_path))
        else: raise ValueError("Don't understand the destination: " + destination)

    # -----------------------------------------------------------------

    def copy_directory(self, path, directory_path, new_name=None):

        """
        This function ...
        :param path:
        :param directory_path:
        :param new_name:
        :return:
        """

        # Create the directory
        dirname = new_name if new_name is not None else fs.name(path)
        copy_path = self.create_directory_in(directory_path, dirname)

        # Copy contents
        self.copy_from_directory(path, copy_path)

    # -----------------------------------------------------------------

    def copy_from_directory(self, from_directory, to_directory, **kwargs):

        """
        This function ...
        :param from_directory:
        :param to_directory:
        :param kwargs:
        :return:
        """

        self.copy_files(self.files_in_path(from_directory, **kwargs), to_directory)

    # -----------------------------------------------------------------

    def copy_files(self, file_paths, directory_path):

        """
        This function ...
        :param file_paths:
        :param directory_path:
        :return:
        """

        for file_path in file_paths: self.copy_file(file_path, directory_path)

    # -----------------------------------------------------------------

    def directory_size(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..units.parsing import parse_quantity

        command = "du -sh '" + path + "'"
        output = self.execute(command)

        string = output[0].split(" ")[0].split("\t")[0].strip()

        if string.endswith("G") or string.endswith("M"): string += "B"
        elif string.lower().endswith("terra") or string.lower().endswith("giga") or string.lower().endswith("mega"): string = string.lower() + "byte"
        else: string += "byte" # guess

        # Parse the quantity and convert to GB
        return parse_quantity(string).to("GB")

    # -----------------------------------------------------------------

    def file_size(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..units.parsing import parse_quantity

        command = "du -sh '" + path + "'"
        output = self.execute(command)

        string = output[0].split(" ")[0].lower() + "byte"
        return parse_quantity(string)

    # -----------------------------------------------------------------

    def file_nbytes(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        command = "du -b '" + path + "' | awk '{ print $1}' | bc"
        output = self.execute(command)
        for line in output:
            if "No such file or directory" in line: raise IOError("File does not exist: '" + path + "'")
        return int(output[0])

    # -----------------------------------------------------------------

    def file_nbits(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return self.file_nbytes(path) * 8

    # -----------------------------------------------------------------

    def to_home_directory(self):

        """
        This function ...
        """

        # Navigate to the home directory
        self.execute("cd ~", output=False)

    # -----------------------------------------------------------------

    def create_directory(self, path, show_output=False, recursive=False):

        """
        This function ...
        :param path:
        :param show_output:
        :param recursive:
        :return:
        """

        # Create the remote directory
        if recursive: command = "mkdir -p '" + path + "'"
        else: command = "mkdir '" + path + "'"

        # Execute
        output = self.execute(command, output=True, show_output=show_output)
        for line in output:
            if "cannot create directory" in line: raise IOError("Cannot create directory '" + path + "'")
            elif "quota exceeded" in line: raise IOError("Cannot create directory '" + path + "': disk quota exceeded")

    # -----------------------------------------------------------------

    def create_directories(self, *paths, **kwargs):

        """
        This function ...
        :param paths:
        :param kwargs:
        :return:
        """

        recursive = kwargs.pop("recursive", False)

        # Create the remote directories
        if recursive: command = "mkdir -p '" + "' '".join(paths) + "'"
        else: command = "mkdir '" + "' '".join(paths) + "'"

        # Execute
        output = self.execute(command, output=False)
        for line in output:
            if "cannot create directory" in line:
                which = line.split("cannot create directory ")[1].split(":")[0]
                raise IOError("Cannot create directory " + which)
            elif "quota exceeded" in line: raise IOError("Cannot create directories: disk quota exceeded")

    # -----------------------------------------------------------------

    def create_directory_in(self, base_path, name, show_output=False, recursive=False):

        """
        This function ...
        :param base_path:
        :param name:
        :param show_output:
        :param recursive:
        :return:
        """

        directory_path = fs.join(base_path, name)
        self.create_directory(directory_path, show_output=show_output, recursive=recursive)
        return directory_path

    # -----------------------------------------------------------------

    def decompress_directory_in_place(self, filepath, remove=False, show_output=False):

        """
        This function ...
        :param filepath:
        :param remove:
        :param show_output:
        :return:
        """

        # Inform the user
        self.info("Decompressing '" + filepath + "' ...")

        # tar.gz file
        if filepath.endswith(".tar.gz"):

            # Determine the path of the directory
            new_path = filepath.split(".tar.gz")[0]

            # fs.create_directory(new_path)
            dir_path = fs.directory_of(new_path)

            # Debugging
            self.debug("New path: '" + new_path + "'")
            self.debug("Decompressing in directory '" + dir_path + "' ...")

            dir_path = fs.join(dir_path, filepath.split(".tar.gz")[0])
            self.create_directory(dir_path)

            # Decompress
            #command = "tar -zxvf " + filepath + " --directory " + dir_path
            command = "tar -zxvf " + filepath + " -C " + dir_path + " --strip-components=1"
            log.debug("Decompress command: '" + command + "'")

            # Execute the command
            self.execute(command, show_output=show_output)

        # Not implemented
        else: raise NotImplementedError("Not implemented yet")

        # Remove the file
        if remove: self.remove_file(filepath, show_output=show_output)

        # Return the new path
        return new_path

    # -----------------------------------------------------------------

    def decompress_directory_to(self, filepath, destination, remove=False, show_output=False):

        """
        This function ...
        :param filepath:
        :param destination:
        :param remove:
        :param show_output:
        :return:
        """

        # Inform the user
        self.info("Decompressing '" + filepath + "' ...")

        # tar.gz file
        if filepath.endswith(".tar.gz"):

            # Debugging
            self.debug("Decompressing to '" + destination + "' ...")

            # Decompress
            command = "tar -zxvf " + filepath + " -C " + destination + " --strip-components=1"
            log.debug("Decompress command: '" + command + "'")

            # Execute the command
            self.execute(command, show_output=show_output)

        # Not implemented
        else: raise NotImplementedError("Not implemented yet")

        # Remove the file, if requested
        if remove: self.remove_file(filepath, show_output=show_output)

    # -----------------------------------------------------------------

    def decompress_file(self, path, new_path):

        """
        This function ...
        :param path:
        :param new_path:
        :return:
        """

        if path.endswith(".bz2"): self.decompress_bz2(path, new_path)
        elif path.endswith(".gz"): self.decompress_gz(path, new_path)
        elif path.endswith(".xz"): self.decompress_xz(path, new_path)
        elif path.endswith(".zip"): self.decompress_zip(path, new_path)
        else: raise ValueError("Unrecognized archive type (must be bz2, gz or zip)")

    # -----------------------------------------------------------------

    def decompress_bz2(self, path, new_path):

        """
        This function ...
        """

        raise NotImplementedError()

    # -----------------------------------------------------------------

    def decompress_gz(self, path, new_path):

        """
        This function ...
        :param path:
        :param new_path:
        :return:
        """

        command = "tar -zxvf " + path + " -C " + new_path
        self.execute(command, show_output=True)

    # -----------------------------------------------------------------

    def decompress_xz(self, path, new_path):

        """
        This function ...
        :param path:
        :param new_path:
        :return:
        """

        command = "tar -xvf " + path + " -C " + new_path
        self.execute(command, show_output=True)

    # -----------------------------------------------------------------

    def decompress_zip(self, path, new_path):

        """
        This function ...
        :param path:
        :param new_path:
        :return:
        """

        command = "unzip " + path + " -d " + new_path
        self.execute(command, show_output=True)

    # -----------------------------------------------------------------

    def write_line(self, filepath, line):

        """
        This function ...
        :param filepath:
        :param line:
        :return:
        """

        #if "'" in line: command = 'echo "' + line + '" > ' + filepath
        #else: command = "echo '" + line + "' > " + filepath

        if strings.contains_single_quotes(line): command = 'echo "' + line.replace("$", "\$").replace('"', r'\"') + '" > ' + filepath
        else: command = "echo '" + line + "' > " + filepath

        output = self.execute(command)

        for line in output:
            if "disk quota exceeded" in line.lower(): raise IOError("Cannot write to the file '" + filepath + "': disk quota exceeded")

    # -----------------------------------------------------------------

    def write_lines(self, filepath, lines):

        """
        This function ...
        :param filepath:
        :param lines:
        :return:
        """

        self.write_line(filepath, lines[0])
        for line in lines[1:]: self.append_line(filepath, line)

    # -----------------------------------------------------------------

    def append_line(self, filepath, line):

        """
        This function ...
        :param filepath:
        :param line:
        :return:
        """

        if "'" in line: command = 'echo "' + line + '" >> ' + filepath
        else: command = "echo '" + line + "' >> " + filepath
        self.execute(command, output=False)

    # -----------------------------------------------------------------

    def append_lines(self, filepath, lines):

        """
        This function ...
        :param filepath:
        :param lines:
        :return:
        """

        for line in lines: self.append_line(filepath, line)

    # -----------------------------------------------------------------

    def read_lines(self, path, add_sep=False):

        """
        This function ...
        :param path:
        :param add_sep:
        :return:
        """

        # Make path absolute
        path = self.absolute_path(path)

        #command = "cat '" + path + "' | while read CMD; do     echo $CMD; done" # DOES WEIRD THINGS
        command = 'while IFS= read -r LINE; do     echo "$LINE"; done < "' + path + '"'

        # 2 tries
        try: output = self.execute(command)
        except RuntimeError:
            log.warning("A runtime error was encountered. Trying again ...")
            output = self.execute(command)

        # Loop over the output lines
        for line in output:
            if add_sep: yield line + "\n"
            else: yield line

    # -----------------------------------------------------------------

    def filter_lines(self, path, string, full=False):

        """
        This function ...
        :param path:
        :param string:
        :param full:
        :return:
        """

        # Make path absolute
        path = self.absolute_path(path)

        # Create the command
        if full: command = 'grep -F "' + string + '" "' + path + '" -x --color=never'
        else: command = 'grep -F "' + string + '" "' + path + '" --color=never'

        # Get the output of the grep
        output = self.execute(command)

        # Return the lines
        return output

    # -----------------------------------------------------------------

    def contains_line(self, path, string, full=False):

        """
        This function ...
        :param path:
        :param string:
        :param full:
        :return:
        """

        # Get the lines that match
        lines = self.filter_lines(path, string, full=full)

        # Matching lines?
        return len(lines) > 0

    # -----------------------------------------------------------------

    def get_last_line_containing(self, path, string):

        """
        This function ...
        :param path:
        :param string:
        :return:
        """

        # Get the matching lines
        lines = self.filter_lines(path, string, full=False)

        # Return the last line
        if len(lines) > 0: return lines[-1]
        else: return None

    # -----------------------------------------------------------------

    @memoize_method
    def get_last_line_containing_memoize(self, path, string):

        """
        This function ...
        :param path:
        :param string:
        :return:
        """

        return self.get_last_line_containing(path, string)

    # -----------------------------------------------------------------

    def get_first_line_containing(self, path, string):

        """
        This function ...
        :param path:
        :param string:
        :return:
        """

        # Get the matching lines
        lines = self.filter_lines(path, string, full=False)

        # Return the first line
        if len(lines) > 0: return lines[0]
        else: return None

    # -----------------------------------------------------------------

    @memoize_method
    def get_first_line_containing_memoize(self, path, string):

        """
        This function ...
        :param path:
        :param string:
        :return:
        """

        return self.get_first_line_containing(path, string)

    # -----------------------------------------------------------------

    def get_lines_between_first_and_last_occurence(self, path, string):

        """
        This function ...
        :param path:
        :param string:
        :return:
        """

        # COMMANDS: WORKS!
        # read first last <<< $(grep -n "/home/anersesi/SKIRT/run/2018-01-13--12-20-31-140924/M83.ski" screenlog.0 | awk -F: 'NR==1 {printf "%d ", $1}; END{print $1}')
        # awk -v f=$first -v l=$last 'NR>=f && NR<=l' screenlog.0

        pass

    # -----------------------------------------------------------------

    def get_lines_between(self, path, start, end):

        """
        This function ...
        :param path:
        :param start:
        :param end:
        :return:
        """

        # awk '/140924\/M83.ski/,/Finished simulation/' screenlog.0

        # Check
        if "'" in start: raise ValueError("Patterns cannot contain single quotes")
        if "'" in end: raise ValueError("Patterns cannot contain single quotes")

        # Set command
        command = "awk '/" + start.replace("/", "\/") + "/,/" + end.replace("/", "\/") + "/' " + path

        # Get the output
        output = self.execute(command)

        # Return the output
        return output

    # -----------------------------------------------------------------

    def get_noccurences(self, path, string):

        """
        Thisn function ...
        :param path:
        :param string:
        :return:
        """

        count = 0
        for line in self.read_lines(path):
            if string in line: count += 1
        return count

    # -----------------------------------------------------------------

    def get_lines(self, path, add_sep=False):

        """
        This function ...
        :param path:
        :param add_sep:
        :return:
        """

        return list(self.read_lines(path, add_sep=add_sep))

    # -----------------------------------------------------------------

    def get_text(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return "\n".join(self.get_lines(path))

    # -----------------------------------------------------------------

    def read_first_lines(self, path, nlines):

        """
        This function ...
        :param path:
        :param nlines:
        :return:
        """

        command = "head -n " + str(nlines) + " " + path

        # Execute
        output = self.execute(command)

        # Return the lines
        return output

    # -----------------------------------------------------------------

    def read_first_line(self, path):

        """
        Thisf unction ...
        :param path:
        :return:
        """

        return self.read_first_lines(path, 1)[0]

    # -----------------------------------------------------------------

    def get_first_line(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return self.read_first_line(path)

    # -----------------------------------------------------------------

    def read_last_lines(self, path, nlines):

        """
        This function ...
        :param path:
        :param nlines
        :return:
        """

        # Execute 'tail'
        output = self.execute("tail -" + str(nlines) + " " + path)

        # Return the lines
        return output

    # -----------------------------------------------------------------

    def read_last_line(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return self.read_last_lines(path, 1)[0]

    # -----------------------------------------------------------------

    def get_last_line(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return self.read_last_line(path)

    # -----------------------------------------------------------------

    def read_lines_reversed(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Expand the path to absolute form
        path = self.absolute_path(path)

        # Launch the command
        command = "sed '1!G;h;$!d' " + path
        output = self.execute(command)

        # Return the output
        for line in output: yield line

    # -----------------------------------------------------------------

    def download_retry(self, origin, destination, timeout=None, new_name=None, compress=False, show_output=False, connect_timeout=90, max_nattempts=3):

        """
        This function ...
        :return:
        """

        success = False
        nattempts = 0

        # Try multiply times
        while not success:

            # Check the number of attempts
            if nattempts == max_nattempts: break

            # Try uploading
            success = self.download(origin, destination, timeout, new_name, compress, show_output, connect_timeout, return_success=True)

            # Increment the number of attempts
            nattempts += 1

        # Return if it was succesful
        return success

    # -----------------------------------------------------------------

    def download(self, origin, destination, timeout=None, new_name=None, compress=False, show_output=False,
                 connect_timeout=90, return_success=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :param connect_timeout:
        :param return_success:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug: show_output = True

        # WITH SPACES:
        # scp documents.zip remote:"\"/var/path/containing some spaces/foo/\""

        # Construct the command string
        copy_command = "scp "
        if self.host.port is not None: copy_command += "-P " + str(self.host.port) + " "
        if connect_timeout is not None: copy_command += "-o ConnectTimeout=" + str(connect_timeout) + " "
        if compress: copy_command += "-C "

        # Add the host address
        #copy_command += self.host.user + "@" + self.host.name + ":"

        origin_type = None

        # If the origin is a string, we assume it represents a single file path or directory path
        if types.is_string_type(origin):

            # Check if the origin represents a file
            #if self.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if self.is_file(origin):

                copy_command += self.host.user + "@" + self.host.name + ":"

                origin_type = "file"
                copy_command += "'" + origin.replace(" ", "\ ") + "' "

            # Check if it represents a directory
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\ ") + "/* " + "-r "
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\\ ") + "/* "
            elif self.is_directory(origin):

                origin_type = "directory"
                #copy_command += "-r '" + origin + "' "
                #copy_command += "'" + origin.replace(" ", "\ ") + "/*' "

                copy_command += "-r "
                copy_command += self.host.user + "@" + self.host.name + ":"
                copy_command += "'" + origin.replace(" ", "\ ") + "/*' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file on the remote host")

        # If the origin is a list, we assume it contains multiple file paths
        elif types.is_sequence(origin):

            # Check whether the files exist remotely
            for file_path in origin:
                if not self.is_file(file_path): raise ValueError("The file " + file_path + " does not exist on the remote host")

            origin_type = "files"

            copy_command += self.host.user + "@" + self.host.name + ":"

            # Escape possible space characters
            #origin = [path.replace(" ", "\\\ ") for path in origin]
            origin_strings = ["'" + path.replace(" ", "\ ") + "'" for path in origin]

            # Add a quotation mark character because the seperate file paths are going to be separated by spaces
            # (the command is going to be of the form scp username@ip.of.server.copyfrom:"file1.log file2.log" "~/yourpathtocopy")
            copy_command += '"'

            # Add the file paths to the command string
            copy_command += " ".join(origin_strings)

            # Add another quotation mark to identify the end of the filepath list
            copy_command += '" '

        # Invalid
        else: raise ValueError("Invalid origin: " + str(origin))

        # Check whether the destination is a directory
        if "." in fs.name(destination): destination_string = "'" + destination + "'" # destination is a file
        elif not fs.is_directory(destination): # destination should be existing directory
            raise ValueError("Destination directory '" + destination + "' does not exist")
        else:
            if new_name is not None: destination_string = "'" + destination + "/" # existing directory
            else: destination_string = "'" + destination + "/'" # existing directory

        copy_command += destination_string
        if new_name is not None: copy_command += new_name + "'"

        # Debugging
        self.debug("Download command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            index = child.expect(['password: ', pexpect.EOF])
            if index == 0: child.sendline(self.host.password)
            else:
                if return_success: return False
                else:
                    lines = child.before.split("\r\n")
                    raise RuntimeError(" ".join(lines))

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        if show_output: child.logfile = sys.stdout

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        else:

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # Execute the command and get the output
        child.expect(pexpect.EOF, timeout=None)
        lines = child.before.split("\r\n") # output lines
        child.close()

        #
        if not show_output:

            # Retrieve file contents -- this will be
            # 'First line.\nSecond line.\n'
            #stdout = child.logfile.getvalue()

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(" ".join(lines))

            # Get the output lines
            #lines = stdout.split("\n")

            # Debugging: show the output of the scp command
            self.debug("Copy stdout: " + str(" ".join(lines)))

        # Check for errors
        for line in lines:
            if "not a regular file" in line: raise ValueError(line)
            if "scp: ambiguous target" in line: raise ValueError(line)
            if "No such file or directory" in line: raise ValueError(line)
            if "No space left on device" in line: raise IOError("Not enough disk space")

        #else:

        # EXTRA CHECK
        if origin_type == "files":

            for path in origin:
                name = fs.name(path)
                local_path = fs.join(destination, name)
                if not fs.is_file(local_path): raise RuntimeError("Something went wrong: file '" + name + "' is missing from destination (" + local_path + ")")

        elif origin_type == "file":

            if fs.is_directory(destination):

                name = fs.name(origin) if new_name is None else new_name
                local_path = fs.join(destination, name)
                if not fs.is_file(local_path): raise RuntimeError("Something went wrong: file '" + name + "' is missing from destination (" + local_path + ")")

            # It must be a file then
            elif not fs.is_file(destination): raise RuntimeError("Something went wrong: file '" + destination + "' is missing")

        elif origin_type == "directory":

            for remote_path in self.files_in_path(origin):
                filename = fs.name(remote_path)
                local_path = fs.join(destination, filename)
                if not fs.is_file(local_path): raise RuntimeError("Something went wrong: file '" + filename + "' is missing (" + local_path + ")")

            for remote_path in self.directories_in_path(origin):
                dirname = fs.name(remote_path)
                local_path = fs.join(destination, dirname)
                if not fs.is_directory(local_path): raise RuntimeError("Something went wrong: directory '" + dirname + "' is missing (" + local_path + ")")

        else: raise ValueError("Invalid origin type")

        # Success: return True
        return True

    # -----------------------------------------------------------------

    def download_file_to(self, filepath, destination, remove=False, new_name=None):

        """
        This function ...
        :param filepath:
        :param destination:
        :param remove:
        :param new_name:
        :return:
        """

        # Determine the local path
        filename = new_name if new_name is not None else fs.name(filepath)
        local_path = fs.join(destination, filename)

        # Download
        self.download(filepath, local_path)

        # Check whether present
        if not fs.is_file(local_path): raise RuntimeError("Something went wrong downloading the file")

        # Now we can safely remove the remote file if requested
        if remove: self.remove_file(filepath)

        # Return the local path
        return local_path

    # -----------------------------------------------------------------

    def download_directory_to(self, dirpath, destination, remove=False, new_name=None):

        """
        This function ...
        :param dirpath:
        :param destination:
        :param remove:
        :param new_name:
        :return:
        """

        # Determine the local path
        dirname = new_name if new_name is not None else fs.name(dirpath)
        local_path = fs.join(destination, dirname)

        # Download
        self.download(dirpath, destination)

        # Check whether present
        if not fs.is_directory(local_path): raise RuntimeError("Something went wrong downloading the directory")

        # Now we can safely remove the remote directory if requested
        if remove: self.remove_directory(dirpath)

        # Return the local path
        return local_path

    # -----------------------------------------------------------------

    def upload_retry(self, origin, destination, timeout=None, new_name=None, compress=False, show_output=False, connect_timeout=90, max_nattempts=3):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :param connect_timeout:
        :param max_nattempts
        :return:
        """

        success = False
        nattempts = 0

        # Try multiply times
        while not success:

            # Check the number of attempts
            if nattempts == max_nattempts: break

            # Try uploading
            success = self.upload(origin, destination, timeout, new_name, compress, show_output, connect_timeout, return_success=True)

            # Increment the number of attempts
            nattempts += 1

        # Return if it was succesful
        return success

    # -----------------------------------------------------------------

    def upload(self, origin, destination, timeout=None, new_name=None, compress=False, show_output=False,
               connect_timeout=90, return_success=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :param connect_timeout:
        :param return_success:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug: show_output = True

        # WITH SPACES:
        # scp documents.zip remote:"\"/var/path/containing some spaces/foo/\""

        # Construct the command string
        copy_command = "scp "
        if self.host.port is not None: copy_command += "-P " + str(self.host.port) + " "
        if connect_timeout is not None: copy_command += " -o ConnectTimeout=" + str(connect_timeout) + " "
        if compress: copy_command += "-C "

        origin_type = None

        # If the origin is a string, we assume it represents a single file path or directory path
        if types.is_string_type(origin):

            # Check if the origin represents a file
            #if fs.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if fs.is_file(origin):
                origin_type = "file"
                #copy_command += "'" + origin.replace(" ", "\ ") + "' "
                #copy_command += "'\\\"" + origin + "\\\"' "
                copy_command += "'" + origin + "' "

            # Check if it represents a directory
            #elif fs.is_directory(origin): copy_command += "-r " + origin.replace(" ", "\\\ ") + "/ "
            elif fs.is_directory(origin):
                origin_type = "directory"
                #copy_command += "-r '" + origin.replace(" ", "\ ") + "/' "
                #copy_command += "-r '\\\"" + origin + "/\\\"' "
                copy_command += "-r '" + origin + "' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file")

        # If the origin is a list, we assume it contains multiple file paths
        elif types.is_sequence(origin):

            # Check whether the files exist locally
            for file_path in origin:
                if not fs.is_file(file_path): raise ValueError("The file " + file_path + " does not exist")

            #origin = ["'" + path.replace(" ", "\ ") + "'" for path in origin]
            #origin = ["'\\\"" + path + "\\\"'" for path in origin]
            #origin = [path for path in origin]
            #origin_strings = ["'" + path.replace(" ", "\ ") + "'" for path in origin]
            origin_strings = ["'" + path + "'" for path in origin]

            origin_type = "files"

            # Add the file paths to the command string
            copy_command += " ".join(origin_strings) + " "

        # Invalid argument
        else: raise ValueError("The origin must be a string or a list of strings")

        # Set destination string depending on whether it is a file or a directory
        if "." in fs.name(destination):
            destination_string = "'" + destination.replace(" ", "\ ") + "'" # destination is a file
        elif not self.is_directory(destination): # destination should be existing directory
            raise ValueError("Destination directory '" + destination + "' does not exist on remote host '" + self.host_id + "'")
        else:
            if new_name is not None: destination_string = "'" + destination.replace(" ", "\ ") + "/" # existing directory
            else: destination_string = "'" + destination.replace(" ", "\ ") + "/'" # existing directory

        # Construct command
        copy_command += self.host.user + "@" + self.host.name + ":" + destination_string
        if new_name is not None: copy_command += new_name + "'"

        # Debugging
        self.debug("Upload command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            index = child.expect(['password: ', pexpect.EOF])
            if index == 0: child.sendline(self.host.password)
            else:
                if return_success: return False
                else:
                    lines = child.before.split("\r\n")
                    raise RuntimeError(" ".join(lines))

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        if show_output: child.logfile = sys.stdout

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        else:

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # Execute the command and get the output
        try:
            child.expect(pexpect.EOF, timeout=None)
            lines = child.before.split("\r\n")
        except pexpect.EOF:
            pass
            lines = None
        child.close()

        if lines is None: lines = child.logfile.getvalue()

        # Show output lines in debug mode
        if not show_output:

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(" ".join(lines))

            # Debugging: show the output of the scp command
            self.debug("Copy stdout: " + str(" ".join(lines)))

        # Check for errors
        for line in lines:
            if "not a regular file" in line: raise ValueError(line)
            if "scp: ambiguous target" in line: raise ValueError(line)
            if "No such file or directory" in line: raise ValueError(line)
            if "No space left on device" in line: raise IOError("Not enough disk space")

        #else:
            #print(sys.stdout)
            #index = 0
            #for line in sys.stdout: print(line)
            #for line in fs.reverse_read_line_impl(sys.stdout):
            #    if index == 4: break
            #    print(line)
            #    index += 1

        # EXTRA check
        # Check whether remote files are present
        if origin_type == "files":

            for filepath in origin:
                filename = fs.name(filepath)
                remote_path = fs.join(destination, filename)
                if not self.is_file(remote_path): raise RuntimeError("Something went wrong: file " + filename + " is not present at destination (" + remote_path + ")")

        elif origin_type == "file":

            if self.is_directory(destination):

                filename = fs.name(origin) if new_name is None else new_name
                remote_path = fs.join(destination, filename)
                if not self.is_file(remote_path): raise RuntimeError("Something went wrong: file " + filename + " is not present at destination (" + remote_path + ")")

            elif not self.is_file(destination): raise RuntimeError("Something went wrong: file " + destination + " is not present")

        elif origin_type == "directory":

            directory_name = fs.name(origin)

            for local_path in fs.files_in_path(origin):

                filename = fs.name(local_path)
                #remote_path = fs.join(destination, filename)
                remote_path = fs.join(destination, directory_name, filename)
                if not self.is_file(remote_path): raise RuntimeError("Something went wrong: file " + filename + " is not present at destination (" + remote_path + ")")

            for local_path in fs.directories_in_path(origin):

                dirname = fs.name(local_path)
                remote_path = fs.join(destination, directory_name, dirname)
                if not self.is_directory(remote_path): raise RuntimeError("Something went wrong: directory " + dirname + " is not present at destination (" + remote_path + ")")

        else: raise ValueError("Invalid origin type: " + str(origin_type))

        # Success: return True
        return True

    # -----------------------------------------------------------------

    def upload_file_to(self, filepath, destination, remove=False, new_name=None, show_output=False, replace=True):

        """
        This function ...
        :param filepath:
        :param destination:
        :param remove:
        :param new_name:
        :param show_output:
        :param replace:
        :return:
        """

        # Determine the remote file path
        filename = new_name if new_name is not None else fs.name(filepath)
        remote_path = fs.join(destination, filename)

        # Already a file on the remote?
        if self.is_file(remote_path):
            if replace: self.remove_file(remote_path)
            else: raise IOError("Already a file on the remote: '" + remote_path + "'")

        # Upload the file
        self.upload(filepath, remote_path, show_output=show_output)

        # Check whether present
        if not self.is_file(remote_path): raise RuntimeError("Something went wrong uploading the file")

        # Now we can safely remove the local file if requested
        if remove: fs.remove_file(filepath)

        # Return the remote file path
        return remote_path

    # -----------------------------------------------------------------

    def upload_files_to(self, filepaths, destination, remove=False, show_output=False, replace=True):

        """
        This function ...
        :param filepaths:
        :param destination:
        :param remove:
        :param show_output:
        :param replace:
        :return:
        """

        # Initialize a list for the remote filepaths
        remote_paths = []

        # Upload each file
        for filepath in filepaths: remote_paths.append(self.upload_file_to(filepath, destination, remove=remove, show_output=show_output, replace=replace))

        # Return the list of remote filepaths
        return remote_paths

    # -----------------------------------------------------------------

    def upload_files_in_path_to(self, path, destination, remove=False, show_output=False, replace=True, **kwargs):

        """
        This function ...
        :param path:
        :param destination:
        :param remove:
        :param show_output:
        :param replace:
        :param kwargs:
        :return:
        """

        # Upload
        return self.upload_files_to(fs.files_in_path(path, **kwargs), destination, remove=remove, show_output=show_output, replace=replace)

    # -----------------------------------------------------------------

    def upload_directory_to(self, dirpath, destination, remove=False, compress=False, show_output=False):

        """
        This function ...
        :param dirpath:
        :param destination:
        :param remove:
        :param compress:
        :param show_output:
        :return:
        """

        # Determine remote directory name
        remote_path = fs.join(destination, fs.name(dirpath))
        self.upload(dirpath, destination, compress=compress, show_output=show_output)

        # Check whether present
        if not self.is_directory(remote_path): raise RuntimeError("Something went wrong uploading the directory")

        # Now we can safely remove the local directory if requested
        if remove: fs.remove_directory(dirpath)

        # Return the remote directory path
        return remote_path

    # -----------------------------------------------------------------

    def synchronize(self, origin, destination, timeout=None, connect_timeout=90, compress=True, show_output=False,
                    delete=False, hidden=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param connect_timeout:
        :param compress:
        :param show_output:
        :param delete:
        :param hidden: also synchronize hidden files
        :return:
        """

        # Check whether connection timeout can be specified
        if connect_timeout is not None:
            from ..tools import terminal
            lines = terminal.execute("rsync --contimeout", return_first=True)
            for line in lines:
                if "unknown option" in line:
                    log.warning("The connection timeout cannot be specified with your version of rsync")
                    connect_timeout = None
                    break

        # Create absolute paths
        origin = fs.absolute_path(origin)

        # Set full path to the destination
        destination = self.absolute_path(destination)

        # Check paths
        if not fs.is_directory(origin): raise IOError("Origin direcory '" + origin + "' does not exist")
        if not self.is_directory(destination): raise IOError("Destination directory '" + destination + "' does not exist")

        # If debugging is enabled, always show the scp output
        if log.is_debug: show_output = True

        # Construct the command string
        command = "rsync -chavP --stats " #--stats user@remote.host:/path/to/copy /path/to/local/storage"
        if compress: command += "-z "
        if connect_timeout is not None: command += "--contimeout=" + str(connect_timeout) + " "
        if delete: command += "--delete "
        if not hidden:
            command += '--exclude ".*/" ' # hidden directories
            command += '--exclude=".*" '  # hidden files

        # Add the origin directory, adding '/' to synchronize its contents
        command += "'" + origin + "/' "

        # Set destination string depending on whether it is a file or a directory
        if "." in fs.name(destination): destination_string = "'" + destination.replace(" ", "\ ") + "'"  # destination is a file
        elif not self.is_directory(destination):  # destination should be existing directory
            raise ValueError("Destination directory '" + destination + "' does not exist on remote host '" + self.host_id + "'")
        else: destination_string = "'" + destination.replace(" ", "\ ") + "/'"  # existing directory

        # Construct command
        command += self.host.user + "@" + self.host.name + ":" + destination_string

        # Debugging
        self.debug("Upload command: " + command)

        # Create the pexpect child instance
        child = pexpect.spawn(command, timeout=timeout)
        if self.host.password is not None:
            index = child.expect(['password: ', pexpect.EOF])
            if index == 0: child.sendline(self.host.password)
            else:
                lines = child.before.split("\r\n")
                raise RuntimeError(" ".join(lines))

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        if show_output: child.logfile = sys.stdout

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        else:

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # Execute the command and get the output
        try:
            child.expect(pexpect.EOF, timeout=None)
            lines = child.before.split("\r\n")
        except pexpect.EOF:
            pass
            lines = None
        child.close()

        # Get output lines
        if lines is None: lines = child.logfile.getvalue()

        # Show output lines in debug mode
        if not show_output:

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(" ".join(lines))

            # Debugging: show the output of the scp command
            self.debug("Synchronize stdout: " + str(" ".join(lines)))

        # Check for errors
        for line in lines:
            # TODO: fill in the error message strings for rsync
            #if "not a regular file" in line: raise ValueError(line)
            #if "scp: ambiguous target" in line: raise ValueError(line)
            #if "No such file or directory" in line: raise ValueError(line)
            #if "No space left on device" in line: raise IOError("Not enough disk space")
            pass

        # Return success
        return True

    # -----------------------------------------------------------------

    def version_of(self, name_or_path, show_output=False):

        """
        This function ...
        :param name_or_path:
        :param show_output:
        :return:
        """

        # Execute
        output = self.execute(name_or_path + " --version", show_output=show_output)

        # Get name
        name = fs.name(name_or_path)

        # Return the line with the version number
        for line in output:
            if line.startswith("[warn]"): continue
            if "warning: setlocale:" in line: continue
            if line.startswith(name): return line.split(name)[1].strip()
            else: return line.strip()

        # Cannot determine version
        return None

    # -----------------------------------------------------------------

    def find_conda(self):

        """
        This function ...
        :return:
        """

        # Find conda path
        if self.is_executable("conda"): conda_executable_path = self.find_executable("conda")
        else: conda_executable_path = None

        # Find conda installation in the home directory
        if conda_executable_path is None:

            conda_path = fs.join(self.home_directory, "miniconda", "bin", "conda")
            if self.is_file(conda_path): conda_executable_path = conda_path

            # Search in scratch path
            if conda_executable_path is None and self.scratch_path is not None:

                conda_path = fs.join(self.scratch_path, "miniconda", "bin", "conda")
                if self.is_file(conda_path): conda_executable_path = conda_path

        # If conda is present
        if conda_executable_path is not None:

            conda_installation_path = fs.directory_of(fs.directory_of(conda_executable_path))
            conda_main_executable_path = conda_executable_path

            # Return the paths
            return conda_installation_path, conda_main_executable_path

        # Conda not present
        else: return None, None

    # -----------------------------------------------------------------

    @property
    def has_conda(self):

        """
        This function ...
        :return:
        """

        return self.is_executable("conda")

    # -----------------------------------------------------------------

    @property
    def conda_version(self):

        """
        This function ...
        :return:
        """

        if not self.has_conda: return None

        output = self.execute("conda -V")
        return output[0]

    # -----------------------------------------------------------------

    def conda_version_at(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        output = self.execute(path + " -V")
        return output[0]

    # -----------------------------------------------------------------

    @property
    def conda_installation_path(self):

        """
        This function ...
        :return:
        """

        conda_exec_path = self.find_executable("conda")

        if conda_exec_path is None:
            try:
                conda_installation_path = self.find_directory_in_path(self.home_directory, contains="conda")
                conda_exec_path = fs.join(conda_installation_path, "bin", "conda")
            except ValueError: pass
        if conda_exec_path is None: conda_exec_path = self.find_file_in_path(self.home_directory, recursive=True, exact_name="conda")

        if fs.name(fs.directory_of(fs.directory_of(fs.directory_of(conda_exec_path)))) == "envs":
            return fs.directory_of(fs.directory_of(fs.directory_of(fs.directory_of(conda_exec_path))))
        else: return fs.directory_of(fs.directory_of(conda_exec_path))

    # -----------------------------------------------------------------

    @property
    def conda_main_executable_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.conda_installation_path, "bin", "conda")
        if fs.is_file(path): return path
        else: return None

    # -----------------------------------------------------------------

    def is_conda_environment(self, name, conda_path="conda"):

        """
        This function ...
        :param name:
        :param conda_path:
        :return:
        """

        return name in self.conda_environments(conda_path=conda_path)

    # -----------------------------------------------------------------

    def conda_environments(self, conda_path="conda"):

        """
        This function ...
        :param conda_path:
        :return:
        """

        if not self.is_executable(conda_path): raise ValueError("The conda executable does not exist: '" + conda_path + "'")

        output = self.execute(conda_path + " env list")
        envs = []
        for line in output:
            if line.startswith("#"): continue
            if line.strip() == "": continue
            envs.append(line.split()[0])
        return envs

    # -----------------------------------------------------------------

    @property
    def conda_environment_for_pts(self):

        """
        This function ...
        :return:
        """

        # Get environment name
        pts_alias = self.resolve_alias("pts")
        pts_python_path = pts_alias.split()[0]

        # If just 'python'
        #if pts_python_path == "python": raise Exception("Cannot determine the conda environment used for pts")
        if pts_python_path == "python": return None

        # Get the environment name
        env_path = fs.directory_of(fs.directory_of(pts_python_path))
        environment_name = fs.name(env_path)
        return environment_name

    # -----------------------------------------------------------------

    @property
    def conda_pts_environment_bin_path(self):

        """
        This function ...
        :return: 
        """

        # Find conda
        conda_installation_path, conda_main_executable_path = self.find_conda()

        # Get env name
        env_name = self.conda_environment_for_pts

        # Set environment bin path
        environment_bin_path = fs.join(conda_installation_path, "envs", env_name, "bin")

        # Return the bin path
        return environment_bin_path

    # -----------------------------------------------------------------

    @property
    def conda_pts_environment_python_path(self):

        """
        This function ...
        :return: 
        """

        # Get env bin path
        bin_path = self.conda_pts_environment_bin_path

        # Set python path
        conda_python_path = fs.join(bin_path, "python")

        # Return the path
        return conda_python_path

    # -----------------------------------------------------------------

    @property
    def conda_pip_name_for_pts(self):

        """
        This function ...
        :return:
        """

        aliases = self.aliases

        #pip_path = None
        name = None
        for alias in aliases:

            target = aliases[alias]

            if not target.endswith("pip"): continue

            #pip_path = None
            name = alias

        return name

    # -----------------------------------------------------------------

    @property
    def conda_jupyter_name_for_pts(self):

        """
        This function ...
        :return:
        """

        aliases = self.aliases

        #jupyter_path = None
        name = None
        for alias in aliases:

            target = aliases[alias]

            if not target.endswith("jupyter"): continue

            name = alias

        # Return the name
        return name

    # -----------------------------------------------------------------

    def conda_active_environment(self, conda_path="conda"):

        """
        This function ...
        :param conda_path:
        :return:
        """

        if not self.is_executable(conda_path): raise ValueError("The conda executable does not exist: '" + conda_path + "'")
        output = self.execute(conda_path + " env list")
        env = None
        for line in output:
            if line.startswith("#"): continue
            if line.strip() == "": continue
            splitted = line.split()
            if splitted[1].strip() == "*":
                env = splitted[0].strip()
                break
        return env

    # -----------------------------------------------------------------

    def conda_active_environment_fast(self, assert_activated=False):

        """
        This function ...
        :param assert_activated
        :return:
        """

        # Execute dummy command
        #self.execute("echo")
        # Send the command
        #self.ssh.sendline("echo")
        #self.ssh.prompt()

        # Get last line
        last = self.ssh.before.split("\r\n")[-1].strip()
        if last.startswith("(") and last.endswith(")"): return last[1:-1]
        elif assert_activated: raise RuntimeError("Conda environment cannot be determined from '" + last + "'")
        else: return None

    # -----------------------------------------------------------------

    def activate_conda_environment(self, environment_name, conda_path="conda", activate_path="activate", show_output=False):

        """
        This function ...
        :param environment_name:
        :param conda_path:
        :param activate_path:
        :param show_output:
        :return:
        """

        previous = self.conda_active_environment(conda_path=conda_path)
        self.execute("source " + activate_path + " " + environment_name, show_output=show_output)
        self.conda_activated = True
        return previous

    # -----------------------------------------------------------------

    def deactivate_conda(self, deactivate_path="deactivate"):

        """
        This function ...
        :param deactivate_path:
        :return:
        """

        self.execute("source " + deactivate_path)
        self.conda_activated = False

    # -----------------------------------------------------------------

    def installed_conda_packages(self, conda_path="conda", environment_name=None):

        """
        This function ...
        :param conda_path:
        :param environment_name:
        :return:
        """

        packages = dict()

        if not self.is_executable(conda_path): raise ValueError("The conda executable does not exist: '" + conda_path + "'")

        # Execute the conda list command
        if environment_name is not None: output = self.execute(conda_path + " list -n " + environment_name)
        else: output = self.execute(conda_path + " list")

        # Loop over the lines
        for line in output:

            #print(line)
            if line.startswith("#"): continue
            if line.strip() == "": continue

            #print(line)

            try:

                name = line.split()[0].strip()
                version = line.split()[1].strip()
                packages[name] = version

            except IndexError:
                log.warning("Unexpected line: '" + line + "'")
                pass

        # return the dictionary
        return packages

    # -----------------------------------------------------------------

    def is_present_conda_package(self, name, environment_name=None, conda_path="conda"):

        """
        This function ...
        :param name:
        :param environment_name:
        :param conda_path:
        :return:
        """

        if not self.is_executable(conda_path): raise ValueError("The conda executable does not exist: '" + conda_path + "'")

        if environment_name is not None: command = conda_path + " list " + name + " --name " + environment_name
        else: command = conda_path + " list " + name

        output = self.execute(command)

        # Find the package
        for line in output:
            if line.startswith("#"): continue
            if line.strip() == "": continue
            if line.split()[0].lower() == name.lower(): return True
        return False

    # -----------------------------------------------------------------

    def available_conda_packages(self, conda_path="conda"):

        """
        This function ...
        :param conda_path:
        :return:
        """

        output = self.execute(conda_path + " search")
        available_packages = []
        for line in output:
            if not line.split(" ")[0]: continue
            available_packages.append(line.split(" ")[0])
        return available_packages

    # -----------------------------------------------------------------

    @property
    def pts_conformity_issues(self):

        """
        This function ...
        :return:
        """

        issues = []

        pts_root_dir_name = fs.name(self.pts_root_path)
        if pts_root_dir_name != "PTS": issues.append("PTS root directory is not called 'PTS'")
        if fs.directory_of(self.pts_root_path) != self.home_directory: issues.append("PTS installation is not located in the home directory")

        if not self.has_conda: issues.append("Conda executable cannot be located based on the PATH")
        else:

            installation_path = self.conda_installation_path
            if fs.name(installation_path) != "miniconda": issues.append("Name of Conda installation directory is not 'miniconda'")
            if fs.directory_of(installation_path) != self.home_directory: issues.append("Conda not installed in home directory")

        if len(issues) == 0: return None
        else: return ", ".join(issues)

    # -----------------------------------------------------------------

    def pts_installation_is_conform(self):

        """
        This function ...
        :return:
        """

        return self.pts_conformity_issues is None

    # -----------------------------------------------------------------

    @property
    def skirt_conformity_issues(self):

        """
        This function ...
        :return:
        """

        issues = []

        skirt_root_dir_name = fs.name(self.skirt_root_path)
        if skirt_root_dir_name != "SKIRT": issues.append("SKIRT root directory is not called 'SKIRT'")

        if fs.directory_of(self.skirt_root_path) != self.home_directory: issues.append("SKIRT installation is not located in the home directory")

        if len(issues) == 0: return None
        else: return ", ".join(issues)

    # -----------------------------------------------------------------

    def skirt_installation_is_conform(self):

        """
        This function ...
        """

        return self.skirt_conformity_issues is None

    # -----------------------------------------------------------------

    @property
    def has_skirt(self):

        """
        This function ...
        :return:
        """

        return self.is_executable("skirt")

    # -----------------------------------------------------------------

    @property
    def skirt_version(self):

        """
        This function ...
        :return:
        """

        # load module if necessary (C libraries etc.)
        if self.has_lmod: self.find_and_load_cpp_compiler()

        output = self.execute("skirt -v")
        skirt_version = None

        for line in output:

            if not "Welcome to" in line: continue
            skirt_version = line.split("Welcome to ")[1].split(" built on")[0] + ")"
            break

        return skirt_version

    # -----------------------------------------------------------------

    @property
    def has_pts(self):

        """
        This function ...
        :return:
        """

        for line in self.execute("pts -v"):

            if "command not found" in line: return False
            elif "No module named pts" in line: return False
            elif "is currently not installed" in line: return False

        return True

    # -----------------------------------------------------------------

    @property
    def pts_version(self):

        """
        This function ...
        :return:
        """

        if self.has_pts:
            output = self.execute("pts --version")
            for line in output:
                if "Welcome to PTS" in line: continue
                return line
            else: return None
        else: return None

    # -----------------------------------------------------------------

    @property
    def system_name(self):

        """
        This function ...
        :return:
        """

        return self.host.system_name

    # -----------------------------------------------------------------

    @property
    def home_directory(self):

        """
        This function ...
        :return:
        """

        # Find out the path to the user's home directory and return it
        return self.resolve_environment_variable("HOME")

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_temp_directories(self):

        """
        This function ...
        :return:
        """

        paths = self.directories_in_path(self.pts_temp_path, returns="path", startswith="session", exact_not_name=self.session_temp_name)
        return list(sorted(paths))

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_temp_names(self):

        """
        This function ...
        :return:
        """

        names = self.directories_in_path(self.pts_temp_path, returns="name", startswith="session", exact_not_name=self.session_temp_name)
        return list(sorted(names))

    # -----------------------------------------------------------------

    @lazyproperty
    def session_temp_directory(self):

        """
        This function ...
        :return:
        """

        # Determine the path to a new temporary directory
        path = fs.join(self.pts_temp_path, time.unique_name("session"))

        # Create the directory using a bash command
        self.create_directory(path)

        # Return the path to the new temporary directory
        return path

    # -----------------------------------------------------------------

    @property
    def session_temp_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.session_temp_directory)

    # -----------------------------------------------------------------

    def new_temp_directory(self):

        """
        This function ...
        :return:
        """

        # Generate the path to a new unique temporary directory
        path = fs.join(self.pts_temp_path, time.unique_name("new"))

        # Create the directory using a bash command
        self.create_directory(path)

        # Return the path to the new temporary directory
        return path

    # -----------------------------------------------------------------

    def absolute_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Replace environment variables
        if "$" in path:

            splitted = path.split("/")
            new_splitted = []
            for part in splitted:
                if part.startswith("$"): new = self.resolve_environment_variable(part[1:])
                else: new = part
                new_splitted.append(new)
            path = fs.join(*new_splitted)

        if path.startswith("~"): return fs.join(self.home_directory, path.split("~/")[1])
        elif path.startswith("/"): return path
        else: return fs.join(self.working_directory, path)

    # -----------------------------------------------------------------

    def absolute_or_in(self, path, in_path):

        """
        This function ...
        :param path:
        :param in_path:
        :return:
        """

        if fs.is_absolute(path): return path
        else: return fs.join(in_path, path)

    # -----------------------------------------------------------------

    def absolute_or_in_home(self, path):
        
        """
        This function ...
        :param path: 
        :return: 
        """

        return self.absolute_or_in(path, self.home_directory)

    # -----------------------------------------------------------------

    def relative_to_home(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return path.split(self.home_directory + "/")[1]

    # -----------------------------------------------------------------

    def relative_to_cwd(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        return path.split(self.working_directory + "/")[1]

    # -----------------------------------------------------------------

    def find_executable(self, name, show_output=False):

        """
        This function ...
        :param name:
        :param show_output:
        :return:
        """

        # Get the output of the 'which' command
        output = self.execute("which " + name, show_output=show_output)

        if len(output) == 0: return None
        else:

            first_line = output[0]

            if "no " + name + " in" in first_line: return None
            else: return first_line

    # -----------------------------------------------------------------

    def is_executable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.find_executable(name) is not None

    # -----------------------------------------------------------------

    @property
    def working_directory(self):

        """
        This function ...
        :return:
        """

        # Find out the path to the current working directory and return it
        #output = self.execute("echo $PWD")
        #return output[0]

        self.ssh.sendline("echo $PWD")
        self.ssh.prompt()
        output = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]
        return output[0]

    # -----------------------------------------------------------------

    @property
    def free_cores(self):

        """
        This function ...
        :return:
        """

        # The calculation doesn't work for multiple nodes
        if self.nodes > 1: raise ValueError("This function cannot be called for systems with multiple nodes")

        # Determine the number of free cores on the only node
        return self.cores_per_node * (1.0 - self.cpu_load)

    # -----------------------------------------------------------------

    @property
    def free_sockets(self):

        """
        This function ...
        :return:
        """

        # The calculation doesn't work for multiple nodes
        if self.nodes > 1: raise ValueError("This function cannot be called for systems with multiple nodes")

        # Determine the number of free cores on the only node
        return self.sockets_per_node * (1.0 - self.cpu_load)

    # -----------------------------------------------------------------

    @property
    def free_memory(self):

        """
        This function uses the 'free' command to get information about the virtual memory usage
        :return:
        """

        # Works on some systems
        output = self.execute("free -t | grep '+ buffers/cache'")

        # There is output
        if len(output) != 0:

            splitted = output[0].split(":")[1].split()
            free = splitted[1]
            # Calculate the free amount of memory in gigabytes
            free = float(free) / 1e6

        else:

            output = self.execute("free -t | grep ^Mem")
            free_without_buffers_cache = float(output[0].split(":")[1].strip().split()[2])

            output = self.execute("free -t | grep 'buffers '")
            if len(output) != 0:

                output = self.execute("free -t | grep ^Mem")
                splitted = output[0].split(":")[1].strip().split()
                buffers = float(splitted[4])
                cache = float(splitted[5])
                free = free_without_buffers_cache + buffers + cache
                free = free / 1e6

            else:

                output = self.execute("free -t | grep buff/cache")
                if len(output) == 0: raise RuntimeError("Unexpected output")

                output = self.execute("free -t | grep ^Mem")
                buffers_and_cache = float(output[0].split(":")[1].strip().split()[4])
                free = free_without_buffers_cache + buffers_and_cache
                free = free / 1e6

        # Return the free amount of virtual memory in gigabytes
        return free * u("GB")

    # -----------------------------------------------------------------

    @property
    def total_space(self):

        """
        This function ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        total = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            if "G" in splitted[1]: total += float(splitted[1].split("G")[0])
            elif "T" in splitted[1]: total += float(splitted[1].split("T")[0]) * 1e3
            elif "M" in splitted[1]: total += float(splitted[1].split("M")[0]) * 1e-3
            elif splitted[1].strip() == "0": pass
            else: raise ValueError("Could not parse line: '" + entry)

        # Return the amount of total memory
        total = total * u("GB")
        if total.value > 999: total = total.to("TB")
        elif total.value < 1: total = total.to("MB")
        return total

    # -----------------------------------------------------------------

    @property
    def used_space(self):

        """
        This fnction ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        used = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            if "G" in splitted[2]: used += float(splitted[2].split("G")[0])
            elif "T" in splitted[2]: used += float(splitted[2].split("T")[0]) * 1e3
            elif "M" in splitted[2]: used += float(splitted[2].split("M")[0]) * 1e-3
            elif splitted[2].strip() == "0": pass
            else: raise ValueError("Could not parse line: '" + entry)

        # Return the amount of used memory
        used = used * u("GB")
        if used.value > 999: used = used.to("TB")
        elif used.value < 1: used = used.to("MB")
        return used

    # -----------------------------------------------------------------

    @property
    def free_space(self):

        """
        This function ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        free = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            if "G" in splitted[3]: free += float(splitted[3].split("G")[0])
            elif "T" in splitted[3]: free += float(splitted[3].split("T")[0]) * 1e3
            elif "M" in splitted[3]: free += float(splitted[3].split("M")[0]) * 1e-3
            elif splitted[3].strip() == "0": pass
            else: raise ValueError("Could not parse line: '" + entry)

        # Return the amount of free memory
        free = free * u("GB")
        if free.value > 999: free = free.to("TB")
        elif free.value < 1: free = free.to("MB")
        return free

    # -----------------------------------------------------------------

    @property
    def free_space_home_directory(self):

        """
        This function returns the usable amount of disk space in the home directory
        :return:
        """

        # Quota command
        if self.host.quota_command is not None:

            # Get the used and max amount of memory for the home directory
            used, total = self.quota[self.home_directory]

            # Calculate the free amount of memory
            free = (total - used).to("GB")

            # Return it in GB
            return free

        # No quota: return free memory on the entire disk
        else: return self.free_space

    # -----------------------------------------------------------------

    @property
    def scheduler(self):

        """
        This property ...
        :return:
        """

        return self.host.scheduler

    # -----------------------------------------------------------------

    @property
    def use_hyperthreading(self):

        """
        This function ...
        :return:
        """

        return self.host.use_hyperthreading

    # -----------------------------------------------------------------

    @property
    def multi_node_communication(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, check whether multi node communication is possible based on
        # the configuration of the current cluster
        if self.scheduler: return self.host.cluster.multi_node_communication

        # If no scheduler is used, raise an error (this function should not get called)
        else: raise RuntimeError("This function should only be called when using a remote with a scheduling system")

    # -----------------------------------------------------------------

    @property
    def tty(self):

        """
        This function ...
        :return:
        """

        output = self.execute("tty")
        return int(output.split("pts/")[1])

    # -----------------------------------------------------------------

    @property
    def ttys(self):

        """
        This function ...
        :return:
        """

        output = self.execute("who | awk '{print $2,$NF}' |grep -v '(:[0-9]'")

        session_numbers = []
        for line in output: session_numbers.append(int(line.split("pts/")[1].split(" ")[0]))
        return session_numbers

    # -----------------------------------------------------------------

    @property
    def cpu_model(self):

        """
        This function ...
        :return:
        """

        #output = self.execute("cat /proc/cpuinfo | grep 'model name'")
        #first_line = output[0]
        #return first_line.split(": ")[1]
        #return first_line

        model_names = []

        output = self.read_lines("/proc/cpuinfo")

        for line in output:

            if line.startswith("model name"): model_names.append(line.split(":")[1])

        return model_names[0]

    # -----------------------------------------------------------------

    @property
    def virtual_memory_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the amount of virtual memory memory per node is defined
        # in the configuration (this is in GB)
        if self.scheduler: return self.host.cluster.memory

        # If no scheduler is used, assume the number of nodes is 1 and get the total virtual memory (total swap)
        else:

            output = self.execute("free -t | grep Swap")
            splitted = output[0].split(":")[1].split()

            # Calculate the free amount of memory in gigabytes
            total_swap = float(splitted[0]) * 1e-6 * u("GB")

            # Return the free amount of virtual memory in gigabytes
            return total_swap

    # -----------------------------------------------------------------

    @property
    def nodes(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of nodes is defined in the host configuration
        if self.scheduler: return self.host.cluster.nodes

        # If no scheduling system is used, assume the system is only concised of one node
        else: return 1

    # -----------------------------------------------------------------

    @property
    def cores_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of cores on the computing nodes is defined in the configuration
        if self.scheduler: return self.host.cluster.cores_per_socket * self.host.cluster.sockets_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            lines = self.execute("lscpu")
            for line in lines:
                if "cpu(s)" in line.lower():
                    value = int(line.split(": ")[1]) / self.threads_per_core
                    from ..tools import numbers
                    if not numbers.is_integer(value): raise ValueError("Something is wrong")
                    return int(value)
            raise RuntimeError("Could not determine the number of cores per node")

    # -----------------------------------------------------------------

    @property
    def threads_per_core(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of threads per core is defined in the configuration
        if self.scheduler: return self.host.cluster.threads_per_core

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            lines = self.execute("lscpu")
            for line in lines:
                if "thread(s) per core" in line.lower(): return int(line.split(": ")[1])
            raise RuntimeError("Could not determine the number of hardware threads per core")

    # -----------------------------------------------------------------

    @property
    def sockets_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of sockets per node is defined in the configuration
        if self.scheduler: return self.host.cluster.sockets_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            lines = self.execute("lscpu")
            for line in lines:
                if "socket(s)" in line.lower(): return int(line.split(": ")[1])
            raise RuntimeError("Could not determine the number of sockets per node")

    # -----------------------------------------------------------------

    @property
    def cores_per_socket(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of cores per socket is defined in the configuration
        if self.scheduler: return self.host.cluster.cores_per_socket

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            lines = self.execute("lscpu")
            for line in lines:
                if "core(s) per socket" in line.lower(): return int(line.split(": ")[1])
            raise RuntimeError("Could not determine the number of cores per socket")

    # -----------------------------------------------------------------

    @property
    def numa_domains(self):

        """
        This function returns the number of NUMA domains (per node)
        :return:
        """

        # If the remote host uses a scheduling system, the number of numa domains is defined in the configuration
        if self.scheduler: return self.host.cluster.numa_domains_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            lines = self.execute("lscpu")
            for line in lines:
                if "numa node(s)" in line.lower(): return int(line.split(": ")[1])
            raise RuntimeError("Could not determine the number of NUMA domains")

    # -----------------------------------------------------------------

    @property
    def numa_cpus(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the CPU ranks
        cpus = [[] for _ in range(self.numa_domains)]

        # Only works for up to 10 domains!
        output = self.execute("lscpu")

        for line in output:

            # Skip irrelevant lines
            if not line.startswith("NUMA node") or line.startswith("NUMA node(s)"): continue

            # Get the NUMA domain
            domain = int(line.split("NUMA node")[1].split(" CPU")[0])

            # Get the list of CPU's
            string = line.split(": ")[1].strip()
            cpu_list = parsing.integer_list(string)

            # Set the CPU list for the current NUMA domain
            cpus[domain] = cpu_list

        # Return the list of CPU lists for each NUMA domain
        return cpus

    # -----------------------------------------------------------------

    def cpus_for_numa_domain(self, domain_index):

        """
        This function ...
        :param domain_index:
        :return:
        """

        return self.numa_cpus[domain_index]

    # -----------------------------------------------------------------

    @property
    def is_multinode(self):

        """
        This function ...
        :return:
        """

        return self.is_executable("pbsnodes")

    # -----------------------------------------------------------------

    def get_node_status(self):

        """
        This function ...
        :return:
        """

        # Execute the pbsnodes command
        output = self.execute("pbsnodes -x")
        string = output[0]

        # Parse the tree
        tree = etree.fromstring(string, parser=etree.XMLParser(remove_blank_text=True))

        # Return the tree
        return tree

    # -----------------------------------------------------------------

    @property
    def node_names(self):

        """
        This function ...
        :return:
        """

        # Get node status
        status = self.get_node_status()

        names = []
        for element in status.getchildren():
            if element.tag != "Node": continue  # doesn't happen
            name = element.xpath("name")[0].text
            names.append(name)

        return names

    # -----------------------------------------------------------------

    @property
    def active_cluster_name(self):

        """
        This function ...
        :return:
        """

        from ..tools import sequences
        return sequences.get_all_equal_value([name.split(".")[1] for name in self.node_names])

    # -----------------------------------------------------------------

    @property
    def cpu_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'top' command to get the current CPU load
        output = self.execute("top -b -n1 | grep 'Cpu(s)' | awk '{print $2 + $4}'")

        # Convert the output into the fraction of CPU that is used
        load = float(output[0]) / 100.0

        # Return the current CPU load
        return load

    # -----------------------------------------------------------------

    @property
    def memory_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        output = self.execute("free -t | grep 'Total'")
        splitted = output[0].split(":")[1].split()

        # Calculate the total and used amount of memory in gigabytes
        total = float(splitted[0]) / 1e6
        used = float(splitted[1]) / 1e6
        #free = float(splitted[2]) / 1e6

        # Return the fraction of virtual memory that is currently in use
        return used / total

    # -----------------------------------------------------------------

    @property
    def platform(self):

        """
        This fucntion ...
        :return:
        """

        lines = []
        lines.append("import platform")
        lines.append("answer = None")
        lines.append("answer = 'MacOS' if platform.system() == 'Darwin' else answer")
        lines.append("answer = 'Windows' if platform.system() == 'Windows' else answer")
        lines.append("answer = 'Linux' if platform.system() == 'Linux' else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        return output[0]

    # -----------------------------------------------------------------

    @property
    def is_macos(self):

        """
        This function ...
        :return:
        """

        return self.operating_system_short == "Darwin"

    # -----------------------------------------------------------------

    @property
    def is_linux(self):

        """
        This function ...
        :return:
        """

        return self.operating_system_short == "Linux"

    # -----------------------------------------------------------------

    @property
    def operating_system_short(self):

        """
        This function ...
        :return:
        """

        output = self.execute("uname -s")
        return output[0]

    # -----------------------------------------------------------------

    @property
    def operating_system_long(self):

        """
        This function ...
        :return:
        """

        output = self.execute("uname -a")
        return output[0]

    # -----------------------------------------------------------------

    @property
    def architecture(self):

        """
        This function ...
        :return:
        """

        # architecture
        lines = []
        lines.append("import platform")
        lines.append("answer = None")
        lines.append("answer = '64bit' if platform.architecture()[0] == '64bit' else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        return output[0]

    # -----------------------------------------------------------------

    @property
    def python_distributions(self):

        """
        This function ...
        :return:
        """

        # Dictionary
        packages = dict()

        # pip list
        output = self.execute("pip list")

        # Loop over the lines in the output
        for entry in output:

            # Known messages that corrupt the output
            if "You are using pip version" in entry: continue
            if "You should consider upgrading via the" in entry: continue

            # Get name and version
            name, version = entry.split(" (")
            version = version[:-1]

            # Add the package name with its version
            packages[name] = version

        # Return the dictionary of python packages with their version numbers
        return packages

    # -----------------------------------------------------------------

    def file_or_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Launch a bash command to check whether the path exists as either a file or a directory on the remote file system
        output = self.execute("if [ -f " + path + " ]; then echo file; elif [ -d " + path + " ]; then echo directory; else echo None; fi")

        # Return the result
        if output[0] == "None": return None
        else: return output[0]

    # -----------------------------------------------------------------

    def evaluate_boolean_expression(self, expression):

        """
        This function ...
        :param expression:
        :return:
        """

        # The commmand
        command = "if [ " + expression + " ]; then echo True; else echo False; fi"
        #print(command)
        # Launch a bash command to check whether the path exists as a directory on the remote file system
        #output = self.execute(command, show_output=False)

        #if len(output) == 0: self.ssh.prompt()

        #print("BEFORE:", self.ssh.before)
        #print("AFTER:", self.ssh.after)

        # Send the command
        self.ssh.sendline(command)

        matched = self.ssh.prompt()

        #print("MATCHED:", matched)
        #print("BEFORE:", list(self.ssh.before))
        #print("AFTER:", list(self.ssh.after))

        output = self.ssh.before.replace(" \r", "").replace("\x1b[K", "").split("\r\n")
        output = output[1:-1]

        #print("OUTPUT:", output)

        return output[0] == "True"

    # -----------------------------------------------------------------

    def is_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        return self.evaluate_boolean_expression("-d '" + path + "'")

    # -----------------------------------------------------------------

    def is_empty(self, directory, ignore_hidden=True):

        """
        This function ...
        :param directory:
        :param ignore_hidden:
        :return:
        """

        files = self.files_in_path(directory, ignore_hidden=ignore_hidden)
        directories = self.directories_in_path(directory, ignore_hidden=ignore_hidden)

        return len(files) + len(directories) == 0

    # -----------------------------------------------------------------

    def is_directory_alternative(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        #print(self.execute("ls " + fs.directory_of(path)))
        paths = self.directories_in_path(fs.directory_of(path))
        return path in paths

    # -----------------------------------------------------------------

    def is_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        #print("path,", path)

        # Launch a bash command to check whether the path exists as a regular file on the remote file system
        return self.evaluate_boolean_expression("-f '" + path + "'")

    # -----------------------------------------------------------------

    def is_subdirectory(self, path, parent_path):

        """
        This function ...
        :param path:
        :param parent_path:
        :return:
        """

        if not self.is_directory(path): raise ValueError("Not a directory: " + path)

        path = self.absolute_path(path)
        parent_path = self.absolute_path(parent_path)
        return path.startswith(parent_path)

    # -----------------------------------------------------------------

    def get_file_hash(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # MacOS
        if self.is_macos:

            command = "md5 '" + path + "'"
            output = self.execute(command)
            if len(output) == 0: raise RuntimeError("No output")
            elif len(output) > 1: raise RuntimeError("Too much output")
            line = output[0]
            # MD5(path) =
            hash = line.split(") = ")[1]

        # Linux
        elif self.is_linux:

            command = "md5sum '" + path + "' | awk '{ print $1 }'"
            output = self.execute(command)
            if len(output) == 0: raise RuntimeError("No output")
            elif len(output) > 1: raise RuntimeError("Too much output")
            hash = output[0]

        # Not supported
        else: raise NotImplementedError("Only MacOS and Linux are supported")

        # Return the hash code
        return hash

    # -----------------------------------------------------------------

    def exists_and_equal_to_local(self, path, local_path):

        """
        This function ...
        :param path:
        :param local_path:
        :return:
        """

        return self.is_file(path) and self.equal_to_local_file(path, local_path)

    # -----------------------------------------------------------------

    def equal_to_local_file(self, path, local_path):

        """
        This function ...
        :param path:
        :param local_path:
        :return:
        """

        # Get hash of local file and of remote file
        hash = self.get_file_hash(path)
        local_hash = fs.get_file_hash(local_path)

        # Return whether the hashes are equal
        return hash == local_hash

    # -----------------------------------------------------------------

    def resolve_environment_variable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.execute("echo $" + name)[0]

    # -----------------------------------------------------------------

    @property
    def username(self):

        """
        This function ...
        :return:
        """

        return self.resolve_environment_variable("USER")

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        return self.host.id

    # -----------------------------------------------------------------

    @property
    def cluster_name(self):

        """
        This function ...
        :return:
        """

        return self.host.cluster_name

    # -----------------------------------------------------------------

    @property
    def cluster_names(self):

        """
        This function ...
        :return:
        """

        return self.host.cluster_names

    # -----------------------------------------------------------------

    @property
    def other_cluster_names(self):

        """
        This function ...
        :return:
        """

        names = self.cluster_names
        names.remove(self.cluster_name)
        return names

    # -----------------------------------------------------------------

    @property
    def pts_git_branches(self):

        """
        This function ...
        :return:
        """

        args = ["git", "branch"]
        output = self.execute(" ".join(args), cwd=self.pts_package_path)

        branches = []
        for line in output.split("\n"):
            if line: branches.append(line.split("* ")[1])
        return branches

    # -----------------------------------------------------------------

    @property
    def pts_git_remotes(self):

        """
        This function ...
        :return:
        """

        args = ["git", "remote"]
        output = self.execute(" ".join(args), cwd=self.pts_package_path)

        remotes = []
        for line in output.split("\n"):
            if line: remotes.append(line)
        return remotes

    # -----------------------------------------------------------------

    def pts_git_remote_url(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #args = ["git", "remote", "show", name]
        #output = self.execute(" ".join(args), cwd=self.pts_package_path)

        #for line in output.split("\n"):
        #    if "Fetch URL" in line: return line.split(": ")[1]

        #raise RuntimeError("Remote '" + name + "' not found!")

        from ..tools import git
        from ..prep.update import fix_repo_url
        try: return git.get_url_repository(self, self.pts_package_path, name)
        except git.AuthenticationError as e:
            url = fix_repo_url(e.url, self.pts_package_path, remote=self)
            return url

    # -----------------------------------------------------------------

    @property
    def skirt_git_branches(self):

        """
        This function ...
        :return:
        """

        args = ["git", "branch"]
        output = self.execute(" ".join(args), cwd=self.skirt_repo_path)

        branches = []
        for line in output.split("\n"):
            if line: branches.append(line.split("* ")[1])
        return branches

    # -----------------------------------------------------------------

    @property
    def skirt_git_remotes(self):

        """
        This function ...
        :return:
        """

        args = ["git", "remote"]
        output = self.execute(" ".join(args), cwd=self.skirt_repo_path)

        remotes = []
        for line in output.split("\n"):
            if line: remotes.append(line)
        return remotes

    # -----------------------------------------------------------------

    def skirt_git_remote_url(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        from ..tools import git
        from ..prep.update import fix_repo_url
        try: return git.get_url_repository(self, self.skirt_repo_path, name)
        except git.AuthenticationError as e:
            url = fix_repo_url(e.url, self.skirt_repo_path, remote=self)
            return url

    # -----------------------------------------------------------------

    @property
    def skirt_root_path(self):

        """
        This function ...
        :return:
        """

        path = self.absolute_path("~/SKIRT")
        return path

    # -----------------------------------------------------------------

    @property
    def skirt_repo_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.skirt_root_path, "git")

    # -----------------------------------------------------------------

    @property
    def pts_root_path(self):

        """
        This function ...
        :return:
        """

        path = self.absolute_path("~/PTS")
        return path

    # -----------------------------------------------------------------

    @property
    def has_pts_root_path(self):

        """
        This function ...
        :return:
        """

        return self.is_directory(self.pts_root_path)

    # -----------------------------------------------------------------

    @property
    def pts_package_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.pts_root_path, "pts")

    # -----------------------------------------------------------------

    @property
    def pts_do_path(self):

        """
        This fnction ...
        :return:
        """

        return fs.join(self.pts_package_path, "do")

    # -----------------------------------------------------------------

    @property
    def pts_main_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.pts_do_path, "__main__.py")

    # -----------------------------------------------------------------

    @property
    def pts_temp_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "temp")
        if not self.is_directory(path) and self.has_pts_root_path: self.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def has_pts_temp_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "temp")
        return self.is_directory(path)

    # -----------------------------------------------------------------

    @property
    def pts_tests_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "test")
        if not self.is_directory(path) and self.has_pts_root_path: self.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @property
    def has_pts_tests_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "test")
        return self.is_directory(path)

    # -----------------------------------------------------------------

    def clear_temp_and_sessions(self):

        """
        This function ...
        :return: 
        """

        if self.has_pts_temp_path: self.clear_pts_temp()
        self.close_all_screen_sessions()
        self.close_all_tmux_sessions()

    # -----------------------------------------------------------------

    def clear_pts_temp(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.info("Clearing temporary data ...")

        # Clear the temporary directory
        self.clear_directory(self.pts_temp_path)

    # -----------------------------------------------------------------

    def clear_pts_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.info("Clearing test data ...")

        # Clear the tests directory
        self.clear_directory(self.pts_tests_path)

    # -----------------------------------------------------------------

    @property
    def pts_run_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "run")
        if not self.is_directory(path): self.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def local_pts_host_run_dir(self):

        """
        This function ...
        :return:
        """

        path = fs.join(introspection.pts_run_dir, self.host.id)
        if not fs.is_directory(path): fs.create_directory(path, recursive=True)
        return path

    # -----------------------------------------------------------------

    @property
    def scratch_path(self):

        """
        This function ...
        :return:
        """

        if self.host.scratch_path is None: return None
        else: return self.absolute_path(self.host.scratch_path)

    # -----------------------------------------------------------------

    @property
    def quota(self):

        """
        This function ...
        :return:
        """

        from ..units.parsing import parse_quantity

        if self.host.quota_command is None: return None

        #print(self.host.quota_command)
        output = self.execute(self.host.quota_command)
        #print(output)

        quotas = dict()

        for line in output:

            if ": used" not in line: continue

            alias = line.split(":")[0].strip()
            if "[" in alias and "]" in alias: alias = alias.split(" [")[0].strip()
            path = self.resolve_environment_variable(alias)

            used = parse_quantity(line.split("used ")[1].split(" (")[0])
            limit = parse_quantity(line.split("quota ")[1].split(" (")[0])

            # Set information
            quotas[path] = (used, limit)

        # Return the quotas
        return quotas

    # -----------------------------------------------------------------

    def retrieve_tasks(self):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Initialize a list to contain the tasks that have been retrieved
        tasks = []

        # Loop over the different entries of the status list
        for path, task_status in self.get_task_status():

            # If a task has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
            # tasks (again) so that its results can be analysed
            if task_status == "retrieved":

                # Open the task file
                task = Task.from_file(path)

                # Add the task to the list
                tasks.append(task)

            # Finished task
            elif task_status == "finished":

                # Open the task file
                task = Task.from_file(path)

                # CHECK WHETHER OUTPUT HAS TO BE DOWNLOADED
                if task.local_output_path is not None:

                    # Debug info
                    self.debug("Retrieving the complete remote output directory ...")
                    self.debug("Local output directory: " + task.local_output_path)
                    self.debug("Remote output directory: " + task.remote_output_path)

                    # Check whether the output directory exists; if not, create it
                    if not fs.is_directory(task.local_output_path): fs.create_directory(task.local_output_path)

                    # Download the PTS task output
                    self.download(task.remote_output_path, task.local_output_path)

                # Local output path not defined
                else: self.warning("Local output path not defined: remote PTS task output will not be retrieved (look on the remote filesystem for the results in '" + task.remote_output_path + "')")

                # Add the retrieved task to the list
                tasks.append(task)

                # If retrieval was succesful, add this information to the task file
                task.retrieved = True
                task.save()

                # Remove the task from the remote
                task.remove_from_remote(self)

        # Return the list of retrieved tasks
        return tasks

    # -----------------------------------------------------------------

    def get_task_status(self):

        """
        This function ..
        :return:
        """

        # Initialize a list to contain the statuses
        entries = []

        # If the remote host does not use a scheduling system
        if not self.scheduler:

            # Search for task files in the local PTS run/host_id directory
            for path in fs.files_in_path(self.local_pts_host_run_dir, extension="task", sort=int):

                # Open the task file
                task = Task.from_file(path)

                # Check whether the task has already been analysed
                if task.analysed: task_status = "analysed"

                # Check whether the task has already been retrieved
                elif task.retrieved: task_status = "retrieved"

                # Not yet retrieved, could be any other state ...
                else:

                    # Get screen name and output path
                    screen_name = task.screen_name
                    output_path = task.remote_output_path
                    log_output_path = task.remote_log_path

                    # Check whether the remote output path exists
                    if not self.is_directory(output_path): task_status = "invalid: remote output directory has been deleted"

                    # The remote output path exists
                    else:

                        # Check whether the log file is present
                        log_path = None
                        for filename in self.files_in_path(log_output_path):
                            if "log" in filename:
                                log_path = fs.join(log_output_path, filename)
                                break

                        # Check whether the report file exists
                        if log_path is not None:

                            if self.is_active_screen(screen_name): task_status = "running"
                            else:

                                # Get the last two lines of the remote log file
                                output = self.read_last_lines(log_path, 2)
                                last_line = output[1]

                                # Check whether the last line states that the task has finished
                                if "Finished" in last_line: task_status = "finished"
                                elif "Error:" in last_line: task_status = "crashed: " + last_line
                                else: task_status = "aborted"

                        # If the log file does not exist, the task has not started yet or has been cancelled
                        else:

                            # The task has not started or it's screen session has been cancelled
                            if self.is_active_screen(screen_name): task_status = "queued"
                            else: task_status = "cancelled"

                # Add the task properties to the list
                entries.append((path, task_status))

        # If the remote has a scheduling system for launching jobs
        #else: raise NotImplementedError("Running PTS tasks on remote schedulers is not supported yet")

        # Return the list of task properties
        return entries

    # -----------------------------------------------------------------

    def _task_ids_in_use(self):

        """
        This function ...
        :return:
        """

        # Check the contents of the local run directory to see which task id's are currently in use
        current_ids = []
        for name in fs.files_in_path(self.local_pts_host_run_dir, extension="task", returns="name"):

            # Get the task ID and add it to the list
            current_ids.append(int(name))

        # Return the list of currently used ID's
        return current_ids

    # -----------------------------------------------------------------

    def _new_task_id(self):

        """
        This function ...
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._task_ids_in_use()

        # Sort the current task ID's and find the lowest 'missing' integer number
        if len(current_ids) > 0:
            current_ids = sorted(current_ids)
            task_id = max(current_ids) + 1
            for index in range(max(current_ids)):
                if current_ids[index] != index:
                    task_id = index
                    break

            # Return the task ID
            return task_id

        # If no task ID's are currently in use, return 0
        else: return 0

# -----------------------------------------------------------------
