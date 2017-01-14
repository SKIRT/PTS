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

# Import astronomical modules
from astropy.utils import lazyproperty
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .host import Host, HostDownException
from .vpn import VPN
from ..tools.logging import log
from ..tools import parsing
from ..tools import filesystem as fs
from ..tools import time
from ..basics.task import Task
from ..tools import introspection
from ..tools.introspection import possible_cpp_compilers, possible_mpi_compilers, possible_mpirun_names
from .python import RemotePythonSession

# -----------------------------------------------------------------

def is_available(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    remote = Remote()

    try: remote.setup(host_id)
    except HostDownException: return False

    del remote
    return True

# -----------------------------------------------------------------

def active_keys():

    """
    This function ...
    :return:
    """

    output = subprocess.check_output(["ssh-add", "-L"])

    names = []
    for line in output.split("\n"):
        if not line: continue
        _, key, path = line.split(" ")
        name = fs.name(path)
        names.append(name)

    return names

# -----------------------------------------------------------------

def add_key(name):

    """
    This function ...
    :param name:
    :return:
    """

    # Check if the key exists
    key_path = fs.absolute_path(fs.join("~/.ssh", name))
    if not fs.is_file(key_path): raise ValueError("The key '" + name + "' does not exist in the .ssh directory")

    # Add the key
    subprocess.call(["ssh-add", key_path])

    # password is asked:
    # Create the pexpect child instance
    #child = pexpect.spawn("ssh-add " + key_path, timeout=30)
    #password = "tokiotokio"
    #if password is not None:
    #    child.expect([': '])
    #    child.sendline(password)

# -----------------------------------------------------------------

class Remote(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
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

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster_name=None, login_timeout=15):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :param login_timeout:
        :return:
        """

        # Create the host object
        self.host = Host(host_id, cluster_name)

        # If a VPN connection is required for the remote host
        if self.host.requires_vpn: self.connect_to_vpn()

        # Check if key is active
        if self.host.key is not None:
            if self.host.key not in active_keys(): add_key(self.host.key)

        # Make the connection
        self.login(login_timeout)

        # Swap cluster
        if self.host.cluster_name is not None: self.swap_cluster(self.host.cluster_name)

        # LMOD_DISABLE_SAME_NAME_AUTOSWAP
        #self.define_environment_variable("LMOD_DISABLE_SAME_NAME_AUTOSWAP", "yes")

    # -----------------------------------------------------------------

    def __del__(self):

        """
        The destructor ...
        :return:
        """

        # Disconnect from the remote host
        self.logout()

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

    def load_module(self, module_name, show_output=False):

        """
        This function ...
        :param module_name:
        :param show_output:
        :return:
        """

        self.execute("module load " + module_name, show_output=show_output)

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

        # Inform the user
        log.info("Unloading all modules ...")

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

    def get_module_versions(self, module_name, exact=True):

        """
        This function ...
        :param module_name:
        :param exact:
        :return:
        """

        if exact: command = "module spider " + module_name
        else: command = "module spider -r " + module_name

        print(command)
        print(self.ssh.before)
        output = self.execute(command, show_output=True)

        versions = []

        triggered = False
        for line in output:

            if triggered:

                if not line.strip(): break

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

        output = self.execute("module list")

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

    def fix_configuration_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking and fixing shell configuration files ...")

        # Check if bashrc exists, if not create it
        bashrc_path = fs.join(self.home_directory, ".bashrc")
        if not self.is_file(bashrc_path): self.touch(bashrc_path)

        # Check whether .profile exists
        profile_path = fs.join(self.home_directory, ".profile")
        bash_profile_path = fs.join(self.home_directory, ".bash_profile")

        # If bash_profile file exists
        if self.is_file(bash_profile_path):

            # Check if it points to the bashrc file
            links_bashrc = ".bashrc" in ";".join(self.read_lines(bash_profile_path))

            # Add link to bashrc
            if not links_bashrc: self._add_bashrc_link(bash_profile_path)

        # If profile file exists
        elif self.is_file(profile_path):

            # Check if it points to the bashrc file
            links_bashrc = ".bashrc" in ";".join(self.read_lines(profile_path))

            # Add link to bashrc
            if not links_bashrc: self._add_bashrc_link(profile_path)

        # Profile and bash_profile files do not exist
        else:

            # Make bash profile
            self.touch(bash_profile_path)

            # Make it link to bashrc
            self._add_bashrc_link(bash_profile_path)

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

    @property
    def has_mpi_compiler(self):

        """
        This function ...
        :return:
        """

        for name in possible_mpi_compilers:
            if self.is_executable(name): return True

    # -----------------------------------------------------------------

    @property
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

    def swap_cluster(self, cluster_name):

        """
        This function ...
        :param cluster_name:
        :return:
        """

        # Check if not already the loaded cluster
        if "cluster/" + cluster_name in self.loaded_modules: return

        # Swap to requested cluster
        self.execute("module swap cluster/" + cluster_name)

    # -----------------------------------------------------------------

    def find_and_load_python(self):

        """
        This function ...
        :return:
        """

        if self.has_lmod: self.load_python_module()
        return self.python_path

    # -----------------------------------------------------------------

    def load_python_module(self):

        """
        This function ...
        :return:
        """

        # Find Python module versions
        versions = self.get_module_versions("Python")

        latest_version = None
        latest_floating_python_version = None
        latest_intel_version = None
        latest_intel_year = None

        for version in versions:

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

        self.load_module(latest_version)  # Load the python module

    # -----------------------------------------------------------------

    def find_and_load_git(self):

        """
        This function ...
        :return:
        """

        # Load intel compiler toolkit if possible
        if self.has_lmod:

            versions = self.get_module_versions("git")

            latest_version = None
            latest_module_name = None
            for version in versions:

                # Get version
                the_version = version.split("/")[1].split("-")[0]

                if latest_version is None or the_version > latest_version:
                    latest_version = the_version
                    latest_module_name = version

            if latest_module_name is None: raise RuntimeError("Git not available from the modules")

            self.load_module(latest_module_name)

        # Get path and version
        git_path = self.find_executable("git")
        git_version = self.version_of("git")

        # Return
        return git_path, git_version

    # -----------------------------------------------------------------

    def find_and_load_cpp_compiler(self):

        """
        This function ...
        :return:
        """

        # Load intel compiler toolkit if possible
        if self.has_lmod: self.load_intel_compiler_toolkit()

        # Search for the compiler and return its path
        return self.cpp_compiler_path

    # -----------------------------------------------------------------

    def load_intel_compiler_toolkit(self):

        """
        This function ...
        :return:
        """

        # Find latest Intel Compiler Toolkit version and load it
        intel_version = self._find_latest_iimpi_version_module()
        if intel_version is None:
            log.warning("Intel Cluster Toolkit Compiler Edition could not be found")
        elif intel_version not in self.loaded_modules: self.load_module(intel_version) # Load the module

    # -----------------------------------------------------------------

    def _find_latest_iimpi_version_module(self):

        """
        This function ...
        :return:
        """

        # Load the Intel Cluster Toolkit Compiler Edition (inludes Intel MPI)
        versions = self.get_module_versions("iimpi")

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

    def find_and_load_mpi_compiler(self):

        """
        This function ...
        :return:
        """

        if self.has_lmod: self.load_intel_compiler_toolkit()
        return self.mpi_compiler_path

    # -----------------------------------------------------------------

    def find_and_load_qmake(self):

        """
        This function ...
        :return:
        """

        # Use modules or not
        if self.has_lmod: return self._check_qt_remote_lmod()
        else: return self._check_qt_remote_no_lmod()

    # -----------------------------------------------------------------

    def _check_qt_remote_lmod(self):

        """
        This function ...
        :return:
        """

        # Find latest Qt5 version
        qt5_version = self._find_latest_qt_version_module()

        # If no version could be found, give an error
        if qt5_version is None: raise RuntimeError("No Intel-compiled version of Qt 5 could be found between the available modules")

        # Load the module
        self.load_module(qt5_version, show_output=True)

        # Get the qmake path
        qmake_path = self.find_executable("qmake")

        # Get the Intel compiler version
        #if "intel" in qt5_version: intel_version = qt5_version.split("intel-")[1].split("-")[0]
        #else: intel_version = qt5_version.split("ictce-")[1].split("-")[0]

        # Return the qmake path
        return qmake_path

    # -----------------------------------------------------------------

    def _check_qt_remote_no_lmod(self):

        """
        This function ...
        :return:
        """

        # Keep the qmake paths in a list to decide later which one we can use
        qmake_paths = []

        # Search for qmake in the home directory
        command = "find " + self.home_directory + "/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.execute(command)

        # Search for qmake in the /usr/local directory
        command = "find /usr/local/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.execute(command)

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
        qmake_path = self.find_executable("qmake")
        if qmake_path is not None: qmake_paths.append(qmake_path)

        latest_version = None

        latest_qmake_path = None

        # Loop over the qmake paths
        for qmake_path in qmake_paths:

            # Get the version
            output = self.execute(qmake_path + " -v")

            qt_version = output[1].split("Qt version ")[1].split(" in")[0]

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

    def _find_latest_qt_version_module(self):

        """
        This function ...
        :return:
        """

        versions = []

        # Check 'Qt5' module versions
        versions += self.get_module_versions("Qt5")

        # Check 'Qt' module versions
        versions += self.get_module_versions("Qt")

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
        log.info("Connecting to vpn service '" + self.host.vpn.service + "' ...")

        # Connect to the VPN service
        self.vpn = VPN(self.host.vpn.service)
        self.vpn.connect(self.host.vpn.user, self.host.vpn.password, self.host.vpn.secret, self.host.vpn.prompt_time_delay)

    # -----------------------------------------------------------------

    def login(self, login_timeout=30):

        """
        This function ...
        :param login_timeout:
        :return:
        """

        # Inform the user
        log.info("Logging in to the remote environment on host '" + self.host.id + "' ...")

        # Try connecting to the remote host
        try: self.connected = self.ssh.login(self.host.name, self.host.user, self.host.password, port=self.host.port, login_timeout=login_timeout)
        except ExceptionPexpect: raise HostDownException()

        # Check whether connection was succesful
        if not self.connected: raise RuntimeError("Connection failed")

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        if log is not None: log.info("Logging out from the remote environment on host '" + self.host_id + "' ...")
        # the conditional statement is because of this error message during destruction at the end of a script:
        # Exception AttributeError: "'NoneType' object has no attribute 'info'" in <bound method SkirtRemote.__del__ of <pts.core.simulation.remote.SkirtRemote object at 0x118628d50>> ignored

        # Disconnect
        if self.connected:

            # TODO: Check if there are python sessions open?

            self.ssh.logout()
            self.connected = False

        # Disconnect from the VPN service if necessary
        #if self.vpn is not None: self.vpn.disconnect()

    # -----------------------------------------------------------------

    def start_screen(self, name, local_script_path, script_destination, screen_output_path=None,
                     keep_remote_script=False, attached=False):

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

    def kill_job(self, id):

        """
        This function ...
        :param id:
        :return:
        """

        # Stop the job with the specified ID
        self.execute("qdel " + str(id), output=False)

    # -----------------------------------------------------------------

    def screen_names(self):

        """
        This function ...
        :return:
        """

        output = self.execute("screen -ls | grep \(")

        # Get the names
        names = []
        for line in output:
            name = line.split(".")[1].split("(")[0].strip()
            names.append(name)

        # Return the screen names
        return names

    # -----------------------------------------------------------------

    def screen_numbers(self):

        """
        This function ...
        :return:
        """

        output = self.execute("screen -ls | grep \(")

        # Get the numbers
        numbers = []
        for line in output:
            number = int(line.split(".")[0])
            numbers.append(number)

        # Return the screen numbers
        return numbers

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

    def download_from_url_to(self, url, filepath, show_output=True):

        """
        This function ...
        :return:
        """

        command = "wget " + url + " -O " + filepath
        self.execute(command, show_output=show_output)

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
        log.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, unique_session_name + ".cfg")

        # Save the configuration file to the temporary directory
        config.saveto(temp_conf_path)

        # Debugging
        log.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

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

            for name in input_dict:

                filename = name + "." + input_dict[name].default_extension

                local_filename = fs.join(temp_path, filename) # local temporary path

                # TODO: save the input item here

                filepath = fs.join(input_path, filename) # remote temporary path

                # TODO: UPLOAD THE INPUT item here

                # TODO: delete the local input item

        ####

        # Debugging
        log.debug("Creating a script for remote execution ...")

        # Determine the path to the remote equivalent of this file
        remote_main_path = fs.join(self.pts_package_path, "do", "__main__.py")

        # Create a bash script
        temp_script_path = fs.join(temp_path, unique_session_name + ".sh")

        # Write the lines to the script file
        with open(temp_script_path, 'w') as script_file:

            #script_file.write("#!/usr/bin/env python\n")
            #script_file.write("# -*- coding: utf8 -*-\n")
            #script_file.write("\n")
            script_file.write("python " + remote_main_path + " --configfile " + remote_conf_path + " " + command + "\n")

        # Execute the script
        # name, local_script_path, script_destination, screen_output_path=None, keep_remote_script=False
        self.start_screen(unique_session_name, temp_script_path, remote_temp_path, screen_output_path=remote_temp_path, keep_remote_script=keep_remote_output)

        # Remove the local script file
        fs.remove_file(temp_script_path)

        # Remove the remote temporary directory
        if keep_remote_output: log.info("Remote output will be placed in '" + remote_output_path + "'")

        # Create a new Task object
        task = Task(command, config.to_string())

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

        # Change the working directory if necessary
        if cwd is not None: original_cwd = self.change_cwd(cwd)
        else: original_cwd = None

        # Send the command
        self.ssh.sendline(command)

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Retrieve the output if requested
        if expect is None: matched = self.ssh.prompt(timeout=timeout)
        else: matched = self.ssh.expect(expect, timeout=timeout)

        # If an extra EOF is used before the actual output line (don't ask me why but I encounter this on the HPC UGent infrastructure), do prompt() again
        if contains_extra_eof: matched = self.ssh.prompt()

        # If the command could not be sent, raise an error
        if not matched and expect_eof and not contains_extra_eof: raise RuntimeError("The command could not be sent")

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Ignore the first and the last line (the first is the command itself, the last is always empty)
        # Trial and error to get it right for HPC UGent login nodes; don't know what is happening
        if contains_extra_eof:
            splitted = self.ssh.before.replace('\x1b[K', '').split("\r\n")
            #print(splitted)
            if splitted[-1] == "": the_output = splitted[output_start:-1]
            else: the_output = splitted[output_start:]
        else:
            splitted = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")
            #print(splitted)
            if splitted[-1] == "": the_output = splitted[output_start:-1]
            else: the_output = splitted[output_start:]

        # Set the working directory back to the original
        if original_cwd is not None: self.change_cwd(original_cwd)

        # Return the output
        if output:

            # Something weird, encountered on Cindy,
            # before is "staggering behind", we are already giving commands further and a 'prompt' has been missed or something like that...
            """
            if len(the_output) == 0:

                new_output = []

                new_output += self.ssh.before.split("\r\n")
                new_output += self.ssh.after.split("\r\n")

                # print(self.ssh.before)
                self.ssh.sendline("\n")
                self.ssh.prompt()
                self.ssh.prompt()
                # print(self.ssh.before)
                # output = self.execute("echo $HOME")
                # print("OUTPUT", output[0])
                #print(self.ssh.before)

                #splitted = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")
                #return output[0]

                new_output += self.ssh.before.split("\r\n")
                new_output += self.ssh.after.split("\r\n")

                print(new_output)

                if splitted[-1] == "": the_output = splitted[output_start:-1]
                else: the_output = splitted[output_start:]
            """

            return the_output

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

        # Loop over the lines
        for line in args:

            # If string
            if isinstance(line, basestring):

                # Send the command
                self.ssh.sendline(line)

            # If tuple
            elif isinstance(line, tuple):

                # Expect
                if len(line) == 3 and line[2]:
                    #index = self.ssh.expect([self.ssh.PROMPT, line[0]]) # this is not working, why?
                    index = self.ssh.expect(["$", line[0]])
                    if index == 0: pass
                    elif index == 1: self.ssh.sendline(line[1])
                    #eof = self.ssh.prompt()
                else:
                    self.ssh.expect(line[0])
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

    def start_python_session(self, assume_pts=True):

        """
        This function ...
        :param assume_pts: assume PTS is present, import some basic PTS tools
        :return:
        """

        # Start python session and return it
        return RemotePythonSession(self, assume_pts=assume_pts)

    # -----------------------------------------------------------------

    def execute_python_interactive(self, lines, show_output=False):

        """
        This function ...
        :param lines:
        :param show_output:
        :return:
        """

        # Create python session
        python = self.start_python_session()

        # Execute the lines
        output = python.send_lines(lines, show_output=show_output)

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

        # Use the 'mv' command to rename the file
        self.execute("mv " + old_path + " " + new_path)

    # -----------------------------------------------------------------

    def remove_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.execute("rm -rf " + path, output=False)

    # -----------------------------------------------------------------

    def remove_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.execute("rm " + path, output=False)

    # -----------------------------------------------------------------

    def change_cwd(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        original_cwd = self.working_directory

        # Try to change the directory, give an error if this fails
        #output = self.execute("cd " + path)
        #if len(output) > 0 and "No such file or directory" in output[0]: raise RuntimeError("The directory does not exist")

        #return original_cwd

        # Send the command
        self.ssh.sendline("cd " + path)
        self.ssh.prompt()
        output = self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]
        if len(output) > 0 and "No such file or directory" in output[0]: raise RuntimeError("The directory does not exist")

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

            command = "ls -R /path | awk '/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\n'
            command += "NF&&f{ print "
            command += 's"/"$0 }'
            command += "'"

            output = self.execute(command)

            paths = output

            # Return the paths
            return paths

        else:

            directories = self.execute("for i in $(ls */); do echo ${i%%/}; done")
            if "cannot access */" in directories[0]: directories = []

            files = self.execute("for f in *; do [[ -d $f ]] || echo $f; done")

            # Get paths
            items = directories + files
            paths = [fs.join(path, name) for name in items]

            # Return the paths
            return paths

    # -----------------------------------------------------------------

    def directories_in_path(self, path, startswith=None, recursive=False):

        """
        This function ...
        :param path:
        :param startswith:
        :param recursive:
        :return:
        """

        # List the directories in the provided path
        if recursive:
            command = "ls -R -d */ | awk '\n/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\nNF&&f{ print '
            command += 's"/"$0 }'
            command += "'"
            output = self.execute(command, cwd=path)
            paths = [dirpath for dirpath in output if self.is_directory(dirpath)]
        else:
            #print(path)
            output = self.execute("for i in $(ls -d */); do echo ${i%%/}; done", cwd=path)
            #print(output)
            if "cannot access */" in output[0]: return []
            paths = [fs.join(path, name) for name in output]

        # Filter
        filtered_paths = []
        for path in paths:
            if startswith is not None:
                name = fs.name(path)
                if not name.startswith(startswith): continue
            filtered_paths.append(path)

        # Return the list of directory paths
        return filtered_paths

    # -----------------------------------------------------------------

    def files_in_path(self, path, recursive=False):

        """
        This function ...
        :param path:
        :param recursive:
        :return:
        """

        # List the files in the provided path
        if recursive:
            command = "ls -R " + path + " | awk '\n/:$/&&f{s=$0;f=0}\n/:$/&&!f{sub(/:$/,"
            command += '"");s=$0;f=1;next}\nNF&&f{ print '
            command += 's"/"$0 }'
            command += "'"
            output = self.execute(command, cwd=path)
            paths = [filepath for filepath in output if self.is_file(filepath)]
        else:
            output = self.execute("for f in *; do [[ -d $f ]] || echo $f; done", cwd=path)
            paths = [fs.join(path, name) for name in output]

        # Return the list of files
        return paths

    # -----------------------------------------------------------------

    def to_home_directory(self):

        """
        This function ...
        """

        # Navigate to the home directory
        self.execute("cd ~", output=False)

    # -----------------------------------------------------------------

    def create_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create the remote directory
        self.execute("mkdir " + path, output=False)

    # -----------------------------------------------------------------

    def create_directories(self, *paths):

        """
        This function ...
        :return:
        """

        # Create the remote directories
        self.execute("mkdir " + " ".join(paths), output=False)

    # -----------------------------------------------------------------

    def create_directory_in(self, base_path, name):

        """
        This function ...
        :param base_path:
        :param name:
        :return:
        """

        directory_path = fs.join(base_path, name)
        self.create_directory(directory_path)
        return directory_path

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

        command = "echo '" + line + "' > " + filepath
        self.execute(command, output=False)

    # -----------------------------------------------------------------

    def append_line(self, filepath, line):

        """
        This function ...
        :param filepath:
        :param line:
        :return:
        """

        command = "echo '" + line + "' >>" + filepath
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

    def read_lines(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Expand the path to absolute form
        path = self.absolute_path(path)

        # Load the text file into a variable
        self.execute("value='cat " + path + "'")

        # Print the variable to the console, and obtain the output
        for line in self.execute('echo "$($value)"'): yield line

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

    def download(self, origin, destination, timeout=60, new_name=None, compress=False, show_output=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug(): show_output = True

        # Construct the command string
        copy_command = "scp "
        if compress: copy_command += "-C "

        # Add the host address
        copy_command += self.host.user + "@" + self.host.name + ":"

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            #if self.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if self.is_file(origin): copy_command += "'" + origin + "' "

            # Check if it represents a directory
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\ ") + "/* " + "-r "
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\\ ") + "/* "
            elif self.is_directory(origin): copy_command += "'" + origin + "/*' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file on the remote host")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist remotely
            for file_path in origin:
                if not self.is_file(file_path): raise ValueError("The file " + file_path + " does not exist on the remote host")

            # Escape possible space characters
            #origin = [path.replace(" ", "\\\ ") for path in origin]
            origin = ["'" + path + "'" for path in origin]

            # Add a quotation mark character because the seperate file paths are going to be separated by spaces
            # (the command is going to be of the form scp username@ip.of.server.copyfrom:"file1.log file2.log" "~/yourpathtocopy")
            copy_command += '"'

            # Add the file paths to the command string
            copy_command += " ".join(origin)

            # Add another quotation mark to identify the end of the filepath list
            copy_command += '" '

        # Add the destination path to the command
        #copy_command += destination.replace(" ", "\\\ ") + "/"
        copy_command += "'" + destination + "/'"
        if new_name is not None: copy_command += new_name

        # Debugging
        log.debug("Copy command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        if not show_output:

            # Temporary file for output of the scp command
            #temp_file_path = tempfile.mktemp()
            #temp_file = open(temp_file_path, 'w')

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        else: child.logfile = sys.stdout

        # Execute the command and get the output
        child.expect(pexpect.EOF, timeout=None)
        child.close()

        if not show_output:

            # Retrieve file contents -- this will be
            # 'First line.\nSecond line.\n'
            stdout = child.logfile.getvalue()

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(stdout)

            # Get the output lines
            lines = stdout.split("\n")

            # Check for messages that signal an error that has occured
            for line in lines:
                if "not a regular file" in line: raise ValueError(line)

            # Debugging: show the output of the scp command
            log.debug("Copy stdout: " + str(" ".join(lines)))

    # -----------------------------------------------------------------

    def upload(self, origin, destination, timeout=60, new_name=None, compress=False, show_output=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug(): show_output = True

        # Construct the command string
        copy_command = "scp "
        if compress: copy_command += "-C "

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            #if fs.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if fs.is_file(origin): copy_command += "'" + origin + "' "

            # Check if it represents a directory
            #elif fs.is_directory(origin): copy_command += "-r " + origin.replace(" ", "\\\ ") + "/ "
            elif fs.is_directory(origin): copy_command += "-r '" + origin + "/' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist locally
            for file_path in origin:
                if not fs.is_file(file_path): raise ValueError("The file " + file_path + " does not exist")

            origin = ["'" + path + "'" for path in origin]

            # Add the file paths to the command string
            copy_command += " ".join(origin) + " "

        # Invalid argument
        else: raise ValueError("The origin must be a string or a list of strings")

        # Add the host address and the destination directory
        #copy_command += self.host.user + "@" + self.host.name + ":" + destination.replace(" ", "\\\ ") + "/"
        copy_command += self.host.user + "@" + self.host.name + ":'" + destination + "/'"
        if new_name is not None: copy_command += new_name

        # Debugging
        log.debug("Copy command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        if not show_output:

            # Temporary file for output of the scp command
            #temp_file_path = tempfile.mktemp()
            #temp_file = open(temp_file_path, 'w')

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        else: child.logfile = sys.stdout

        # Execute the command and get the output
        try:
            child.expect(pexpect.EOF, timeout=None)
        except pexpect.EOF:
            pass
        child.close()

        if not show_output:

            # Retrieve file contents -- this will be
            # 'First line.\nSecond line.\n'
            stdout = child.logfile.getvalue()

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(stdout)

            # Get the output lines
            lines = stdout.split("\n")

            # Check for messages that signal an error that has occured
            for line in lines:
                if "not a regular file" in line: raise ValueError(line)

            # Debugging: show the output of the scp command
            log.debug("Copy stdout: " + str(" ".join(lines)))

    # -----------------------------------------------------------------

    def version_of(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Execute
        output = self.execute(name + " --version")

        # Return the relevant portion of the output
        return output[0]

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

        return self.is_executable("pts")

    # -----------------------------------------------------------------

    @property
    def pts_version(self):

        """
        This function ...
        :return:
        """

        if self.has_pts: return self.execute("pts --version")[0]
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

    def find_executable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #print(name)

        # Get the output of the 'which' command
        output = self.execute("which " + name)

        #print(output)

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
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        #output = self.execute("free -t | grep 'Total'")
        output = self.execute("free -t | grep buffers/cache")
        splitted = output[0].split(":")[1].split()

        # Calculate the free amount of memory in gigabytes
        free = float(splitted[1]) / 1e6

        # Return the free amount of virtual memory in gigabytes
        return free

    # -----------------------------------------------------------------

    @property
    def free_space(self):

        """
        This function ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        total = 0.0
        used = 0.0
        free = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            total += float(splitted[1].split("G")[0]) if "G" in splitted[1] else float(splitted[1].split("T")[0]) * 1e3
            used += float(splitted[2].split("G")[0]) if "G" in splitted[2] else float(splitted[2].split("T")[0]) * 1e3
            free += float(splitted[3].split("G")[0]) if "G" in splitted[3] else float(splitted[3].split("T")[0]) * 1e3

        # Return the amount of free memory in gigabytes
        return free

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
            total_swap = float(splitted[0]) * 1e-6 * Unit("GB")

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

            # Use the 'lscpu' command to obtain the total number of CPU's (=hardware threads!)
            output = self.execute("lscpu | grep '^CPU(s)'")
            cpus = int(float(output[0].split(":")[1]))

            # Return the number of physical cores
            return cpus / self.threads_per_core

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

            # Use the 'lscpu' command to get the number of hardware threads per core
            output = self.execute("lscpu | grep '^Thread(s) per core'")
            threads_per_core = int(float(output[0].split(":")[1]))

            # Return the amount of hyperthreads or 'hardware' threads per physical core
            return threads_per_core

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

            # Use the 'lscpu' command to get the number of NUMA domains
            output = self.execute("lscpu | grep '^Socket(s):'")
            nsockets = int(output[0])

            # Return the number of sockets
            return nsockets

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

            # Use the 'lscpu' command
            output = self.execute("lscpu | grep '^Core(s) per socket:'")
            ncores = int(output[0])

            # Return the number of cores per socket
            return ncores

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

            # Use the 'lscpu' command to get the number of NUMA domains
            output = self.execute("lscpu | grep '^NUMA node(s):'")
            numa_domains = int(output[0])

            # Return the number of NUMA domains
            return numa_domains

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
        free = float(splitted[2]) / 1e6

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
    def operating_system(self):

        """
        This function ...
        :return:
        """

        output = self.execute("uname")
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

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        output = self.execute(command)
        return output[0] == "True"

    # -----------------------------------------------------------------

    def is_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        return self.evaluate_boolean_expression("-d " + path)

    # -----------------------------------------------------------------

    def is_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Launch a bash command to check whether the path exists as a regular file on the remote file system
        return self.evaluate_boolean_expression("-f " + path)

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
    def pts_package_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.pts_root_path, "pts")

    # -----------------------------------------------------------------

    @property
    def pts_temp_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "temp")
        if not self.is_directory(path): self.create_directory(path)
        return path

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
                    log.debug("Retrieving the complete remote output directory ...")
                    log.debug("Local output directory: " + task.local_output_path)
                    log.debug("Remote output directory: " + task.remote_output_path)

                    # Check whether the output directory exists; if not, create it
                    if not fs.is_directory(task.local_output_path): fs.create_directory(task.local_output_path)

                    # Download the PTS task output
                    self.download(task.remote_output_path, task.local_output_path)

                # Local output path not defined
                else: log.warning("Local output path not defined: remote PTS task output will not be retrieved (look on the remote filesystem for the results in '" + task.remote_output_path + "')")

                # Add the retrieved task to the list
                tasks.append(task)

                # If retrieval was succesful, add this information to the task file
                task.retrieved = True
                task.save()

                ## REMOVE REMOTE OUTPUT IF REQUESTED
                if task.remove_remote_output:

                    # Remove the temporary PTS directory if it contains the output directory
                    if self.is_subdirectory(task.remote_output_path, task.remote_temp_pts_path): self.remove_directory(task.remote_temp_pts_path)
                    else:
                        # Remove the output directory and the temporary directory seperately
                        self.remove_directory(task.remote_output_path)
                        self.remove_directory(task.remote_temp_pts_path)

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
                                output = self.execute("tail -2 " + log_path)
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
