#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.installation_commands List the commands necessary to install SKIRT/PTS on a remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.remote.modules import Modules
from pts.core.basics.log import log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import introspection
from pts.core.tools import git
from pts.core.prep.installation import qt_url, qt_configure_options
from pts.core.basics.map import Map

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host ID", choices=find_host_ids())
definition.add_flag("private", "install SKIRT from private (development) repository")

# Create the configuration
config = parse_arguments("installation_commands", definition)

# -----------------------------------------------------------------

installation_commands = defaultdict(list)

# -----------------------------------------------------------------

# Create the remote
remote = Remote()
if not remote.setup(config.remote):
    log.error("The remote host is not available")
    exit()

# -----------------------------------------------------------------

# Check the modules
modules = Modules(remote)

# -----------------------------------------------------------------

# Qt is present
if modules.paths["qt"] is not None:

    # Get the version
    version = modules.versions["qt"]

    # Check if module has to be loaded
    if modules.names["qt"] is not None:
        installation_commands["qt"].append(Map(line="module load " + modules.names["qt"], comment="Load qt module (" + version + ")"))

    installation_commands["qt"].append(Map(comment="Qmake is present (" + version + "), no installation required"))

else:

    installation_commands["qt"].append(Map(comment="Installation of Qt"))

    #temp_path = remote.home_directory
    #path = fs.join(temp_path, "qt.tar.gz")

    filepath = "~/qt.tar.gz"
    download_command = "wget " + qt_url + " -O " + filepath

    installation_commands["qt"].append(Map(line=download_command, comment="Download the installer"))

    #decompress_path = self.remote.create_directory_in(temp_path, "Qt-install")

    decompress_path = "~/Qt-install"
    decompress_command = "tar -zxvf " + filepath + " -C " + decompress_path

    installation_commands["qt"].append(Map(line=decompress_command, comment="Decompress the installer"))

    installation_commands["qt"].append(Map(line="cd " + decompress_path + "/qt-everywhere-opensource-..", comment="Navigate to the installation directory"))

    configure_command = "./configure " + " ".join(qt_configure_options)
    make_command = "make"
    install_command = "make install"

    installation_commands["qt"].append(Map(lines=[configure_command, make_command, install_command], comment="build Qt"))

# -----------------------------------------------------------------

installation_commands["directories"].append(Map(comment="Create the necessary directories"))
installation_commands["directories"].append(Map(line="cd"))
installation_commands["directories"].append(Map(line="mkdir SKIRT", comment="Create the SKIRT root directory"))
installation_commands["directories"].append(Map(lines=["cd SKIRT", "mkdir git run release debug"], comment="Create subdirectories"))

# -----------------------------------------------------------------

if config.private:
    url = introspection.private_skirt_https_link
    needs_login = True
else:
    url = introspection.public_skirt_https_link
    needs_login = False

# Decompose
host, user_or_organization, repo_name, _, __ = git.decompose_repo_url(url)

# Find the account file for the repository host (e.g. github.ugent.be)
if introspection.has_account(host): username, password = introspection.get_account(host)
elif needs_login:
    username = "[your_username]"
    password = "[your_password]"
else: username = password = None

# Compose HTTPS link
url = git.compose_https(host, user_or_organization, repo_name, username, password)

# -----------------------------------------------------------------

if modules.names["git"] is not None:
    version = modules.versions["clone"]
    installation_commands["clone"].append(Map(line="module load " + modules.names["git"], comment="Load git module (" + version + ")"))

# -----------------------------------------------------------------

skirt_repo_path = "~/SKIRT/git"

# Set the clone command
command = "git clone " + url + " " + skirt_repo_path

# Add commands
installation_commands["clone"].append(Map(command, comment="Clone the repository"))

# Module purge
if modules.names["git"] is not None:
    installation_commands["clone"].append(Map(line="module purge", comment="Unload all modules (to avoid conflicts)"))

# Clone
#self.remote.execute(command, show_output=log.is_debug)

# Get the git version
#git_version = git.get_short_git_version(self.skirt_repo_path, self.remote)

# Show the git version
#log.info("The git version to be installed is '" + self.git_version + "'")

# Determine SKIRT and FitSKIRT main paths
#skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
#fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")

# Add
#comment = "Added by the Python Toolkit for SKIRT (PTS)"
#self.remote.add_to_environment_variable("PATH", skirt_main_path, comment=comment, in_shell=True)
#self.remote.add_to_environment_variable("PATH", fitskirt_main_path, comment=comment, in_shell=True)

# -----------------------------------------------------------------

skirt_main_path = "~/SKIRT/release/SKIRTmain"
#fitskirt_main_path = "~/SKIRT/release/FitSKIRTmain"

variable_name = "PATH"

# SKIRT
value = skirt_main_path
export_command = "export " + variable_name + "=" + value + ":$" + variable_name
installation_commands["environment"].append(Map(line=export_command, comment="make the SKIRT executable detectable"))

# FitSKIRT
#value = fitskirt_main_path
#export_command = "export " + variable_name + "=" + value + ":$" + variable_name
#installation_commands["environment"].append(Map(line=export_command, comment="make the FitSKIRT executable detectable"))

# Set the path to the main SKIRT executable
#self.skirt_path = fs.join(skirt_main_path, "skirt")

# Success
#log.success("SKIRT was successfully downloaded")

# -----------------------------------------------------------------

cpp_path = modules.paths["cpp"]
mpi_path = modules.paths["mpi"]

if modules.names["cpp"] is not None:

    installation_commands["compilers"].append(Map(line="module load " + modules.names["cpp"], comment="Load the C++ compiler module"))

if modules.names["mpi"] is not None:

    if modules.names["mpi"] != modules.names["cpp"]:
        installation_commands["compilers"].append(Map(line="module load " + modules.names["mpi"], comment="Load the MPI module"))

# -----------------------------------------------------------------

installation_commands["build"].append(Map(line="cd " + skirt_repo_path, comment="Navigate to the SKIRT repository directory"))

qmake_path = modules.paths["qmake"]
make_make_command = qmake_path + " BuildSKIRT.pro -o ../release/Makefile CONFIG+=release"
nthreads = remote.cores_per_socket
make_command = "make -j " + str(nthreads) + " -w -C ../release"

# Debugging
log.debug("Make commands:")
log.debug(" 1) " + make_make_command)
log.debug(" 2) " + make_command)
#log.debug("in directory " + skirt_repo_path)

# Configure
#output = remote.execute(make_make_command, show_output=log.is_debug, cwd=skirt_repo_path)

# Overwrite the git version
#git_version_content = 'const char* git_version = " ' + git_version + ' " ;'
#git_version_path = fs.join(skirt_repo_path, "SKIRTmain", "git_version.h")
#write_command = 'echo "' + git_version_content + '" > ' + git_version_path
#remote.execute(write_command)

# Make
#output = remote.execute(make_command, show_output=log.is_debug, cwd=skirt_repo_path)

installation_commands["build"].append(Map(lines=[make_make_command, make_command], comment="build SKIRT"))

# -----------------------------------------------------------------

for step in ["qt", "directories", "clone", "build"]:

    if step not in installation_commands: continue

    print(step)
    print("")

    for command in installation_commands[step]:

        print("")

        if "comment" in command: print("# " + command.comment)

        if "line" in command: print(command.line)

        if "lines" in command:
            for line in command.lines: print(line)

        print("")

# -----------------------------------------------------------------
