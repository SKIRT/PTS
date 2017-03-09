#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.introspection Contains some useful variables that store SKIRT and PTS installation directories
#  and provides functions for checking the presence and use of SKIRT and PTS dependencies.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import getpass
from os import devnull
import warnings
import sys
import imp
import inspect
import socket
import subprocess
from operator import itemgetter, methodcaller
from collections import defaultdict
from contextlib import contextmanager
from distutils.spawn import find_executable
from importlib import import_module

# Import the relevant PTS classes and modules
from . import filesystem as fs
from ..basics.map import Map

# -----------------------------------------------------------------

# Private repository links
private_skirt_ssh_link = "git@github.ugent.be:SKIRT/SKIRT.git"
private_skirt_https_link = "https://github.ugent.be/SKIRT/SKIRT.git"

# Public repository links
public_skirt_https_link = "https://github.com/SKIRT/SKIRT.git"

# -----------------------------------------------------------------

# Private repository links
private_pts_ssh_link = "git@github.ugent.be:SKIRT/PTS.git"
private_pts_https_link = "https://github.ugent.be/SKIRT/PTS.git"

# Public repository links
public_pts_https_link = "https://github.com/SKIRT/PTS.git"

# -----------------------------------------------------------------

possible_cpp_compilers = ["icc", "clang++", "c++", "cxx", "g++", "gcc", "gxx"]
possible_mpi_compilers = ["mpiicpc", "mpicxx", "mpiCC", "mpic++"]
possible_mpirun_names = ["mpirun", "orterun", "mpiexec"]

# -----------------------------------------------------------------

subprojects = ["core", "magic", "eagle", "evolve", "modeling", "dustpedia"]

# The path to the root PTS directory
#pts_root_dir = inspect.getfile(inspect.currentframe()).split("/pts")[0] # BUG: this won't give the correct result when the PTS root directory is named 'pts' in lowercase
pts_root_dir = fs.directory_of(fs.directory_of(fs.directory_of(fs.directory_of(inspect.getfile(inspect.currentframe())))))

# The path to the PTS package directory (PTS/pts)
pts_package_dir = fs.join(pts_root_dir, "pts")

# -----------------------------------------------------------------

# CREATION OF PTS DIRECTORIES EXTERNAL TO THE REPOSITORY

# The path to the PTS run directory (PTS/run)
pts_run_dir = fs.create_directory_in(pts_root_dir, "run")

# The path to the PTS user directory (PTS/user)
pts_user_dir = fs.create_directory_in(pts_root_dir, "user")

# The path to the PTS ext directory (PTS/ext)
pts_ext_dir = fs.create_directory_in(pts_root_dir, "ext")

# The path to the PTS temp directory (PTS/temp)
pts_temp_dir = fs.create_directory_in(pts_root_dir, "temp")

# The path to the PTS test directory (PTS/tests)
pts_tests_dir = fs.create_directory_in(pts_root_dir, "tests")

# The path to the PTS user/hosts directory
pts_user_hosts_dir = fs.create_directory_in(pts_user_dir, "hosts")

# The path to the PTS user/accounts directory
pts_user_accounts_dir = fs.create_directory_in(pts_user_dir, "accounts")

# The path to the PTS user/config directory
pts_user_config_dir = fs.create_directory_in(pts_user_dir, "config")

# -----------------------------------------------------------------

# Path to the versions file
pts_versions_path = fs.join(pts_package_dir, "versions.txt")

def get_constricted_versions():
    versions = dict()
    with open(pts_versions_path, 'r') as fh:
        for line in fh:
            line = line[:-1]
            name,version = line.split("==")
            versions[name] = version
    return versions

# -----------------------------------------------------------------

# Path to the names file
pts_names_path = fs.join(pts_package_dir, "names.txt")

def get_package_names():
    names = dict()
    with open(pts_names_path, 'r') as fh:
        for line in fh:
            line = line[:-1]
            import_name, real_name = line.split(":")
            names[import_name] = real_name
    return names

# -----------------------------------------------------------------

# path to installation commands file
pts_repositories_path = fs.join(pts_package_dir, "repositories.txt")

def get_package_repositories():
    urls = dict()
    with open(pts_repositories_path, 'r') as fh:
        for line in fh:
            line = line[:-1]
            import_name, url = line.split(" ")
            urls[import_name] = url
    return urls

# -----------------------------------------------------------------

# The path to the PTS do directory containing launchable scripts (PTS/pts/do)
pts_do_dir = fs.join(pts_package_dir, "do")

# The path to the PTS executable
pts_executable_path = fs.join(pts_do_dir, "__main__.py")

# The path to the main directory for a given PTS subproject
def pts_subproject_dir(subproject): return fs.join(pts_package_dir, subproject)

# The path to the 'config' directory for a given PTS subproject
def pts_config_dir(subproject): return fs.join(pts_package_dir, subproject, "config")

# The path to the 'dat' directory for a given PTS subproject
def pts_dat_dir(subproject): return fs.join(pts_package_dir, subproject, "dat")

# -----------------------------------------------------------------

def clear_temp():

    """
    This function ...
    :return:
    """

    fs.clear_directory(pts_temp_dir)

# -----------------------------------------------------------------

def clear_tests():

    """
    This function ...
    :return:
    """

    fs.clear_directory(pts_tests_dir)

# -----------------------------------------------------------------

def username():

    """
    This function ...
    :return:
    """

    return getpass.getuser()

# -----------------------------------------------------------------

def pts_version():

    """
    This function ...
    :return:
    """

    label = subprocess.check_output(["git", "describe", "--tags"], cwd=pts_package_dir)
    return label[:-1]

# -----------------------------------------------------------------

def pts_update_date():

    """
    This function ...
    :return:
    """

    # Try getting the time of the last pull
    try:

        if is_macos():

            command = "stat -f '%Sm' $(git rev-parse --show-toplevel)/.git/FETCH_HEAD"
            output = subprocess.check_output(command, cwd=pts_package_dir, shell=True)
            return output.split("\n")[0]

        else: return "update date unknown" # TODO: fix this

    # git pull has not been performed probably, this is a fresh clone
    except subprocess.CalledProcessError:

        # Try getting the time of the creation of the PTS files
        pts_init_file = fs.join(pts_package_dir, "__init__.py")
        return fs.creation_date(pts_init_file).strftime('%Y-%m-%d %H:%M:%S')

# -----------------------------------------------------------------

def has_account(service):

    """
    This function ...
    :param service:
    :return:
    """

    filepath = fs.join(pts_user_accounts_dir, service + ".txt")
    return fs.is_file(filepath)

# -----------------------------------------------------------------

def get_account(service):

    """
    This function ...
    :param service:
    :return:
    """

    import numpy as np
    filepath = fs.join(pts_user_accounts_dir, service + ".txt")
    username, password = np.loadtxt(filepath, dtype=str)
    return username, password

# -----------------------------------------------------------------

# The path to the SKIRT executable
skirt_path = find_executable("skirt")

# The path to the root SKIRT directory
skirt_root_dir = skirt_path.split("/release")[0] if skirt_path is not None else None

# The path to the SKIRT repository
skirt_repo_dir = fs.join(skirt_root_dir, "git") if skirt_path is not None else None

# The path to the SKIRT release directory
skirt_release_dir = fs.join(skirt_root_dir, "release") if skirt_path is not None else None

# The path to the SKIRT run directory
skirt_run_dir = fs.join(skirt_root_dir, "run") if skirt_path is not None else None

# -----------------------------------------------------------------

@contextmanager
def suppress_stdout():

    """
    This function ...
    :return:
    """

    with open(devnull, "w") as dvnll:
        old_stdout = sys.stdout
        sys.stdout = dvnll
        try:
            yield
        finally:
            sys.stdout = old_stdout

# -----------------------------------------------------------------

def python_version_long():

    """
    This function ...
    :return:
    """

    import pexpect

    # Launch interactive python session
    child = pexpect.spawn("python")
    child.expect(">>>")
    output = child.before.split("\r\n")

    # Close python again
    child.sendline("exit()")
    child.expect(pexpect.EOF)

    distribution_and_version = output[0].split("|")[0].strip()
    if "|" in output[0]:
        architecture = output[0].split("|")[1].strip()
        return distribution_and_version + " " + architecture
    else: return distribution_and_version

# -----------------------------------------------------------------

def python_version_short():

    """
    This function ...
    :return:
    """

    return str(sys.version_info.major) + "." + str(sys.version_info.minor) + "." + str(sys.version_info.micro)

# -----------------------------------------------------------------

def in_python_virtual_environment():

    """
    THis function ...
    :return:
    """

    if hasattr(sys, "real_prefix"): return True
    elif "conda" in sys.prefix.lower(): return True
    elif "canopy" in sys.prefix.lower(): return True
    elif "epd" in sys.prefix.lower(): return True
    elif "enthought" in sys.prefix.lower(): return True
    else: return False

# -----------------------------------------------------------------

def is_pts_developer():

    """
    This function ...
    :return:
    """

    from . import git

    if len(pts_git_remotes()) > 1: return True
    elif git.transform_to_simple_https(pts_git_remote_url(pts_git_remotes()[0])) not in [private_pts_https_link, public_pts_https_link]: return True
    else:

        output = subprocess.check_output(["git", "status"], cwd=pts_package_dir)

        for line in output:

            if "modified" in line: return True
            if "deleted" in line: return True

        return False

# -----------------------------------------------------------------

def is_skirt_developer():

    """
    This function ...
    :return:
    """

    from . import git

    if len(skirt_git_remotes()) > 1: return True
    elif git.transform_to_simple_https(skirt_git_remote_url(skirt_git_remotes()[0])) not in [private_skirt_https_link, public_skirt_https_link]: return True
    else:

        output = subprocess.check_output(["git", "status"], cwd=skirt_repo_dir)

        for line in output:

            if "modified" in line: return True
            if "deleted" in line: return True

        return False

# -----------------------------------------------------------------

def pts_git_branches():

    """
    This function ...
    :return:
    """

    args = ["git", "branch"]
    output = subprocess.check_output(args, cwd=pts_package_dir)

    branches = []
    for line in output.split("\n"):
        if line: branches.append(line.split("* ")[1])
    return branches

# -----------------------------------------------------------------

def pts_git_remotes():

    """
    This function ...
    :return:
    """

    args = ["git", "remote"]
    output = subprocess.check_output(args, cwd=pts_package_dir)

    remotes = []
    for line in output.split("\n"):
        if line: remotes.append(line)
    return remotes

# -----------------------------------------------------------------

def pts_git_remote_url(name):

    """
    This function ...
    :param name:
    :return:
    """

    args = ["git", "remote", "show", name]
    output = subprocess.check_output(args, cwd=pts_package_dir)

    for line in output.split("\n"):
        if "Fetch URL" in line: return line.split(": ")[1]

    raise RuntimeError("Remote '" + name + "' not found!")

# -----------------------------------------------------------------

def pts_git_official_remote():

    """
    This function ...
    :return:
    """

    remote_name = None

    for name in pts_git_remotes():

        url = pts_git_remote_url(name)

        if url == private_pts_https_link: remote_name = name
        elif url == private_pts_ssh_link: remote_name = name
        elif remote_name is None and url == public_pts_https_link: remote_name = name

    return remote_name

# -----------------------------------------------------------------

def skirt_git_branches():

    """
    This fucntion ...
    :return:
    """

    if skirt_repo_dir is None: raise RuntimeError("SKIRT is not installed (or not found) locally!")

    args = ["git", "branch"]
    output = subprocess.check_output(args, cwd=skirt_repo_dir)

    branches = []
    for line in output.split("\n"):
        if line: branches.append(line.split("* ")[1])
    return branches

# -----------------------------------------------------------------

def skirt_git_remotes():

    """
    This fucntion ...
    :return:
    """

    if skirt_repo_dir is None: raise RuntimeError("SKIRT is not installed (or not found) locally!")

    args = ["git", "remote"]
    output = subprocess.check_output(args, cwd=skirt_repo_dir)

    remotes = []
    for line in output.split("\n"):
        if line: remotes.append(line)
    return remotes

# -----------------------------------------------------------------

def skirt_git_remote_url(name):

    """
    This function ...
    :param name:
    :return:
    """

    if skirt_repo_dir is None: raise RuntimeError("SKIRT is not installed (or not found) locally!")

    args = ["git", "remote", "show", name]
    output = subprocess.check_output(args, cwd=skirt_repo_dir)

    for line in output.split("\n"):
        if "Fetch URL" in line: return line.split(": ")[1]

    raise RuntimeError("Remote '" + name + "' not found!")

# -----------------------------------------------------------------

def skirt_git_official_remote():

    """
    This function ...
    :return:
    """

    remote_name = None

    for name in skirt_git_remotes():

        url = skirt_git_remote_url(name)

        if url == private_skirt_https_link: remote_name = name
        elif url == private_skirt_ssh_link: remote_name = name
        elif remote_name is None and url == public_skirt_https_link: remote_name = name

    return remote_name

# -----------------------------------------------------------------

def is_macos():

    """
    This function ...
    :return:
    """

    return operating_system_short() == "Darwin"

# -----------------------------------------------------------------

def is_linux():

    """
    This function ...
    :return:
    """

    return operating_system_short() == "Linux"

# -----------------------------------------------------------------

def operating_system_short():

    """
    This function ...
    :return:
    """

    return subprocess.check_output(["uname", "-s"]).split("\n")[0]

# -----------------------------------------------------------------

def operating_system_long():

    """
    This function ...
    :return:
    """

    return subprocess.check_output(["uname", "-a"]).split("\n")[0]

# -----------------------------------------------------------------

def bashrc_path():

    """
    This function ...
    :return:
    """

    # Check if bashrc exists, if not create it if we are on Linux
    bashrc_path = fs.join(fs.home(), ".bashrc")
    if is_linux() and not fs.is_file(bashrc_path): fs.touch(bashrc_path)

    # Return the path
    return bashrc_path

# -----------------------------------------------------------------

def profile_path():

    """
    This function ...
    :return:
    """

    # Check whether .profile exists, create it if we are on MacOS
    profile_path = fs.join(fs.home(), ".profile")
    if is_macos() and not fs.is_file(profile_path): fs.touch(profile_path)

    # Return the path
    return profile_path

# -----------------------------------------------------------------

def bash_profile_path():

    """
    This function ...
    :return:
    """

    bash_profile_path = fs.join(fs.home(), ".bash_profile")
    return bash_profile_path

# -----------------------------------------------------------------

def shell_configuration_path():

    """
    This function ...
    :return:
    """

    if is_macos(): return profile_path()
    elif is_linux(): return bashrc_path()
    else: raise NotImplemented("System must be running MacOS or Linux")

# -----------------------------------------------------------------

def other_configuration_paths():

    """
    This function ...
    :return:
    """

    paths = []

    if is_macos():
        if fs.is_file(bashrc_path()): paths.append(bashrc_path())
        if fs.is_file(bash_profile_path()): paths.append(bash_profile_path())
    elif is_linux():
        if fs.is_file(bash_profile_path()): paths.append(bash_profile_path())
        if fs.is_file(profile_path()): paths.append(bash_profile_path())
    else: raise NotImplementedError("System must be running MacOS or Linux")

    return paths

# -----------------------------------------------------------------

def host_name():

    """
    This function ...
    :return:
    """

    return socket.gethostname()

# -----------------------------------------------------------------

def skirt_is_present():

    """
    This function ...
    :return:
    """

    return skirt_path is not None

# -----------------------------------------------------------------

def qmake_path():

    """
    This function ...
    :return:
    """

    try: output = subprocess.check_output(["which", "qmake"]).split("\n")[0]
    except subprocess.CalledProcessError: return None

    return output

# -----------------------------------------------------------------

def qmake_is_present():

    """
    This function ...
    :return:
    """

    return qmake_path() is not None

# -----------------------------------------------------------------

def has_qmake():

    """
    This function ...
    :return:
    """

    return qmake_is_present()

# -----------------------------------------------------------------

def qmake_version():

    """
    This function ...
    :return:
    """

    # Execute
    output = subprocess.check_output(qmake_path() + " --version", shell=True).split("\n")[0]
    return output

# -----------------------------------------------------------------

def has_conda():

    """
    This function ...
    :return:
    """

    return is_existing_executable("conda")

# -----------------------------------------------------------------

def conda_installation_path():

    """
    This function ...
    :return:
    """

    conda_exec_path = executable_path("conda")
    if fs.name(fs.directory_of(fs.directory_of(fs.directory_of(conda_exec_path)))) == "envs":
        return fs.directory_of(fs.directory_of(fs.directory_of(fs.directory_of(conda_exec_path))))
    else: return fs.directory_of(fs.directory_of(conda_exec_path))

# -----------------------------------------------------------------

def conda_main_executable_path():

    """
    This function ...
    :return:
    """

    path = fs.join(conda_installation_path(), "bin", "conda")
    if fs.is_file(path): return path
    else: return None

# -----------------------------------------------------------------

def pts_conformity_issues():

    """
    This function ...
    :return:
    """

    issues = []

    pts_root_dir_name = fs.name(pts_root_dir)
    if pts_root_dir_name != "PTS": issues.append("PTS root directory is not called 'PTS'")
    if fs.directory_of(pts_root_dir) != fs.home(): issues.append("PTS installation is not located in the home directory")

    if not has_conda(): issues.append("Conda executable cannot be located based on the PATH")
    else:

        installation_path = conda_installation_path()
        if fs.name(installation_path) != "miniconda": issues.append("Name of Conda installation directory is not 'miniconda'")
        if fs.directory_of(installation_path) != fs.home(): issues.append("Conda not installed in home directory")

    if len(issues) == 0: return None
    else: return ", ".join(issues)

# -----------------------------------------------------------------

def pts_installation_is_conform():

    """
    This function ...
    :return:
    """

    return pts_conformity_issues() is None

# -----------------------------------------------------------------

def skirt_conformity_issues():

    """
    This function ...
    :return:
    """

    issues = []

    skirt_root_dir_name = fs.name(skirt_root_dir)
    if skirt_root_dir_name != "SKIRT": issues.append("SKIRT root directory is not called 'SKIRT'")

    if fs.directory_of(skirt_root_dir) != fs.home(): issues.append("SKIRT installation is not located in the home directory")

    if len(issues) == 0: return None
    else: return ", ".join(issues)

# -----------------------------------------------------------------

def skirt_installation_is_conform():

    """
    This function ...
    """

    return skirt_conformity_issues() is None

# -----------------------------------------------------------------

def remote_host_ids():

    """
    This function ...
    :return:
    """

    # Search for files that define remote host configurations
    hosts_directory = fs.join(pts_user_dir, "hosts")
    if not fs.is_directory(hosts_directory): fs.create_directory(hosts_directory, recursive=True)

    # Initialize a list to contain the host ids
    ids = []

    # Loop over the configuration files in the hosts directory
    for name in fs.files_in_path(hosts_directory, extension="cfg", returns="name"):

        # Skip the template configuration file
        if name == "template": continue

        # Add the host name to the list of host ids
        ids.append(name)

    # Return the list of host ids
    return ids

# -----------------------------------------------------------------

def simulations_files_for_host(host_id):

    """
    This function checks whether there are simulation files corresponding to this host ID
    :param host_id:
    :return:
    """

    # Determine the path to the SKIRT run subdirectory for the specified host
    host_run_dir = fs.join(skirt_run_dir, host_id)

    # Return the list of simulation file paths corresponding to the specified host
    return fs.files_in_path(host_run_dir, extension="sim")

# -----------------------------------------------------------------

def is_existing_executable(name):

    """
    This function ...
    :param name:
    :return:
    """

    try:
        dvnll = open(devnull)
        subprocess.Popen(name, stdout=dvnll, stderr=dvnll).communicate()
        return True
    except: return False

# -----------------------------------------------------------------

def has_cpp_compiler():

    """
    This function ...
    :return:
    """

    for name in possible_cpp_compilers:
        if is_existing_executable(name): return True
    return False

# -----------------------------------------------------------------

def cpp_compiler_path():

    """
    This function ...
    :return:
    """

    for name in possible_cpp_compilers:
        path = find_executable(name)
        if path is not None: return path
    return None

# -----------------------------------------------------------------

def cpp_compiler_version():

    """
    This function ...
    :return:
    """

    # Execute
    output = subprocess.check_output(cpp_compiler_path() + " --version", shell=True).split("\n")
    return output[0]

# -----------------------------------------------------------------

def has_mpi():

    """
    This function ...
    :return:
    """

    for name in possible_mpirun_names:
        if is_existing_executable(name): return True
    return False

# -----------------------------------------------------------------

def has_mpi_compiler():

    """
    This function ...
    :return:
    """

    for name in possible_mpi_compilers:
        if is_existing_executable(name): return True
    return False

# -----------------------------------------------------------------

def mpi_compiler_path():

    """
    This function ...
    :return:
    """

    for name in possible_mpi_compilers:
        path = find_executable(name)
        if path is not None: return path
    return None

# -----------------------------------------------------------------

def mpi_compiler_version():

    """
    This function ...
    :return:
    """

    # Execute
    output = subprocess.check_output(mpi_compiler_path() + " --version", shell=True).split("\n")
    return output[0]

# -----------------------------------------------------------------

def mpi_version():

    """
    This function ...
    :return:
    """

    # Execute 'mpirun --version' and get the output
    process = subprocess.Popen(["mpirun", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0]

    # Return the relevant portion of the output
    for line in output:
        if "OpenMPI" in line: return "OpenMPI " + line.split(") ")[1]
        elif "Intel" in line: return "Intel MPI " + line.split("Version ")[1].split(" Build")[0]
    return " ".join(output) # unknown MPI version

# -----------------------------------------------------------------

def executable_path(name):

    """
    This function ...
    :param name:
    :return:
    """

    try:
        path = subprocess.check_output(["which", name]).split("\n")[0]
        return path
    except subprocess.CalledProcessError: return None

# -----------------------------------------------------------------

def has_skirt():

    """
    This function ...
    :return:
    """

    return executable_path("skirt") is not None

# -----------------------------------------------------------------

def skirt_version():

    """
    This function ...
    :return:
    """

    # Execute skirt with incorrect argument list and get its output
    process = subprocess.Popen([skirt_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0]

    # Return the relevant portion of the output
    for line in output.splitlines():
        if "SKIRT v" not in line: continue
        return line.split("Welcome to ")[1].split(" built on")[0] + ")"
    return None

# -----------------------------------------------------------------

def get_pip_versions():

    """
    This function ...
    :return:
    """

    packages = dict()

    # Launch the 'pip freeze' command and get the output
    output = subprocess.check_output(["pip", "list"])

    # Loop over the different package names and record the version number
    for entry in output.split("\n"):

        # Skip empty strings
        if not entry: continue

        # Get the package name and version
        name, version = entry.split(" (")
        version = version.split(")")[0]

        # Add them to the dictionary
        packages[name.lower()] = version

    # Return the dictionary
    return packages

# -----------------------------------------------------------------

@contextmanager
def ignore_site_packages_paths():

    """
    This function ...
    :return:
    """

    paths = sys.path[:]
    # remove working directory so that all
    # local imports fail
    if fs.cwd() in sys.path:
        sys.path.remove(fs.cwd())
    # remove all third-party paths
    # so that only stdlib imports will succeed
    sys.path = list(set(filter(
        None,
        filter(lambda i: all(('site-packages' not in i,
                              'python' in i or 'pypy' in i)),
               map(methodcaller('lower'), sys.path))
    )))
    yield
    sys.path = paths

# -----------------------------------------------------------------

def is_std_lib(module):

    """
    This function ...
    :param module:
    :return:
    """

    if not module:
        return False

    if module in sys.builtin_module_names:
        return True

    with ignore_site_packages_paths():
        imported_module = sys.modules.pop(module, None)
        try:
            with suppress_stdout(): import_module(module)
        except (ImportError, SystemExit): return False
        except AttributeError: return False # pip on peggy
        else:
            return True
        finally:
            if imported_module:
                sys.modules[module] = imported_module

# -----------------------------------------------------------------

def get_internal_imports():

    """
    This function ...
    :return:
    """

    imports = dict()

    # Loop over all python files in the PTS repository
    for filepath, filename in fs.files_in_path(pts_package_dir, extension="py", returns=["path", "name"], recursive=True):

        # Get the import lines for this file and set them
        import_lines = get_internal_imports_file(filepath)
        imports[filepath] = import_lines

    # Return the dictionary
    return imports

# -----------------------------------------------------------------

def get_internal_imports_file(filepath):

    """
    This function ...
    :return:
    """

    import_lines = []

    # Loop over the lines
    for line in fs.read_lines(filepath):

        # Strip whitespace
        line = line.strip()

        # Look for an 'import yyy' or 'from yyy import zzz' statement
        if not (line.startswith("import ") or (line.startswith("from ") and "import" in line)): continue

        # Set absolute import
        if "absolute_import" in line:
            has_absolute_import = True
            continue

        splitted = line.split()

        if splitted[1].startswith("pts") or splitted[1].startswith("."):
            import_lines.append(line)

    # Return the lines
    return import_lines

# -----------------------------------------------------------------

def get_internal_dependencies(debug=False):

    """
    This function ...
    :param debug:
    :return:
    """

    # Create an empty dictionary
    modules = defaultdict(set)

    # Get internal imports for all files
    import_lines = get_internal_imports()

    # Loop over the python files
    for filepath in import_lines:

        modules_file = set()

        internal_import_lines = import_lines[filepath]

        # Loop over the lines where something internally is imported
        for line in internal_import_lines:

            # Get the path to the modules that are being imported in the current line
            internal_module_paths = get_modules(line, filepath)

            # Add the modules
            for module_path in internal_module_paths: modules_file.add(module_path)

        # Set the modules
        modules[filepath] = modules_file

    # Return the dictionary of modules
    return modules

# -----------------------------------------------------------------

def get_all_dependencies():

    """
    This function ...
    :return:
    """

    # Create an empty dictionary to contain the required modules together with the places of use
    modules = defaultdict(set)

    # Recursively loop over all python files inside the PTS package directory
    for filepath, filename in fs.files_in_path(pts_package_dir, recursive=True, extension="py", returns=["path", "name"]):

        has_absolute_import = False
        modules_file = []

        # Read the lines of the script file
        for line in open(filepath, 'r'):

            # Strip whitespace and end-of-line characters
            line = line[:-1]
            line = line.strip()

            # If the "HIDE_DEPENDS" keyword is encountered, skip this file
            if "HIDE_DEPENDS" in line: break

            # Look for an 'import yyy' or 'from yyy import zzz' statement
            if line.startswith("import ") or (line.startswith("from ") and "import" in line):

                # Set absolute import
                if "absolute_import" in line: has_absolute_import = True

                # Get the name of the module(s)
                if line.startswith("import") and "," in line:
                    splitted = [part.strip().split(".")[0] for part in line.split(",")]
                    modules_line = [splitted[0].split(" ")[1]] + splitted[1:]
                else: modules_line = [line.split()[1].split(".")[0]]

                # Loop over the modules imported at this line
                for module in modules_line:

                    if skip_module(module): continue

                    # Add the module name to the list
                    if module: modules_file.append(module)

        # If absolute import is used, add all the imported modules (they must be external packages)
        if has_absolute_import:

            # Loop over found imported modules for this file
            for module_name in modules_file: modules[module_name].add(filepath)

        # If absolute import is not used, check whether some of the module names are just local (in the same directory)
        else:

            directory_path = fs.directory_of(filepath)
            local_module_names_in_dir = [filename[:-3] for filename in fs.files_in_path(directory_path, returns="name") if "__init__" not in filename]

            for module_name in modules_file:

                if module_name in local_module_names_in_dir: continue

                modules[module_name].add(filepath)

    # Return the modules
    return modules

# -----------------------------------------------------------------

def add_dependencies(dependencies, script_path, encountered_internal_modules, prefix="", debug=False):

    """
    This function ...
    :param dependencies:
    :param script_path:
    :param encountered_internal_modules:
    :param prefix:
    :param debug:
    :return:
    """

    # Skip files that are not python scripts
    if not script_path.endswith(".py"): raise ValueError("Not a valid script path")

    # Read the lines of the script file
    import_lines = []
    for line in open(script_path, 'r'):

        # Strip whitespace
        line = line.strip()

        # If the current line does not contain an 'import yyy' or 'from yyy import zzz' statement, skip it
        if not (line.startswith("import ") or (line.startswith("from ") and "import" in line)): continue

        import_lines.append(line)

    if debug:
        print("Import statements found in " + script_path + ":")
        for import_line in import_lines:
            print(" - " + import_line[:-1])

    # Loop over the lines where something is imported
    for line in import_lines:

        # Get the path to the modules that are being imported in the current line
        modules = get_modules(line, script_path)

        for module in modules:

            # Check if the imported module refers to a PTS module or an external package
            if module.startswith("/"): # a PTS module

                if module in encountered_internal_modules: continue
                else:
                    encountered_internal_modules.add(module)
                    add_dependencies(dependencies, module, encountered_internal_modules, prefix=prefix+"  ", debug=debug)

            else: dependencies[module].add(script_path)

# -----------------------------------------------------------------

def installed_python_packages():

    """
    This function ...
    :return:
    """

    # Initialize dictionary to contain the package names and version numbers
    packages = dict()

    import pip

    # Get all python distributions
    distributions = pip.get_installed_distributions()

    # Loop over the distributions
    for distribution in distributions:

        # Get name and version
        top_level_meta_data = list(distribution._get_metadata('top_level.txt'))
        import_name = top_level_meta_data[0] if len(top_level_meta_data) > 0 else distribution.project_name
        version = str(distribution.parsed_version)

        # possible other interesting properties of an entry in the distributions list:
        # .egg_name()
        # .as_requirement()
        # .parsed_version
        # .has_version()
        # .project_name
        # .py_version
        # .requires()

        # Add entry to the dictionary
        packages[import_name] = version

    # Return the dictionary
    return packages

# -----------------------------------------------------------------

def is_present_package(package):

    """
    This function ...
    :return:
    """

    try:
        with suppress_stdout(): imp.find_module(package)
        return True
    except (ImportError, SystemExit): return False

# -----------------------------------------------------------------

def get_modules(import_statement, script_path, return_unresolved=False, debug=False):

    """
    This function ...
    :param import_statement:
    :param script_path:
    :param return_unresolved:
    :param debug:
    :return:
    """

    unresolved = []

    multiple = None

    if "," in import_statement:

        if import_statement.startswith("from"):

            splitted = import_statement.split()

            imported = [splitted[3][:-1]]

            for more in splitted[4:]:

                if "," in more: more = more[:-1]
                imported.append(more)

        else:

            #print(import_statement)

            splitted = import_statement.split(",")
            firstsplitted = splitted[0].split()

            if "." in firstsplitted[1]: fr, wh = firstsplitted[1].split(".")
            else: fr, wh = firstsplitted[1].strip(), None

            if wh is not None: imported = [wh.strip()]
            else: imported = []

            #imported = [firstsplitted]

            for more in splitted[1:]:

                if "." in more: fr2, wh2 = more.split(".")
                else: fr2, wh2 = more.strip(), None
                fr2 = fr2.strip()
                #print(fr2)
                #print(wh2)
                #print(splitted)
                #print(imported)
                if fr2 != fr:
                    if multiple is None:
                        multiple = []
                        multiple.append(fr)
                    multiple.append(fr2)
                    #pass #raise RuntimeError("Cannot proceed: " + fr2 + " is not equal to " + fr)
                if wh2 is not None: imported.append(wh2.strip())

            splitted = ["from", fr]

    else:

        splitted = import_statement.split()

        #print(splitted)
        if len(splitted) <= 2: imported = []
        else: imported = [splitted[3]]

    #print(multiple)
    #print(splitted)
    #print(imported)

    if multiple is None:

        #which = []
        which = defaultdict(set)

        # Check if this line denotes a relative import statement
        if splitted[1].startswith("."):

            after_dots = splitted[1].lstrip(".")

            number_of_dots = len(splitted[1]) - len(after_dots)

            # Determine the path to the PTS subpackage
            subpackage_dir = script_path
            for i in range(number_of_dots):
                subpackage_dir = fs.directory_of(subpackage_dir)

            subpackage_name = after_dots.split(".")[0]
            subpackage_path = fs.join(subpackage_dir, after_dots.replace(".", "/"))

            for name in imported:

                #print(subpackage_path)

                #if debug: print(name)
                module_path = which_module(subpackage_path, name)
                if debug: print(subpackage_path, name, ":", module_path)
                if module_path is not None:
                    #print(fs.strip_extension(fs.name(module_path)), name)
                    if fs.strip_extension(fs.name(module_path)) == name: which[module_path] = None
                    else:
                        if name == "*": name = None
                        if name == "__version__": name = None
                        which[module_path].add(name)
                else: unresolved.append((subpackage_path, name))

        # Absolute import of a pts class or module
        elif splitted[1].startswith("pts"):

            parts = splitted[1].split(".")[1:]

            subpackage_dir = pts_package_dir
            for part in parts:
                subpackage_dir = fs.join(subpackage_dir, part)

            for name in imported:
                #print(name)
                module_path = which_module(subpackage_dir, name)
                #print(module_path)
                if debug: print(subpackage_dir, name, ":", module_path)
                if module_path is not None:
                    if fs.strip_extension(fs.name(module_path)) == name: which[module_path] = None
                    else: which[module_path].add(name)
                else: unresolved.append((subpackage_dir, name))

        # MPL toolkits
        #elif splitted[1].startswith("mpl_toolkits"): pass # skip mpl_toolkits
        # Pylab
        #elif splitted[1].startswith("pylab"): pass # skip pylab

        # Skip modules
        elif skip_module(splitted[1]): pass

        # External module
        else:

            # Get the name of the module
            module = splitted[1].split(".")[0]

            imported_names = []
            if len(splitted) > 2 and splitted[2] == "import":
                for name in imported: imported_names.append(name)
            else: imported_names = ["__init__"]

            for name in imported_names: which[module].add(name)

    else:

        which = defaultdict(set)
        for name in multiple: which[name].add("__init__")

    # Return the list of modules
    if return_unresolved: return which, unresolved
    else: return which

# -----------------------------------------------------------------

def which_module(subpackage, name):

    """
    This function ...
    :param subpackage:
    :param name:
    :return:
    """

    # Find file
    if name.islower() and fs.is_file(fs.join(subpackage, name + ".py")): return fs.join(subpackage, name + ".py")

    #
    elif fs.is_file(subpackage + ".py"): return subpackage + ".py"

    # Init file
    elif fs.is_file(fs.join(subpackage, "__init__.py")): return fs.join(subpackage, "__init__.py")

    else: #raise ValueError("Don't know how to get further with " + subpackage + " and " + name)

        warnings.warn("Incorrect import: " + name + " from " + subpackage)
        return None

# -----------------------------------------------------------------

def get_scripts():

    """
    This function ...
    :return:
    """

    # Loop over the directories within the 'do' subpackage
    scripts = []
    for item_path, item in fs.directories_in_path(pts_do_dir, returns=["path", "name"]):

        # Loop over the files in the directory
        for filepath, filename in fs.files_in_path(item_path, not_contains="__init__", extension="py", returns=["path", "name"]): scripts.append((item, filename + ".py"))

    # Return the sorted list of script names
    return sorted(scripts, key=itemgetter(1))

# -----------------------------------------------------------------

def find_matches_scripts(name, scripts):

    """
    This function ...
    :param name:
    :param scripts:
    :return:
    """

    # Get a list of the script names that match the first command line argument, if there is one
    if "/" in name:

        matches = []
        dir_name = name.split("/")[0]
        script_name = name.split("/")[1]

        # Loop over all found items
        for item in scripts:
            if item[0] == dir_name and item[1].startswith(script_name): matches.append(item)

        return matches

    # Return the list of matching scripts
    elif name is not None: return filter(lambda item: item[1].startswith(name), scripts)
    else: return []

# -----------------------------------------------------------------

def find_matches_tables(name, tables):

    """
    This function ...
    :param name:
    :param tables:
    :return:
    """

    # ...
    if "/" in name:

        matches = []
        dir_name = name.split("/")[0]
        script_name = name.split("/")[1]

        for subproject in tables:

            if dir_name != subproject: continue

            table = tables[subproject]

            for i in range(len(table["Command"])):

                command_name = table["Command"][i]
                if command_name.startswith("*"): command_name = command_name[1:]

                if command_name.startswith(script_name): matches.append((subproject, i))

        return matches

    elif name is not None:

        matches = []

        for subproject in tables:

            table = tables[subproject]

            for i in range(len(table["Command"])):

                command_name = table["Command"][i]
                if command_name.startswith("*"): command_name = command_name[1:]

                if command_name.startswith(name): matches.append((subproject, i))

        return matches

    else: return []

# -----------------------------------------------------------------

def get_arguments_tables():

    """
    This function ...
    :return:
    """

    #import numpy as np

    tables = dict()

    # Loop over the subprojects
    for subproject in subprojects:

        table_path = fs.join(pts_subproject_dir(subproject), "commands.dat")
        if not fs.is_file(table_path): continue

        # Get the columns
        #commands, configuration, where, method, description = np.genfromtxt(table_path, delimiter=" | ", dtype=str, unpack=True)
        commands = []
        configuration = []
        where = []
        method = []
        description = []

        # Numpy-less implementation
        with open(table_path, 'r') as table:
            for line in table:
                if line.startswith("#"): continue
                line = line[:-1]
                if not line: continue
                splitted = line.split(" | ")
                commands.append(splitted[0])
                configuration.append(splitted[1])
                where.append(splitted[2])
                method.append(splitted[3])
                description.append(splitted[4])

        # Fix
        if isinstance(commands, basestring):
            commands = [commands]
            configuration = [configuration]
            where = [where]
            method = [method]
            description = [description]

        # Table
        table = {"Command": commands, "Configuration": configuration, "Path": where, "Configuration method": method, "Description": description}
        tables[subproject] = table

    # Return the tables
    return tables

# -----------------------------------------------------------------

def skip_module(name, path=None):

    """
    This function ...
    :return:
    """

    if name == "core": return True
    if name == "__future__": return True
    if name == "__init__": return True
    if name == "__main__": return True
    if name == "runner": return True
    if name == "run_queue": return True
    if name == "enable_qch_mathjax": return True
    if name == "mpl_toolkits": return True
    if name == "pylab": return True
    if name == "pts": return True

    # If path is known
    if path is not None:

        if path.endswith("eagle/config.py"): return True
        if path.endswith("eagle/collections.py"): return True
        if path.endswith("eagle/database.py"): return True
        if path.endswith("eagle/galaxy.py"): return True
        if path.endswith("eagle/plotresults.py"): return True
        if path.endswith("eagle/runner.py"): return True
        if path.endswith("eagle/scheduler.py"): return True
        if path.endswith("eagle/skirtrun.py"): return True
        if name == "fit2BB_Md": return True

    return False

# -----------------------------------------------------------------

def all_modules():

    """
    This function ...
    :return:
    """

    # Loop over all module files
    for path, name in fs.files_in_path(pts_package_dir, extension="py", recursive=True, returns=["path", "name"]):

        # Skip?
        if skip_module(name, path): continue

        # Get subproject
        subproject = path.split("pts/")[1].split("/")[0]

        is_config = path.split(subproject + "/")[1].split("/")[0] == "config"

        #if subproject != "do" and not is_config:

        if subproject == "do": continue
        if is_config: continue

        # Load the module, get number of functions
        relpath = "pts." + path.split("pts/")[1].replace("/", ".")[:-3]
        # import_name = relpath.split(".")[-1]
        # relpath = relpath.split("." + import_name)[0]
        # print(import_name)

        try: module = import_module(relpath)
        except ImportError: module = None

        # Return module path and the module itself (if it could be imported)
        yield path, subproject, module

# -----------------------------------------------------------------

def classes_in_module(module):

    """
    This function ...
    :param module:
    :return:
    """

    md = module.__dict__
    return [md[c] for c in md if (isinstance(md[c], type) and md[c].__module__ == module.__name__)]

# -----------------------------------------------------------------

def base_classes(cls):

    """
    This function ...
    :param cls:
    :return:
    """

    return inspect.getmro(cls)[1:]

# -----------------------------------------------------------------

def is_leaf_class(cls, hierarchy):

    """
    This function ...
    :param cls:
    :param hierarchy:
    :return:
    """

    for other_cls in hierarchy:

        base_class_list = hierarchy[other_cls]

        if cls in base_class_list: return False

    return True

# -----------------------------------------------------------------

def all_classes():

    """
    This function ...
    :return:
    """

    # Loop over all modules
    for path, subproject, module in all_modules():
        if module is None: continue
        #for name, obj in inspect.getmembers(module):
        #    if inspect.isclass(obj): yield obj

        for cls in classes_in_module(module): yield cls

# -----------------------------------------------------------------

def all_configurable_classes():

    """
    This function ...
    :return:
    """

    from ..basics.configurable import Configurable

    for cls in all_classes():
        if not issubclass(cls, Configurable): continue
        yield cls

# -----------------------------------------------------------------

def all_concrete_configurable_classes():

    """
    This function ...
    :return:
    """

    classes = all_configurable_classes()

    # Construct hierarchy
    hierarchy = dict()
    for cls in classes:
        hierarchy[cls] = base_classes(cls)

    for cls in all_configurable_classes():

        #print(cls.__name__, inspect.isabstract(cls))

        #if cls.__metaclass__ == ABCMeta: continue
        #if inspect.isabstract(cls): continue
        #yield cls

        if is_leaf_class(cls, hierarchy): yield cls

# -----------------------------------------------------------------

def resolve_command_configurable(command):

    """
    This function ...
    :param command:
    :return:
    """

    tables = get_arguments_tables()
    return resolve_command_tables(command, tables)

# -----------------------------------------------------------------

def resolve_command_tables(command, tables):

    """
    This function ...
    :param command:
    :param tables:
    :return:
    """

    # Find matches
    script_name = command
    table_matches = find_matches_tables(script_name, tables)

    # No match or ambigious command
    if len(table_matches) == 0: raise ValueError("The command '" + command + "' could not be resolved")
    elif len(table_matches) > 1: raise ValueError("The command is ambigious")

    # Return the match
    subproject, index = table_matches[0]
    return resolve_from_match(subproject, tables[subproject], index)

# -----------------------------------------------------------------

def resolve_from_match(subproject, table, index):

    """
    This function ...
    :param subproject:
    :param table:
    :param index:
    :return:
    """

    # Get properties
    command_name = table["Command"][index]
    hidden = False
    if command_name.startswith("*"):
        hidden = True
        command_name = command_name[1:]
    command_description = table["Description"][index]
    class_path_relative = table["Path"][index]
    class_path = "pts." + subproject + "." + class_path_relative
    module_path, class_name = class_path.rsplit('.', 1)

    configuration_method_table = table["Configuration method"][index]
    # subproject_path = pts_subproject_dir(subproject)

    # Determine the configuration module path
    configuration_name = table["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    # Return
    return Map(subproject=subproject, command_name=command_name, hidden=hidden, description=command_description,
               module_path=module_path, class_name=class_name, configuration_method=configuration_method_table,
               configuration_module_path=configuration_module_path)

# -----------------------------------------------------------------

def get_class(module_path, class_name):

    """
    This function ...
    :param module_path
    :param class_name:
    :return:
    """

    # Get the class of the configurable of which an instance has to be created
    module = import_module(module_path)
    try: cls = getattr(module, class_name)
    except AttributeError:
        raise Exception("The class '" + class_name + "' could not be found in the module '" + module_path + "'")

    # Return the class
    return cls

# -----------------------------------------------------------------

def get_configuration_definition(configuration_module_path):

    """
    This function ...
    :param configuration_module_path:
    :return:
    """

    # Get the configuration module
    try: configuration_module = import_module(configuration_module_path)
    except ImportError: raise RuntimeError("The configuration module '" + configuration_module_path + "' was not found")
    definition = getattr(configuration_module, "definition")

    # Return the definition
    return definition

# -----------------------------------------------------------------

def get_configuration_definition_pts_not_yet_in_pythonpath(configuration_module_path):

    """
    This function ...
    :param configuration_module_path:
    :return:
    """

    name = configuration_module_path.split(".")[-1]
    relative_import_path = ".".join(configuration_module_path.split("pts.")[1].split(".")[:-1])
    filepath = fs.join(pts_package_dir, relative_import_path.replace(".", "/"), name + ".py")

    # Load the test module
    configuration_module = imp.load_source(name, filepath)

    # Get the definition
    return configuration_module.definition

# -----------------------------------------------------------------

def try_importing_module(path, in_globals=False):

    """
    This function ...
    :param path:
    :param in_globals:
    :return:
    """

    import imp
    import importlib

    # Determine the name of the module
    name = path.split(".")[-1]

    # Internal
    if "pts" in path:

        # Determine relative import path
        rel_import_path = path.split("pts.")[1]

        #print(path)

        # Determine file path
        file_path = fs.join(pts_package_dir, rel_import_path.replace(".", "/") + ".py")
        #module_path = fs.join(pts_package_dir, rel_import_path.replace(".", "/"))

        #print(module_path)

        # Load the module
        #print(file_path)
        #module = imp.load_source(name, file_path)

        module = importlib.import_module(path)

    # External
    else:
        if "." in path:
            #import importlib
            #splitted = path.split(".")
            #file, filename, (suffix, mode, type) = imp.find_module(splitted[0])
            #module_path = fs.join(filename, *(splitted[1:] + ["__init__.py"]))
            #module_path = fs.join(filename, *splitted[1:])
            #module = imp.load_module(splitted[-1], None, module_path, ('', '', 5)) # 5 means pacakge directory (so no __init__.py)
            module = importlib.import_module(path)
        else:
            file, filename, (suffix, mode, type) = imp.find_module(path)
            module = imp.load_module(name, file, filename, (suffix, mode, type))

    # Add the module to the globals
    if in_globals: globals()[name] = module

    # Return the module
    return module

# -----------------------------------------------------------------

def try_importing_class(name, path, in_globals=False):

    """
    This function ...
    :param name:
    :param path:
    :param in_globals:
    :return:
    """

    # Import the module
    module = try_importing_module(path)

    # Get the class
    cls = getattr(module, name)

    #print(cls)

    # Add the class to the globals
    if in_globals: globals()[name] = cls

    #print(globals()[name])

    # Return the class
    return cls

# -----------------------------------------------------------------

def lazy_isinstance(instance, classname, path):

    """
    This function ...
    :param instance:
    :param classname:
    :param path:
    :return:
    """

    cls = try_importing_class(classname, path)
    return isinstance(instance, cls)

# -----------------------------------------------------------------

def lazy_call(function_name, module_name, *args, **kwargs):

    """
    This function ...
    :param function_name:
    :param module_name:
    :param args:
    :param kwargs:
    :return:
    """

    module = try_importing_module(module_name)
    func = getattr(module, function_name)
    return func(*args, **kwargs)

# -----------------------------------------------------------------
