#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.find_files Find files containing a certain string in the current working directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Remote
definition.add_positional_optional("remote", "string", "remote host ID", choices=find_host_ids())
definition.add_positional_optional("remote_path", "string", "path of remote directory (absolute or relative w.r.t. home directory)")

# Add optional settings
definition.add_optional("contains", "string", "a string that should be contained in the file names")
definition.add_optional("not_contains", "string", "a string that should not be contained in the file names")
definition.add_optional("extension", "string", "the file extension")
definition.add_optional("exact_name", "string", "exact filename (for recursive search)")
definition.add_optional("exact_not_name", "string", "not files with exactly this name")
definition.add_optional("directory_not_contains", "string_list", "directory_not_contains")
definition.add_optional("directory_exact_not_name", "string_list", "directory_exact_not_name")

# Add flags
definition.add_flag("recursive", "search recursively", False)
definition.add_flag("full", "show the full paths", False)
definition.add_flag("remove", "remove the found files", False)

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("find_files", definition, description="Find files containing a certain string in the current working directory")

# -----------------------------------------------------------------

# REMOTELY
if config.remote is not None:

    # Create remote
    remote = Remote(host_id=config.remote)

    # Determine path
    if config.remote_path is not None: find_path = remote.absolute_or_in_home(config.remote_path)
    else: find_path = remote.home_directory

    #print(find_path)

    #print(remote.items_in_path(find_path, recursive=True))

    # Loop over the files
    paths = remote.files_in_path(find_path, contains=config.contains, not_contains=config.not_contains, extension=config.extension, recursive=config.recursive, exact_name=config.exact_name, exact_not_name=config.exact_not_name)
    nfiles = len(paths)

    if nfiles == 0: log.warning("No files found")
    else:
        log.info(str(nfiles) + " files found")
        for path in paths:
            if config.full: print(path)
            else: print(path.split(find_path)[1])
            if config.remove: remote.remove_file(path)

# LOCALLY
else:

    # Determine path
    find_path = fs.cwd()

    # Loop over the files
    paths = fs.files_in_path(find_path, contains=config.contains, extension=config.extension, recursive=config.recursive,
                             exact_name=config.exact_name, exact_not_name=config.exact_not_name,
                             directory_not_contains=config.directory_not_contains, directory_exact_not_name=config.directory_exact_not_name)
    nfiles = len(paths)

    if nfiles == 0: log.warning("No files found")
    else:
        log.info(str(nfiles) + " files found")
        for path in paths:
            if config.full: print(path)
            else: print(path.split(find_path)[1])
            if config.remove: fs.remove_file(path)

# -----------------------------------------------------------------
