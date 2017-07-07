#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.clear_temp Clear the complete PTS temporary directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import time

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("remote", "string", "remote host on which to clear the temporary directory", choices=find_host_ids())
definition.add_optional("names", "string_list", "remove temporary directories matching these names (e.g. 'installation' matches 'installation_2017-05-03--08-45-44-410', 'installation_2017-05-03--13-52-28-371', etc. and also exact matches")

# Create setter
config = parse_arguments("clear_temp", definition)

# -----------------------------------------------------------------

# Remote
if config.remote is not None:

    # Create the remote
    remote = Remote(host_id=config.remote)

    # Check whether names are given
    if config.names is not None:

        # Loop over the directories in the temporary directory
        for path, name in remote.directories_in_path(remote.pts_temp_path, returns=["path", "name"]):

            # Contains timestamp
            if time.is_unique_name(name):

                # Without timestamp in list
                basename = time.get_name_from_unique_name(name)
                if basename in config.names: remote.remove_directory(path)

                # Full name in list
                if name in config.names: remote.remove_directory(path)

            # Full name in list
            elif name in config.names: remote.remove_directory(path)

    # Clear all temporary data
    else: remote.clear_directory(remote.pts_temp_path)

# Local
else: fs.clear_directory(introspection.pts_temp_dir)

# -----------------------------------------------------------------
