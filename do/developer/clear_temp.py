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
from pts.core.tools.logging import setup_log
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("remote", "string", "remote host on which to clear the temporary directory", choices=find_host_ids())

# Create setter
setter = ArgumentConfigurationSetter("clear_temp")
config = setter.run(definition)

# -----------------------------------------------------------------

# Set logger
level = "DEBUG" if config.debug else "INFO"
log = setup_log(level)

# -----------------------------------------------------------------

# Remote
if config.remote is not None:

    remote = Remote()
    if not remote.setup(config.remote): raise RuntimeError("Could not connect to the remote")

    # Clear
    remote.clear_directory(remote.pts_temp_path)

# Local
else: fs.clear_directory(introspection.pts_temp_dir)

# -----------------------------------------------------------------
