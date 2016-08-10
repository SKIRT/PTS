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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("contains", "string", "a string that should be contained in the file names")

# Add optional
definition.add_optional("extension", "string", "the file extension")

# Add flag
definition.add_flag("recursive", "search recursively")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("find_files", "Find files containing a certain string in the current working directory")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(config.path, time.unique_name("status") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting find_files ...")

# -----------------------------------------------------------------

for path in fs.files_in_path(fs.cwd(), contains=config.contains, extension=config.extension, recursive=config.recursive): 
    print(path.split(fs.cwd())[1])

# -----------------------------------------------------------------
