#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add required arguments
definition.add_required("ski", "file_path", "name/path of the ski file")
definition.add_required("fski", "file_path", "name/path of the fski file")

# Input and output
definition.add_optional("input", "directory_path", "input directory for the simulation(s)", letter="i")
definition.add_optional("output", "directory_path", "output directory for the simulation(s)", fs.cwd(), letter="o", convert_default=True)

# Add positional arguments
definition.add_positional_optional("remote", "string", "the remote host on which to run the simulation (if none is specified, the simulation is run locally", choices=find_host_ids())

# -----------------------------------------------------------------
