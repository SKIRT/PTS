#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("ski", "file_path", "ski file to be used (if not specified, all ski files in the current directory will be used)")

# Add optional
definition.add_optional("remote", "string", "the remote host on which to launch the simulations", "nancy", choices=find_host_ids())

# -----------------------------------------------------------------
