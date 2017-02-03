#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.jobscript_template Generate a job script template.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.remote.modules import Modules
from pts.core.tools.logging import setup_log
from pts.core.tools import formatting as fmt
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import introspection
from pts.core.tools import git

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host ID", choices=find_host_ids())

# Create the configuration
setter = ArgumentConfigurationSetter("installation_commands")
config = setter.run(definition)

# -----------------------------------------------------------------

# Setup log
if config.debug: log = setup_log("DEBUG")
else: log = setup_log()

# -----------------------------------------------------------------



# -----------------------------------------------------------------
