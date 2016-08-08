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

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add positional optional
definition.add_positional_optional("remote", "string", "the name of the remote host for which to show/retrieve simulations and tasks", choices=find_host_ids())

# From remotesynchronizer.cfg:
# A dictionary that contains the ID's of simulations that have to be deleted from the synchronization in the future
#ids: None
# A list containing the name of each status for which corresponding simulations should be cleared
#statuses: None
# A list of simulation ID's to relaunch
#relaunch: None

# -----------------------------------------------------------------
