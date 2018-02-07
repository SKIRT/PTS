#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids, find_hosts
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import introspection

# -----------------------------------------------------------------

all_host_ids = find_host_ids()
all_hosts = find_hosts()
nhosts = len(all_hosts)
has_hosts = nhosts > 0

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The remote hosts
if has_hosts: definition.add_positional_optional("hosts", "host_list", "remote hosts", choices=all_host_ids, default=all_hosts)
else: definition.add_fixed("hosts", "remote hosts", [])

# Add optional
definition.add_optional("pts_repo_name", "string", "PTS repository name to deploy remotely", "origin", choices=introspection.pts_git_remotes())
definition.add_optional("skirt_repo_name", "string", "SKIRT repository name to deploy remotely", "origin", choices=introspection.skirt_git_remotes())

# Add flags
definition.add_flag("local", "also deploy locally", True)
definition.add_flag("skirt", "deploy SKIRT", True)
definition.add_flag("pts", "deploy PTS", True)
definition.add_flag("check", "check versions after deployment", True)
definition.add_flag("one_attempt", "only perform one attempt at connecting to a remote")

# Also update the dependencies
definition.add_flag("update_dependencies", "update the dependencies if possible", False)

# Add optional
definition.add_optional("pts_on", "string_list", "hosts on which PTS should be installed (None means all)", choices=find_host_ids())

# -----------------------------------------------------------------

# Dangerous stuff
definition.add_flag("clean", "do a completely clean install (remove existing installations) (use with care!)", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------
