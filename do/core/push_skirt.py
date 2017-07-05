#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.push_skirt Add all, commit and push to remote repositories, and pull again on remote systems.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import subprocess

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import introspection
from pts.core.prep.update import SKIRTUpdater

# -----------------------------------------------------------------

# Configuration definition
definition = ConfigurationDefinition()

# Required parameters
definition.add_required("message", "string", "the commit message")

# Optional parameters
definition.add_positional_optional("git_remotes", "string_list", "the remote(s) to commit to", default="origin")
definition.add_optional("host_id", "string", "the remote on which to pull the changes")

# -----------------------------------------------------------------

# Create configuration
config = parse_arguments("push_changes", definition, description="Add all, commit and push to remote repositories, and pull again on remote systems")

# -----------------------------------------------------------------

# Get the path to the local SKIRT repository
skirt_repo_path = introspection.skirt_repo_dir

# Stage all changes
subprocess.call(["git", "add", "--all"], cwd=skirt_repo_path)

# Make a commit with the specified message
subprocess.call(["git", "commit", "-m", config.message], cwd=skirt_repo_path)

# Push to the specified remotes
for remote in config.git_remotes: subprocess.call(["git", "push", remote, "master"], cwd=skirt_repo_path)

# If a remote host id is specified
if config.host_id is not None:

    # Create the SKIRT updater
    updater = SKIRTUpdater(config)

    # Run the updater
    updater.run()

# -----------------------------------------------------------------
