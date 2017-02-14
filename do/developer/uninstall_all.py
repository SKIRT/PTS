#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.uninstall_all Uninstall PTS and/or SKIRT on remote hosts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.tools.logging import log
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# Add optional settings
definition.add_optional("remotes", "string_list", "remote host(s) on which to uninstall", default=find_host_ids(), choices=find_host_ids())
definition.add_positional_optional("skirt_and_or_pts", "string_tuple", "SKIRT and/or PTS", default=["skirt", "pts"], choices=["skirt", "pts"])

# Add flags
definition.add_flag("conda", "also remove conda installation")
definition.add_flag("qt", "also remove Qt installation")

# Create setter
setter = InteractiveConfigurationSetter("uninstall_all", description="uninstall PTS and/or SKIRT on remote hosts", add_logging=True, add_cwd=True)
config = setter.run(definition, prompt_optional=True)

# -----------------------------------------------------------------

def uninstall_skirt(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    skirt_root_path = remote.skirt_root_path
    if not remote.is_directory(skirt_root_path): return

    # Inform the user
    log.info("Removing SKIRT from the remote host '" + remote.host_id + "' ...")

    # Remove the entire directory
    remote.remove_directory(skirt_root_path)

    # Remove lines from shell configuration file
    comment = "For SKIRT and FitSKIRT, added by PTS (Python Toolkit for SKIRT)"
    #terminal.remove_aliases_and_variables_with_comment(comment)
    remote.remove_aliases_and_variables_with_comment(comment)

# -----------------------------------------------------------------

def uninstall_pts(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    pts_root_path = remote.pts_root_path
    if not remote.is_directory(pts_root_path): return

    # Inform the user
    log.info("Removing PTS from the remote host '" + remote.host_id + "' ...")

    # Remove the entire directory
    remote.remove_directory(pts_root_path)

    # Remove lines from shell configuration file
    comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
    #terminal.remove_aliases_and_variables_with_comment(comment)
    remote.remove_aliases_and_variables_with_comment(comment)

# -----------------------------------------------------------------

def uninstall_conda(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Determine path of miniconda installation
    installation_path = fs.join(remote.home_directory, "miniconda")
    if not remote.is_directory(installation_path): return

    # Inform the user
    log.info("Removing Conda from the remote host '" + remote.host_id + "' ...")

    # Remove the directory
    remote.remove_directory(installation_path)

    # Remove lines from shell configuration file
    comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
    remote.remove_aliases_and_variables_with_comment(comment)

# -----------------------------------------------------------------

def uninstall_qt(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    pass

# -----------------------------------------------------------------

# Loop over the remotes
for host_id in config.remotes:

    # Create remote
    remote = Remote()
    if not remote.setup(host_id):
        log.warning("Remote host '" + host_id + "' is offline")
        continue

    # Locate SKIRT directory
    if "skirt" in config.skirt_and_or_pts: uninstall_skirt(remote)

    # Locate PTS directory
    if "pts" in config.skirt_and_or_pts: uninstall_pts(remote)

    # Uninstall conda
    if config.conda: uninstall_conda(remote)

    # Uninstall Qt
    if config.qt: uninstall_qt(remote)

# -----------------------------------------------------------------
