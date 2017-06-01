#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.setup Contains a function that performs a setup.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

def setup(command_name, cwd):

    """
    This function ...
    :param command_name:
    :param cwd:
    :return: 
    """

    # When launching a seperate modeling command, add entry to history
    if command_name == "setup" or command_name == "model": return

    # Load the history
    from .component.component import load_modeling_history
    history = load_modeling_history(cwd)

    # Add entry
    history.add_entry_and_save(command_name)

# -----------------------------------------------------------------

def finish(command_name, cwd):

    """
    This function ...
    :param command_name: 
    :param cwd:
    :return: 
    """

    if command_name == "setup" or command_name == "model": return

    # Load the history
    from .component.component import load_modeling_history
    history = load_modeling_history(cwd)

    # Mark end
    history.mark_end_and_save()

# -----------------------------------------------------------------
