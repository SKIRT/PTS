#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.network Network-related functions

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from subprocess import call, check_output

# -----------------------------------------------------------------

def dns_ips():

    """
    This function ...
    :return:
    """

    # Alternative: cat /etc/resolv.conf ?

    command = "scutil --dns | grep 'nameserver\[[0-9]*\]' | sort | uniq"
    output = check_output(command, shell=True)

    entries = [entry.strip() for entry in output.split("\n") if entry.strip()]

    return [entry.split(" : ")[1] for entry in entries]

# -----------------------------------------------------------------

def dns_search_domains():

    """
    This function ...
    :return:
    """

    # Alternative: cat /etc/resolv.conf ?

    command = "scutil --dns | grep 'search domain\[[0-9]*\]' | sort | uniq"
    output = check_output(command, shell=True)

    entries = [entry.strip() for entry in output.split("\n") if entry.strip()]

    return [entry.split(" : ")[1] for entry in entries]

# -----------------------------------------------------------------
