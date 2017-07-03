#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.loops Functions related to loops.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

def repeat(target, ntimes, **kwargs):

    """
    This function ...
    :param target:
    :param ntimes:
    :param kwargs:
    :return:
    """

    for _ in range(ntimes): target(**kwargs)

# -----------------------------------------------------------------

def repeat_check(target, ntimes, **kwargs):

    """
    This function ...
    :param target: 
    :param ntimes: 
    :param kwargs: 
    :return: 
    """

    for iteration in range(ntimes):
        if not target(**kwargs): raise RuntimeError("Iteration " + str(iteration) + " was not succesful")

# -----------------------------------------------------------------

def repeat_multiple(*targets, **kwargs):

    """
    This function ...
    :param targets:
    :param kwargs:
    :return:
    """

    ntimes = kwargs.pop("ntimes")

    # Repeat
    for _ in range(ntimes):
        for target in targets: target()

# -----------------------------------------------------------------
