#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.utils Contains remote utilities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class HostDownException(Exception):

    """
    This exception should be raised when connection to a host is not possible because it is temporarily down
    """

    def __init__(self, host_id):

        """
        This function ...
        :param host_id:
        """

        # Generate the message
        message = "The remote host '" + host_id + "' is down"

        # Set the message
        self.message = message

        # Call the constructor of the base class
        super(HostDownException, self).__init__(message)

# -----------------------------------------------------------------

class DetachedCalculation(Exception):
        
    """
    This exception should be raised when the calculation is continued remotely in detached mode, so that the
    algorithm should be pauzed/aborted and its state saved.
    """

    def __init__(self, cls, function_name):

        """
        The constructor ...
        :param cls:
        :param function_name:
        """

        # Generate the message
        message = "The calculation is continued in a detached session"

        # without this you may get DeprecationWarning
        self.message = message

        # Special attributes

        # The class and function name
        self.cls = cls
        self.function_name = function_name

        # Call the constructor of the base class
        super(DetachedCalculation, self).__init__(message)

# -----------------------------------------------------------------
