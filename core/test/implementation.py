#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.implementation Contains the TestImplementation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict
from abc import abstractmethod, ABCMeta

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import introspection
from ..basics.configurable import Configurable
from ..launch.pts import launch_local, launch_remote, RemoteInstance
from ..basics.log import log

# -----------------------------------------------------------------

tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

class TestImplementation(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        kwargs["unlisted"] = True
        super(TestImplementation, self).__init__(*args, **kwargs)

        # The test path
        self.path = None

        # The runnable components
        self.components = OrderedDict()
        self.input_dicts = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        self.path = kwargs.pop("path")

    # -----------------------------------------------------------------

    def run_command(self, command, remote=None):

        """
        This function ...
        :param command:
        :param remote:
        :return:
        """

        # Get command properties
        the_command = command.command
        description = command.description
        settings_dict = command.settings
        input_dict = command.input_dict
        cwd = command.cwd

        # Set the cwd if not specified
        if not command.cwd_specified: cwd = self.path

        # Determine the output path
        #output_path = fs.absolute_path(fs.join(self.path, cwd))
        output_path = fs.absolute_path(fs.absolute_or_in(cwd, self.path))

        #print(self.path)
        #print(cwd)
        #print(fs.absolute_or_in(cwd, self.path))

        # Launch locally or remotely
        if remote is not None:
            output = launch_remote(remote, the_command, settings_dict, input_dict)
            return RemoteInstance()
        else: return launch_local(the_command, settings_dict, input_dict, description=description, cwd=output_path, debug=log.is_debug)

# -----------------------------------------------------------------
