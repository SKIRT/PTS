#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.task Contains the Task class.

# -----------------------------------------------------------------

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..tools import serialization
from .configuration import Configuration
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class Task(object):

    """
    This class ...
    """

    def __init__(self, command, config_string):

        """
        The constructor ...
        :param command:
        :param config_string:
        """

        # The task ID and name
        self.id = None
        self.name = None

        # The command and the configuration in string format
        self.command = command
        self.config_string = config_string

        # The task file path
        self.path = None

        # Screen session name and remote screen output path
        self.screen_name = None
        self.remote_screen_output_path = None

        # Flag indicating whether the output of this task has been retrieved or not
        self.retrieved = False

        # The path to the remote output directory
        self.remote_output_path = None

        # Flag indicating whether we want to keep remote output (not used yet)
        self.keep_remote_output = False

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the task object from file
        task = serialization.load(path)

        # Set the path of the task file
        task.path = path

        # Return the task object
        return task

    # -----------------------------------------------------------------

    @lazyproperty
    def config(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_string(self.config_string)

    # -----------------------------------------------------------------

    @lazyproperty
    def output_path(self):

        """
        This function ...
        :return:
        """

        return self.config.path

    # -----------------------------------------------------------------

    @property
    def screen_output_path(self):

        """
        This function ...
        :return:
        """

        if self.remote_screen_output_path is not None:
            name = "screenlog.0"
            path = fs.join(self.remote_screen_output_path, name)
            return path
        else: return None

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether a path is defined for the simulation file
        if self.path is None: raise RuntimeError("The task file does not exist yet")

        # Save to the original path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Set the new path
        self.path = path

        # Serialize and dump the task object
        serialization.dump(self, self.path, method="pickle")

# -----------------------------------------------------------------
