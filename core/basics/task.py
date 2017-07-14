#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.task Contains the Task class.

# -----------------------------------------------------------------

# Import standard modules
import importlib

# Import the relevant PTS classes and modules
from ..tools import serialization
from .configuration import Configuration
from ..tools import filesystem as fs
from pts.core.tools.utils import lazyproperty

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

        # The remote host ID and cluster name
        self.host_id = None
        self.cluster_name = None

        # The task ID and name
        self.id = None
        self.name = None

        # The command and the configuration in string format
        self.command = command
        self.config_string = config_string

        # The task file path
        self.path = None

        # The remote temporary PTS directory made for this task
        self.remote_temp_pts_path = None

        # Screen session name and remote screen output path
        self.screen_name = None
        self.remote_screen_output_path = None

        # Flag indicating whether the output of this task has been retrieved or not
        self.retrieved = False

        # Flag indicating whether the task has been analysed or not
        self.analysed = False

        # The path to the local and remote output directory
        self.local_output_path = None
        self.remote_output_path = None

        # Flag indicating whether we want to remove remote and local output (not used yet)
        self.remove_remote_output = True # AFTER RETRIEVAL
        self.remove_local_output = False # AFTER ANALYSIS

        # The paths to the task analysers
        self.analyser_paths = []

        # The analysis info
        self.analysis_info = dict()

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
    def remote_log_path(self):

        """
        This function ...
        :return:
        """

        return self.config.log_dir_path()

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

    def add_analyser(self, clspath):

        """
        This function ...
        :param clspath:
        :return:
        """

        self.analyser_paths.append(clspath)

    # -----------------------------------------------------------------

    @property
    def analyser_classes(self):

        """
        This function ...
        :return:
        """

        # The list of classes
        classes = []

        # Loop over the class paths
        for class_path in self.analyser_paths:

            module_path, class_name = class_path.rsplit('.', 1)

            # Get the class of the configurable of which an instance has to be created
            module = importlib.import_module(module_path)
            cls = getattr(module, class_name)

            # Add the class to the list of classes
            classes.append(cls)

        # Return the list of classes
        return classes

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Loop over the analyser classes that are defined for this task
        for analyser_class in self.analyser_classes:

            # Create an instance of the analyser class
            analyser = analyser_class.for_task(self)

            # Run the analyser
            analyser.run()

        # Set analysed flag to True
        self.analysed = True
        if self.path is not None: self.save()

        # Remove the local output if requested
        if self.remove_local_output: fs.remove_directory(self.local_output_path)

    # -----------------------------------------------------------------

    def remove_from_remote(self, remote, full=False):

        """
        This function ..
        :param remote:
        :param full:
        :return:
        """

        ## REMOVE REMOTE OUTPUT IF REQUESTED
        if self.remove_remote_output or full:

            # Remove the temporary PTS directory if it contains the output directory
            if remote.is_subdirectory(self.remote_output_path, self.remote_temp_pts_path) and remote.is_directory(self.remote_temp_pts_path): remote.remove_directory(self.remote_temp_pts_path)
            else:
                # Remove the output directory and the temporary directory seperately
                if remote.is_directory(self.remote_output_path): remote.remove_directory(self.remote_output_path)
                if remote.is_directory(self.remote_temp_pts_path): remote.remove_directory(self.remote_temp_pts_path)

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
