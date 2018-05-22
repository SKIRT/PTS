#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.fill Contains the SkiFiller class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import subprocess
import difflib

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..simulation.skifile import SkiFile
from ..tools import formatting as fmt
from ..basics.log import log
from ..simulation.logfile import LogFile
from ..basics.map import Map
from ..tools import introspection, time

# -----------------------------------------------------------------

class SkiFiller(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SkiFiller, self).__init__(*args, **kwargs)

        # The ski file
        self.ski = None

        # The temporary directory path
        self.temp_path = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Run SKIRT to generate the parameters file
        self.generate_parameters()

        # Show the settings that were added
        self.show()

        # Writing
        self.write()

        # Remove temporary files
        self.cleanup()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "ski" in kwargs: self.config.ski_path = kwargs.pop("ski")
        elif "ski_path" in kwargs: self.config.ski_path = kwargs.pop("ski_path")

        if isinstance(self.config.ski_path, SkiFile): self.ski = self.config.ski_path
        elif isinstance(self.config.ski_path, basestring): self.ski = SkiFile(self.config.ski_path)
        else: raise ValueError("Invalid option for 'ski_path': must be ski file or ski file path")

        # Create temporary directory and set the path
        self.temp_path = introspection.create_temp_dir(time.unique_name("fill"))

    # -----------------------------------------------------------------

    def generate_parameters(self):

        """
        This function ...
        :return:
        """

        # Save the ski file in the temporary directory
        ski_path = fs.join(self.temp_path, self.ski.prefix + ".ski")
        self.ski.saveto(ski_path, update_path=False)

        # Create the command
        command = ["skirt", "-p", ski_path, "-o", self.temp_path]

        # Run SKIRT
        if log.is_debug: subprocess.call(command)
        else: subprocess.call(command, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))

        # Load the parameters file
        parameters_path = fs.join(self.temp_path, self.ski.prefix + "_parameters.xml")
        parameters = SkiFile(parameters_path)

        # Get the differences between the files

        a = self.ski.to_lines()
        b = parameters.to_lines()

        for line in difflib.context_diff(a, b):
            print(line)

        print(parameters)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print(fmt.green + "")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the
        #if

        pass

    # -----------------------------------------------------------------

    def cleanup(self):

        """
        This fucntion ...
        :return:
        """

        # Remove the temporary directory
        fs.remove_directory(self.temp_path)

# -----------------------------------------------------------------
