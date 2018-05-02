#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.status Contains the ModelingStatus class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
#from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import formatting as fmt
from ..component.component import load_modeling_history, get_modeling_type
from .steps import commands_for_modeling_type
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class ModelingStatus(list):
        
    """
    This class ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(ModelingStatus, self).__init__()

        # Set the modeling path
        self.modeling_path = modeling_path

        # Set the modeling type
        self.modeling_type = get_modeling_type(modeling_path)

        # Data
        self.ntotal = 0
        self.nfinished = 0

        # Refresh
        self.refresh()

    # -----------------------------------------------------------------

    def refresh(self):

        """
        This function ...
        :return:
        """

        # Clear
        self.clear()

        # Get the history
        history = load_modeling_history(self.modeling_path)
        history.clean()

        # Loop over the commands
        for command in commands_for_modeling_type(self.modeling_type):

            self.ntotal += 1

            if history.is_finished(command):
                self.append((command, "finished"))
                self.nfinished += 1
            elif command in history: self.append((command, "started"))
            else: self.append((command, "not_started"))

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the list
        del self[:]

        # Reset attributes
        self.ntotal = 0
        self.nfinished = 0

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def fraction_finished(self):

        """
        THis function ...
        :return:
        """

        return float(self.nfinished) / float(self.ntotal)

    # -----------------------------------------------------------------

    @property
    def fraction_unfinished(self):

        """
        This function ...
        :return:
        """

        return 1. - self.fraction_finished

    # -----------------------------------------------------------------

    @property
    def commands(self):

        """
        This function ...
        :return:
        """

        for command, status in self: yield command

    # -----------------------------------------------------------------

    @property
    def statuses(self):

        """
        This function ...
        :return:
        """

        for command, status in self: yield status

    # -----------------------------------------------------------------

    @property
    def colors(self):

        """
        This function ...
        :return:
        """

        #colours = []
        for command, status in self:
            if status == "finished": color = "green"
            elif status == "started": color = "yellow"
            elif status == "not_started": color = "red"
            else: color = "red"
            yield color
            #colours.append(color)
        #return colours

    # -----------------------------------------------------------------

    def show(self, name=True):

        """
        This function ...
        :return:
        """

        print("")
        if name:
            print(fmt.underlined + self.name + fmt.reset)
            print("")

        # Loop over the command
        for command, status in self:

            # Print status
            if status == "finished": print(fmt.green + " - " + fmt.bold + command + fmt.reset_bold + ": finished" + fmt.reset)
            elif status == "started": print(fmt.yellow + " - " + fmt.bold + command + fmt.reset_bold + ": started" + fmt.reset)
            elif status == "not_started": print(fmt.red + " - " + fmt.bold + command + fmt.reset_bold + ": not started" + fmt.reset)
            else: print(fmt.red + " ! " + command + ": status not recognized: " + status + fmt.reset)

        print("")

        print(tostr(self.fraction_finished * 100, ndigits=3, round=True) + "% completed")
        print("")

# -----------------------------------------------------------------
