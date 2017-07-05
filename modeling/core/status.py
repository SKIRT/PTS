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

class ModelingStatus(object):
        
    """
    This class ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        """

        # Set the modeling path
        self.modeling_path = modeling_path

        # Set the modeling type
        self.modeling_type = get_modeling_type(modeling_path)

        # Data
        self.status = []
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

        # Reset attributes
        self.status = []
        self.ntotal = 0
        self.nfinished = 0

        # Get the history
        history = load_modeling_history(self.modeling_path)
        history.clean()

        # Loop over the commands
        for command in commands_for_modeling_type(self.modeling_type):

            self.ntotal += 1

            if history.finished(command):
                self.status.append((command, "finished"))
                self.nfinished += 1
            elif command in history: self.status.append((command, "started"))
            else: self.status.append((command, "not_started"))

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

    def show(self):

        """
        This function ...
        :return:
        """

        print("")
        print(fmt.underlined + self.name + fmt.reset)
        print("")

        # Loop over the command
        for command, status in self.status:

            # Print status
            if status == "finished": print(fmt.green + " - " + command + ": finished" + fmt.reset)
            elif status == "started": print(fmt.yellow + " - " + command + ": started" + fmt.reset)
            elif status == "not_started": print(fmt.red + " - " + command + ": not started" + fmt.reset)
            else: print(fmt.red + " ! " + command + ": status not recognized: " + status + fmt.reset)

        print("")

        print(tostr(self.fraction_finished * 100, ndigits=3, round=True) + "% completed")
        print("")

# -----------------------------------------------------------------
