#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.history Contains the ModelingHistory class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import time, tables
from .steps import single_commands

# -----------------------------------------------------------------

class ModelingHistory(SmartTable):
    
    """
    This class...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Command"] = (str, None, "name of the modeling command")
    _column_info["Start time"] = (str, None, "timestamp for start of command")
    _column_info["End time"] = (str, None, "timestamp for end of command")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if len(args) > 0: from_astropy = True
        else: from_astropy = False

        # Call the constructor of the base class
        super(ModelingHistory, self).__init__(*args, **kwargs)

        # Add column info
        if not from_astropy: self.add_all_column_info(self._column_info)

        # Clean
        # DOESN'T WORK HERE: TypeError: could not convert reader output to ModelingHistory class.
        # SO.. CLEAN() HAS TO BE CALLED EXPLICITLY
        #self.clean()

    # -----------------------------------------------------------------

    def clean(self):

        """
        This function ...
        :return:
        """

        # Setup the table
        self._setup()

        # Loop over all commands
        for command in single_commands:

            indices = tables.find_indices(self, command)

            # More than one occurence?
            if len(indices) > 0:

                keep = None

                # Only keep the one that finished?
                # Loop from last
                for index in reversed(indices):
                    if not self["End time"].mask[index]:
                        keep = index
                        break

                # Finished entry was found, remove all other except this
                if keep is not None:
                    to_remove = copy.copy(indices)
                    to_remove.remove(keep)
                    self.remove_rows(to_remove)

                # No finished found, only keep last unfinished
                else:
                    keep = indices[-1]
                    to_remove = copy.copy(indices)
                    to_remove.remove(keep)
                    self.remove_rows(to_remove)

    # -----------------------------------------------------------------

    def add_entry(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Set the values
        values = [command, time.timestamp(), None]

        # Add a row to the table
        self.add_row(values)

    # -----------------------------------------------------------------

    def add_entry_and_save(self, command):

        """
        Thi function ...
        :param command:
        :return:
        """

        self.add_entry(command)
        self.save()

    # -----------------------------------------------------------------

    def remove_entry(self, command):

        """
        This function removes
        :param command: 
        :return: 
        """

        # Find the index of the corresponding row
        index = tables.find_index(self, command)
        self.remove_row(index)

    # -----------------------------------------------------------------

    def remove_entry_and_save(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        self.remove_entry(command)
        self.save()

    # -----------------------------------------------------------------

    def remove_entries(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        indices = tables.find_indices(self, command)
        self.remove_rows(indices)

    # -----------------------------------------------------------------

    def remove_entries_and_save(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        self.remove_entries(command)
        self.save()

    # -----------------------------------------------------------------

    @property
    def commands(self):

        """
        This function ...
        :return:
        """

        return list(self["Command"])

    # -----------------------------------------------------------------

    @property
    def unique_commands(self):

        """
        This function ...
        :return:
        """

        return list(set(self["Command"]))

    # -----------------------------------------------------------------

    @property
    def finished_commands(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def __contains__(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        return command_name in self.commands

    # -----------------------------------------------------------------

    def is_finished_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return not self["End time"].mask[index] # not masked

    # -----------------------------------------------------------------

    def is_finished(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        if command_name not in self.commands: return False
        else:
            index = tables.find_index(self, command_name)
            return self.is_finished_index(index)

    # -----------------------------------------------------------------

    def has_finished_commands(self, command_names):

        """
        This function ...
        :param command_names:
        :return:
        """

        for command in command_names:
            if not self.is_finished(command): return False
        return True

    # -----------------------------------------------------------------

    def has_finished_all(self, *command_names):

        """
        This function ...
        :param command_names:
        :return:
        """

        return self.has_finished_commands(command_names)

    # -----------------------------------------------------------------

    def has_finished_any(self, *command_names):

        """
        This function ...
        :param command_names:
        :return:
        """

        for command in command_names:
            if self.is_finished(command): return True
        return False

    # -----------------------------------------------------------------

    @property
    def finished_maps(self):

        """
        This function ...
        :return:
        """

        return self.has_finished_all("make_old_stellar_maps", "make_young_stellar_maps", "make_ionizing_stellar_maps", "make_dust_map")

    # -----------------------------------------------------------------

    @property
    def finished_maps_commmands(self):

        """
        This function ...
        :return:
        """

        from ..maps.component import maps_commands

        finished = []
        for command in maps_commands:
            if self.is_finished(command): finished.append(command)

        return finished

    # -----------------------------------------------------------------

    @property
    def has_configured_fit(self):

        """
        This function ...
        :return:
        """

        #return "configure_fit" in self
        return self.is_finished("configure_fit")

    # -----------------------------------------------------------------

    @property
    def has_initialized_fit(self):

        """
        This function ...
        :return:
        """

        #if "initialize_fit_sed" in self: return True
        #elif "initialize_fit_galaxy" in self: return True
        #else: return False

        if self.is_finished("initialize_fit_sed"): return True
        elif self.is_finished("initialize_fit_galaxy"): return True
        else: return False

    # -----------------------------------------------------------------

    def mark_end(self):

        """
        This function ...
        :return:
        """

        timestamp = time.timestamp()

        self._resize_string_column("End time", timestamp)

        # Set the value
        self["End time"].mask[-1] = False
        self["End time"][-1] = timestamp

    # -----------------------------------------------------------------

    def mark_end_and_save(self):

        """
        Thi function ...
        :return:
        """

        self.mark_end()
        self.save()

    # -----------------------------------------------------------------

    def register(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        return RegisterScope(self, cls_or_instance)

    # -----------------------------------------------------------------

    def register_start(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        # Get command name
        command_name = cls_or_instance.command_name()

        # Add entry
        self.add_entry(command_name)

    # -----------------------------------------------------------------

    def register_start_and_save(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        self.register_start(cls_or_instance)
        self.save()

    # -----------------------------------------------------------------

    def register_end(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        # Get command name
        command_name = cls_or_instance.command_name()

        # Check last
        if self["Command"][-1] != command_name: raise ValueError("Last command is '" + self["Command"][-1] + "', not '" + command_name + "'")

        # Mark end
        self.mark_end()

    # -----------------------------------------------------------------

    def register_end_and_save(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        self.register_end(cls_or_instance)
        self.save()

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Clean before saving
        self.clean()

        # Call the function of the base class
        super(ModelingHistory, self).saveto(path)

# -----------------------------------------------------------------

class RegisterScope(object):

    """
    This class ...
    """

    def __init__(self, history, cls_or_instance):
        self.history = history
        self.cls_or_instance = cls_or_instance

    def __enter__(self):
        self.history.register_start_and_save(self.cls_or_instance)
        return None

    def __exit__(self, exc_type, exc_value, traceback):
        #print(exc_type, exc_value, traceback)
        if exc_type is not None: return False # not succesful
        else: self.history.register_end_and_save(self.cls_or_instance) # succesful
        return None

# -----------------------------------------------------------------
