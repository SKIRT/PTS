#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.ensemble Contains the RemotesEnsemble and SKIRTRemotesEnsemble class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools.utils import LazyDictionary
from .remote import Remote
from .host import load_host
from ..simulation.remote import SKIRTRemote
from ..tools import types
from ..tools.utils import memoize_method

# -----------------------------------------------------------------

class RemotesEnsemble(LazyDictionary):

    """
    This class ...
    """

    evaluator = Remote

    # -----------------------------------------------------------------

    def __init__(self, host_ids=None):

        """
        This function ...
        :param host_ids:
        """

        # Call the constructor of the base class
        super(RemotesEnsemble, self).__init__(self.evaluator)

        # Set the host IDS
        if host_ids is not None: self.add_host_ids(host_ids)

    # -----------------------------------------------------------------

    def is_used(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.evaluated[name]

    # -----------------------------------------------------------------

    def add_host_ids(self, host_ids):

        """
        This function ...
        :param host_ids:
        :return:
        """

        # Dictinoary
        if types.is_dictionary(host_ids):

            # Add them under a custom name
            for name in host_ids: self.add_host_id(host_ids[name], name=name)

        # Sequence
        elif types.is_sequence(host_ids) or types.is_tuple(host_ids):

            # Add host IDs with their IDs as key
            for host_id in host_ids: self.add_host_id(host_id)

        # Invalid
        else: raise ValueError("Invalid type for host IDs: " + str(type(host_ids)))

    # -----------------------------------------------------------------

    def add_host_id(self, host_id, name=None):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        # Set name
        if name is None: name = host_id

        # Set the item
        self.set(name, host_id)

        # Return the name
        return name

    # -----------------------------------------------------------------

    def remove_host_id(self, name):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Check name
        if name not in self.names: raise ValueError("No remote under the name '" + name + "' in the ensemble")

        # Remove
        del self[name]

    # -----------------------------------------------------------------

    def add_remote(self, remote, name=None):

        """
        This function ...
        :param remote:
        :param name:
        :return:
        """

        # Check whether remote instance
        if not isinstance(remote, Remote): raise ValueError("First argument must be Remote instance")

        # Get name
        if name is None: name = remote.host_id

        # Check whether not yet present
        if name in self.names: raise ValueError("Already an entry for name '" + name + "'")

        # Set the value
        self.set_value(name, remote)

    # -----------------------------------------------------------------

    def set_remote(self, remote, name=None):

        """
        This function ...
        :param remote:
        :param name:
        :return:
        """

        # Check whether remote instance
        if not isinstance(remote, Remote): raise ValueError("First argument must be Remote instance")

        # Get name
        if name is None: name = remote.host_id

        # Check whether present
        if name not in self.names: raise ValueError("No entry for name '" + name + "'")

        # Check whether not used yet
        if self.is_used(name): raise ValueError("Remote already set for name '" + name + "'")

        # Set
        self.set_value(name, remote)

    # -----------------------------------------------------------------

    def set_or_add_remote(self, remote, name=None):

        """
        This function ...
        :param remote:
        :param name:
        :return:
        """

        # Check whether remote instance
        if not isinstance(remote, Remote): raise ValueError("First argument must be Remote instance")

        # Get name
        if name is None: name = remote.host_id

        # Set
        if name in self.names: self.set_remote(remote, name=name)

        # Add
        else: self.add_remote(remote, name=name)

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.keys()

    # -----------------------------------------------------------------

    @memoize_method
    def get_host_id(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get raw
        value = self.get_raw(name)

        # Return
        if types.is_string_type(value): return value
        else: return value.host_id # assume it is a Remote instance

    # -----------------------------------------------------------------

    @memoize_method
    def get_host(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get raw
        value = self.get_raw(name)

        # Return
        if types.is_string_type(value): return load_host(value)
        else: return value.host

    # -----------------------------------------------------------------

    def is_scheduler(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_host(name).scheduler

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function ...
        :return:
        """

        return [self.get_host_id(name) for name in self.names]

    # -----------------------------------------------------------------

    @property
    def hosts(self):

        """
        This function ...
        :return:
        """

        return [self.get_host(name) for name in self.names]

    # -----------------------------------------------------------------

    @property
    def nremotes(self):

        """
        This function ...
        :return:
        """

        return len(self.names)

    # -----------------------------------------------------------------

    @property
    def has_remotes(self):

        """
        This function ...
        :return:
        """

        return self.nremotes > 0

    # -----------------------------------------------------------------

    @property
    def has_single_remote(self):

        """
        This function ...
        :return:
        """

        return self.nremotes == 1

    # -----------------------------------------------------------------

    @property
    def single_name(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_remote: raise ValueError("Doesn't have a single host ID")
        return self.names[0]

    # -----------------------------------------------------------------

    @property
    def single_host_id(self):

        """
        This function ...
        :return:
        """

        return self.get_host_id(self.single_name)

    # -----------------------------------------------------------------

    @property
    def single_host(self):

        """
        This function ...
        :return:
        """

        return self.get_host(self.single_name)

    # -----------------------------------------------------------------

    @property
    def single_remote(self):

        """
        This function ...
        :return:
        """

        return self[self.single_name]

# -----------------------------------------------------------------

class SKIRTRemotesEnsemble(RemotesEnsemble):

    """
    This class ...
    """

    evaluator = SKIRTRemote

    # -----------------------------------------------------------------

    def __init__(self, host_ids=None):

        """
        The constructor ...
        :param host_ids:
        """

        # Create evaluator for the lazy dicts below
        screens_eval = lambda name: self[name].screen_states() if not self.is_scheduler(name) else None
        jobs_eval = lambda name: self[name].get_jobs_status() if self.is_scheduler(name) else None

        # Initialize dictionaries for the screen states and jobs status
        self.screens = LazyDictionary(screens_eval)
        self.jobs = LazyDictionary(jobs_eval)

        # Call the constructor of the base class
        super(SKIRTRemotesEnsemble, self).__init__(host_ids=host_ids)

    # -----------------------------------------------------------------

    def add_host_id(self, host_id, name=None):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        # Add, get the name
        name = super(SKIRTRemotesEnsemble, self).add_host_id(host_id, name=name)

        # Add to screens and jobs
        self.screens.set(name, name)
        self.jobs.set(name, name)

    # -----------------------------------------------------------------

    def remove_host_id(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Remove
        super(SKIRTRemotesEnsemble, self).remove_host_id(name)

        # Remove from screens and jobs
        del self.screens[name]
        del self.jobs[name]

    # -----------------------------------------------------------------

    def reset_screens(self, names=None):

        """
        This function ...
        :param names: names of the hosts to reset
        :return:
        """

        # Reset all to just the names
        if names is None: names = self.names
        for name in names: self.screens.set(name, name)

    # -----------------------------------------------------------------

    def reset_jobs(self, names=None):

        """
        This function ...
        :param names: names of the hosts to reset
        :return:
        """

        # Reset all to just the names
        if names is None: names = self.names
        for name in names: self.jobs.set(name, name)

# -----------------------------------------------------------------
