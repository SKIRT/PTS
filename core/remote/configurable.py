#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.configurable Contains the RemotesConfigurable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from pts.core.remote.remote import Remote
from pts.core.basics.log import log
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

class RemotesConfigurable(Configurable):
    
    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RemotesConfigurable, self).__init__(*args, **kwargs)

        # The remotes
        self.remotes = []

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function ...
        :return:
        """

        ids = []

        # Loop over the remotes
        for remote in self.remotes: ids.append(remote.host_id)

        # Return the IDs
        return ids

    # -----------------------------------------------------------------

    @property
    def has_remotes(self):

        """
        This function ...
        :return:
        """

        return len(self.remotes) > 0

    # -----------------------------------------------------------------

    @property
    def nremotes(self):

        """
        This function ...
        :return:
        """

        return len(self.remotes)

    # -----------------------------------------------------------------

    @property
    def single_remote(self):

        """
        This function ...
        :return:
        """

        if self.nremotes != 1: raise RuntimeError("Do not contain exactly one remote")
        return self.remotes[0]

    # -----------------------------------------------------------------

    @property
    def single_host_id(self):

        """
        THis function ...
        :return:
        """

        return self.single_remote.host_id

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RemotesConfigurable, self).setup(**kwargs)

        # If remotes are passed
        if "remotes" in kwargs:

            # Get the remotes
            self.remotes = kwargs.pop("remotes")

            # Check if they are setup
            for remote in self.remotes:
                if not remote.connected: raise RuntimeError("Remotes must be connected")

        # Remotes are not passed
        else:

            # Gather host IDs
            #if self.config.not_remotes is not None: host_ids = [host_id for host_id in self.config.host_ids if host_id not in self.config.not_remotes]
            #else: host_ids = self.config.host_ids

            # Loop over the different hosts
            #for host_id in host_ids:
            for host in self.config.hosts:

                # Create the remote object
                remote = Remote(log_conda=kwargs.pop("log_conda", False))

                # Login to the remote
                #if not remote.setup(host_id, one_attempt=self.config.one_attempt, cluster_name=self.config.clustername):
                if not remote.setup(host, one_attempt=self.config.one_attempt):
                    log.warning("Remote host '" + host.id + "' is down: skipping")
                    continue

                # Add the remote to the list
                self.remotes.append(remote)

# -----------------------------------------------------------------
