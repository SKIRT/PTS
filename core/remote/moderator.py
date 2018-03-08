#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.ensemble Contains the RemotesEnsemble class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .remote import is_available
from ..basics.log import log
from .host import load_host
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

class PlatformModerator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Host ids of the available hosts
        self.available_host_ids = set()

        # Attributes
        self.local = []
        self.single = dict()
        self.ensemble = dict()

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        names_list = []
        names_list += self.local
        names_list += self.single.keys()
        names_list += self.ensemble.keys()
        return names_list

    # -----------------------------------------------------------------

    def has_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name in self.local: return True
        elif name in self.single: return True
        elif name in self.ensemble: return True
        else: return False

    # -----------------------------------------------------------------

    def add_local(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_name(name): raise ValueError("Already set the allowed platform(s) for procedure '" + name + "'")
        self.local.append(name)

    # -----------------------------------------------------------------

    def add_single(self, name, host_ids):

        """
        This function ...
        :return:
        """

        if host_ids is None: raise ValueError("Cannot pass None as host IDs, use add_local")

        if self.has_name(name): raise ValueError("Already set the allowed platform(s) for procedure '" + name + "'")
        self.single[name] = host_ids

    # -----------------------------------------------------------------

    def add_ensemble(self, name, host_ids):

        """
        This function ...
        :return:
        """

        if host_ids is None: raise ValueError("Cannot pass None as host IDs, use add_local")

        if self.has_name(name): raise ValueError("Already set the allowed platform(s) for procedure '" + name + "'")
        self.ensemble[name] = host_ids

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding available hosts ...")

        # Loop over the names
        encountered_host_ids = set()
        for name in self.names:

            # Find available hosts from host_ids list
            for host_id in self.allowed_host_ids(name):

                # Check if already encountered
                if host_id in encountered_host_ids: continue
                encountered_host_ids.add(host_id)

                # Check if available
                if is_available(host_id):

                    log.debug("Host '" + host_id + "' is available")
                    self.available_host_ids.add(host_id)

                # Not available
                else: log.debug("Host '" + host_id + "' is not available")

        # No available host in the list of preferred host ids
        if len(self.available_host_ids) == 0 and self.has_allowed_host_ids_all():
            raise RuntimeError("None of the preferred hosts are available at this moment")

    # -----------------------------------------------------------------

    def allowed_host_ids(self, name):

        """
        This function ...
        :param name: name of the procedure
        :return:
        """

        if name in self.local: return []
        elif name in self.single: return self.single[name]
        elif name in self.ensemble: return self.ensemble[name]
        else: raise ValueError("No allowed platforms defined for procedure '" + name + "'")

    # -----------------------------------------------------------------

    def has_allowed_host_ids(self, name):

        """
        This function ...
        :return:
        """

        return len(self.allowed_host_ids(name)) > 0

    # -----------------------------------------------------------------

    def has_allowed_host_ids_all(self):

        """
        This function ...
        :return:
        """

        for name in self.names:
            if self.has_allowed_host_ids(name): return True
        return False

    # -----------------------------------------------------------------

    def host_id_for_single(self, name):

        """
        This function returns the ID of the host used for the procedure 'name' where only one remote should be used
        :return:
        """

        # Loop over the preferred hosts
        for host_id in self.allowed_host_ids(name):

            # Get the first available host
            if host_id in self.available_host_ids: return host_id

        # No host avilable
        return None

    # -----------------------------------------------------------------

    def host_ids_for_ensemble(self, name, none_if_none=False):

        """
        This function ...
        :return:
        """

        host_ids = []

        # Loop over the specified hosts
        for host_id in self.allowed_host_ids(name):

            # Add each available host to the list
            if host_id in self.available_host_ids: host_ids.append(host_id)

        # Return the list of host IDs
        if len(host_ids) == 0 and none_if_none: return None
        else: return host_ids

    # -----------------------------------------------------------------

    @property
    def all_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = set()

        # Loop over the names in 'single'
        for name in self.single:
            host_id = self.host_id_for_single(name)
            if host_id is not None: host_ids.add(host_id)

        # Loop over the names in 'ensemble'
        for name in self.ensemble: host_ids |= set(self.host_ids_for_ensemble(name))

        # Return the list of host IDs
        return list(host_ids)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_hosts(self):

        """
        This function ...
        :return:
        """

        return [load_host(host_id) for host_id in self.all_host_ids]

    # -----------------------------------------------------------------

    def single_is_local(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.host_id_for_single(name) is None

    # -----------------------------------------------------------------

    def ensemble_is_local(self, name):

        """
        This function ...
        :return:
        """

        return len(self.host_ids_for_ensemble(name)) == 0

    # -----------------------------------------------------------------

    @property
    def all_is_local(self):

        """
        This function ...
        :return:
        """

        return len(self.all_host_ids) == 0

    # -----------------------------------------------------------------

    @property
    def any_remote(self):

        """
        This function ...
        :return:
        """

        return len(self.all_host_ids) > 0

    # -----------------------------------------------------------------

    def host_ids_non_schedulers_ensemble(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        host_ids = []
        for host_id in self.host_ids_for_ensemble(name):

            host = load_host(host_id)

            # Skip schedulers
            if host.scheduler: continue
            else: host_ids.append(host_id)

        return host_ids

    # -----------------------------------------------------------------

    def host_ids_schedulers_ensemble(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        host_ids = []
        for host_id in self.host_ids_for_ensemble(name):

            host = load_host(host_id)

            # Skip non-schedulers
            if not host.scheduler: continue
            else: host_ids.append(host_id)

        return host_ids

    # -----------------------------------------------------------------

    def clear_all_hosts(self):

        """
        This function ...
        :return:
        """

        from .remote import Remote

        # Loop over the hosts
        for host_id in self.all_host_ids:

            # Inform the user
            log.info("Clearing remote '" + host_id + "' ...")

            # Setup the remote
            remote = Remote()
            if not remote.setup(host_id):
                log.warning("Could not connect to remote host '" + host_id + "'")
                continue

            # Clear temporary data and close sessions
            remote.clear_temp_and_sessions()

# -----------------------------------------------------------------
