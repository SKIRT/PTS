#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.vpn Contains the VPN class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import time
import warnings
from subprocess import call, check_output

# -----------------------------------------------------------------

class VPN(object):

    """
    This class ...
    """

    def __init__(self, service):

        """
        This function ...
        :param service:
        :return:
        """

        # Check if the specified VPN service exists
        if service not in get_vpn_services(): raise ValueError("The specified VPN service is not configured on this system")

        # Set the service name
        self.service = service

    # -----------------------------------------------------------------

    def connect(self, user_name=None, password=None, secret=None, delay=None):

        """
        This function ...
        # From: http://apple.stackexchange.com/questions/128297/how-to-create-a-vpn-connection-via-terminal
        :param user_name:
        :param password:
        :param secret:
        :param delay:
        :return:
        """

        # Show a warning if the connection already exists
        if self.connected:
            warnings.warn("The VPN connection is already active")
            return

        command = ["scutil", "--nc", "start", self.service]

        if user_name is not None: command += ["--user", user_name]
        if password is not None: command += ["--password", password]
        if secret is not None: command += ["--secret", secret]

        call(command)

        if delay is not None: time.sleep(delay) # give the user the time to fill in the password (5 seconds)
        # possible solution for this: http://www.proposedsolution.com/solutions/vpn-ipsec-prompting-saved-password/
        # Tried, but didn't work ...

        if not self.connected: raise RuntimeError("Could not connect to the VPN service")

    # -----------------------------------------------------------------

    def disconnect(self):

        """
        This function ...
        :return:
        """

        if not self.connected:
            warnings.warn("The VPN connection is not active")
            return

        command = ["scutil", "--nc", "stop", self.service]
        call(command)

        if self.connected: raise RuntimeError("Could not disconnect from the VPN service")

    # -----------------------------------------------------------------

    @property
    def connected(self):

        """
        This function ...
        :return:
        """

        return self.status == "connected"

    # -----------------------------------------------------------------

    @property
    def status(self):

        """
        This function ...
        :return:
        """

        command = ["scutil", "--nc", "status", self.service]
        output = check_output(command)

        # Return the status
        return output.split("\n")[0].lower()

# -----------------------------------------------------------------

# From: http://apple.stackexchange.com/questions/128297/how-to-create-a-vpn-connection-via-terminal
def get_vpn_services():

    """
    This function ...
    :return:
    """

    vpns_string = check_output(["scutil", "--nc", "list"]) # lists all VPN services
    vpns = re.findall('"(.+)"', vpns_string) # service names are double-quoted

    return vpns

# -----------------------------------------------------------------
