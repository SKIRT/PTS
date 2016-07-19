#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.hostconfigurer Contains the HostConfigurer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools.logging import log
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, write_config

# -----------------------------------------------------------------

# The path to the hosts configuration directory
hosts_config_path = fs.join(introspection.pts_config_dir("core"), "hosts")

# Determine the path to the template host configuration file
raw_template_path = fs.join(hosts_config_path, "template.cfg")

# Determine the path to the user hosts directory
hosts_directory = fs.join(introspection.pts_user_dir, "hosts")

# -----------------------------------------------------------------

class HostConfigurer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(HostConfigurer, self).__init__(config)

        # The configuration setter
        self.setter = InteractiveConfigurationSetter("host")

        # The configuration definition
        self.definition = None

        # The host configuration
        self.host_config = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the host template configuration
        self.load_template()

        # 3. Prompt for the settings
        self.prompt()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(HostConfigurer, self).setup()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Determine template path (raw template for new host or one of the preconfigured templates)
        if self.config.preconfigured is None: template_path = raw_template_path
        else: template_path = fs.join(hosts_config_path, self.config.preconfigured + ".cfg")

        # Load the configuration template
        self.definition = ConfigurationDefinition.from_file(template_path)

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Create the configuration
        self.config = self.setter.run(self.definition)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the host configuration
        self.write_config()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the host file
        path = fs.join(hosts_directory, self.config.name + ".cfg")

        # Write the host configuration
        write_config(self.host_config, path)

# -----------------------------------------------------------------
