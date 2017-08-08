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
from ..basics.log import log
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from ..tools import formatting as fmt
from ..basics.map import Map

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HostConfigurer, self).__init__(*args, **kwargs)

        # The configuration setter
        self.setter = InteractiveConfigurationSetter("host", add_cwd=False, add_logging=False)

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

        # Show
        self.show()

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
        self.definition = ConfigurationDefinition.from_file(template_path, write_config=False)

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Create the configuration
        self.host_config = self.setter.run(self.definition, prompt_optional=True)

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
        self.host_config.saveto(path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the settings ...")

        kwargs = dict()
        kwargs["round"] = True
        kwargs["scientific"] = True
        kwargs["ndigits"] = 3

        print("")

        # Print in 2 columns
        with fmt.print_in_columns(3, delimiter="   ", indent=" ", tostr_kwargs=kwargs) as print_row:

            # Print row for each variable
            for name in self.host_config.keys():

                if name.startswith("_"): continue

                # Mapping
                if isinstance(self.host_config[name], Map):

                    prefix = name
                    for nname in self.host_config[name].keys():

                        #print(name + "/" + nname)

                        value = self.host_config[name][nname]
                        the_name = prefix + "/" + nname
                        print_row(the_name, ":", value)

                # Other setting
                else: print_row(name, ":", self.host_config[name])

        print("")

# -----------------------------------------------------------------
