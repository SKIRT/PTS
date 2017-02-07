#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.test Contains the PTSTest class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem as fs
from ..basics.map import Map
from ..basics.configuration import DictConfigurationSetter
from ..tools import formatting as fmt
from ..tools import strings

# -----------------------------------------------------------------

class PTSTest(object):

    """
    This class ...
    """

    def __init__(self, name, description, setup_function, test_function, output_path, keep=False):

        """
        This function ...
        :param name:
        :param description:
        :param setup_function:
        :param test_function:
        :param keep:
        """

        # Properties of this test
        self.name = name
        self.description = description
        self.setup_function = setup_function
        self.test_function = test_function
        self.output_path = output_path
        self.keep = keep

        # The runnable components
        self.components = OrderedDict()
        self.input_dicts = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2. Show info
        self.info()

        # 3. Perform
        self.perform()

        # 4. Check
        self.check()

        # 5. Show
        self.show()

        # 6. Write
        self.write()

        # 7. Clear
        if not self.keep: self.clear()

    # -----------------------------------------------------------------

    def add_component(self, name, cls, configuration_module_path, settings_dict, output_path, input_dict):

        """
        This function ...
        :param name:
        :param cls:
        :param configuration_module_path:
        :param settings_dict:
        :param output_path:
        :param input_dict:
        :return:
        """

        # Add the component to the dictionary
        self.components[name] = Map(cls=cls, conf_path=configuration_module_path, settings=settings_dict, output_path=output_path, input_dict=input_dict)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the test ...")

        # Execute setup function
        self.setup_function(self.output_path)

    # -----------------------------------------------------------------

    def info(self):

        """
        This function ...
        :return:
        """

        prefix = "    "
        width = 50

        # Print the info
        fmt.print_empty()
        fmt.print_filled("-", prefix=prefix, length=width)
        fmt.print_border("|", prefix=prefix, length=width)
        fmt.print_centered_around_border("PTS test '" + self.name + "'", "|", prefix=prefix, length=width)
        fmt.print_border("|", prefix=prefix, length=width)
        for line in strings.split_in_lines(self.description, length=40, as_list=True): fmt.print_centered_around_border(line, "|", prefix=prefix, length=width)
        fmt.print_border("|", prefix=prefix, length=width)
        fmt.print_filled("-", prefix=prefix, length=width)
        fmt.print_empty()

    # -----------------------------------------------------------------

    def perform(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the test ...")

        # Loop over the components, invoke their run function
        for name in self.components:

            # Debugging
            log.debug("Setting configuration for component '" + name + "' ...")

            configuration_module_path = self.components[name].conf_path
            output_path = self.components[name].output_path
            settings_dict = self.components[name].settings
            cls = self.components[name].cls
            input_dict = self.components[name].input_dict

            # Change working directory
            fs.change_cwd(output_path)

            configuration_module = importlib.import_module(configuration_module_path)
            definition = getattr(configuration_module, "definition")

            # Parse the configuration
            setter = DictConfigurationSetter(settings_dict, name, description=None)
            config = setter.run(definition)

            # Set working directory (output directory)
            config.path = output_path

            # Create the class instance, configure it with the configuration settings
            inst = cls(config)

            # Debugging
            log.debug("Executing component '" + name + "' ...")

            # Run with input
            inst.run(**input_dict)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the test output ...")

        # Execute test function
        self.test_function()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing the output of the " + self.name + " test ...")

        # Remove the entire output directory
        fs.remove_directory(self.output_path)

# -----------------------------------------------------------------
