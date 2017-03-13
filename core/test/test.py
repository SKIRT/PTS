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
from ..tools import time

# -----------------------------------------------------------------

class PTSTest(object):

    """
    This class ...
    """

    #def __init__(self, name, description, setup_function, test_function, output_path, keep=False, open_output=False):
    def __init__(self, name, description, implementation, test_function, output_path, keep=False, open_output=False):

        """
        This function ...
        :param name:
        :param description:
        :param test_function:
        :param keep:
        :param open_output:
        """

        # Properties of this test
        self.name = name
        self.description = description
        #self.setup_function = setup_function
        self.test_function = test_function
        self.output_path = output_path
        self.keep = keep
        self.open_output = open_output

        # The test implementation
        self.implementation = implementation

        # The runnable components
        #self.components = OrderedDict()
        #self.input_dicts = dict()

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

        # 7. Open the output directory for inspection
        if self.open_output: self.open()

        # 8. Clear
        if not self.keep: self.clear()

    # -----------------------------------------------------------------

    #def add_component(self, name, cls, configuration_module_path, settings_dict, output_path, input_dict, description, finish=None):

        #"""
        #This function ...
        #:param name:
        #:param cls:
        #:param configuration_module_path:
        #:param settings_dict:
        #:param output_path:
        #:param input_dict:
        #:param description:
        #:param finish:
        #:return:
        #"""

        # Add the component to the dictionary
        #self.components[name] = Map(cls=cls, conf_path=configuration_module_path, settings=settings_dict, output_path=output_path, input_dict=input_dict, description=description, finish=finish)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the test ...")

        # Execute setup function
        #self.setup_function(self.output_path)

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
        fmt.print_centered_around_border("****", "|", prefix=prefix, length=width)
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

        # Run the test implemenation
        self.implementation.run()

    # -----------------------------------------------------------------

    def perform_old(self):

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

            # Get properties of the component
            configuration_module_path = self.components[name].conf_path
            output_path = self.components[name].output_path
            settings_dict = self.components[name].settings
            cls = self.components[name].cls
            input_dict = self.components[name].input_dict
            description = self.components[name].description
            finish_function = self.components[name].finish

            # Change working directory
            fs.change_cwd(output_path)

            # Get the configuration definition
            configuration_module = importlib.import_module(configuration_module_path)
            definition = getattr(configuration_module, "definition")

            # Parse the configuration
            setter = DictConfigurationSetter(settings_dict, name, description=description)
            config = setter.run(definition)

            # Set working directory (output directory)
            config.path = output_path

            # Create the class instance, configure it with the configuration settings
            inst = cls(config)

            # Inform the user
            log.info("Executing component '" + name + "': " + description + " ...")

            # Run with input
            inst.run(**input_dict)

            # If finish function is defined, call it and pass the instance for which it finishes up
            if finish_function is not None: finish_function(inst)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the test output ...")

        # Execute test function
        self.test_function(self.output_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def open(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Opening the test output directory ...")

        # Open the output directory
        fs.open_directory(self.output_path, wait=True)

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
