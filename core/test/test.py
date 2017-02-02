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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class PTSTest(object):

    """
    This class ...
    """

    def __init__(self, name, description, setup_function, test_function, output_path):

        """
        This function ...
        :param name:
        :param description:
        :param setup_function:
        :param test_function:
        """

        # Properties of this test
        self.name = name
        self.description = description
        self.setup_function = setup_function
        self.test_function = test_function
        self.output_path = output_path

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

        # Perform
        self.perform()

        # Check
        self.check()

        # Show
        self.show()

        # Write
        self.write()

        # Clear
        self.clear()

    # -----------------------------------------------------------------

    def add_component(self, name, component, input_dict=None):

        """
        This function ...
        :param name:
        :param component:
        :param input_dict:
        :return:
        """

        self.components[name] = component
        if input_dict is not None: self.input_dicts[name] = input_dict

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Execute setup function
        self.setup_function()

    # -----------------------------------------------------------------

    def info(self):

        """
        This function ...
        :return:
        """

        print("PTS TEST")
        print(self.name)
        print(self.description)

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
            log.debug("Executing component '" + name + "' ...")

            # Get input
            input_dict = self.input_dicts[name] if name in self.input_dicts else dict()

            # Run with input
            self.components[name].run(**input_dict)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

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
        log.info("")

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        fs.remove_directory(self.output_path)

# -----------------------------------------------------------------
