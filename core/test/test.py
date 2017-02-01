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

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

class PTSTest(object):

    """
    This class ...
    """

    def __init__(self, name, description, setup_function, test_function):

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

        # The runnable components
        self.components = []

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        """

        # Setup
        self.setup()

        # Show info
        self.info()

        # Perform
        self.perform()

        # Check
        self.check()

        # Show
        self.show()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def add_component(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        self.components.append(component)

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
        for component in self.components:

            # Debugging
            log.debug("Executing component '" + component.name + "' ...")

            # Run
            component.run()

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
