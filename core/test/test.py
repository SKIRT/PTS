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
from ..basics.log import log
from ..tools import filesystem as fs
from ..tools import formatting as fmt
from ..tools import strings

# -----------------------------------------------------------------

class PTSTest(object):

    """
    This class ...
    """

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
        self.test_function = test_function
        self.output_path = output_path
        self.keep = keep
        self.open_output = open_output

        # The test implementation
        self.implementation = implementation

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
        self.implementation.run(path=self.output_path)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the test output ...")

        # Execute test function
        if self.test_function is not None: self.test_function(self.output_path)

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

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

# -----------------------------------------------------------------
