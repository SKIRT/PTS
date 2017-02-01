#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.pts Contains the PTSTestSuite class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter

# -----------------------------------------------------------------

def subprojects_with_tests():

    """
    This function ...
    :return:
    """

    subprojects = []
    for subproject in introspection.subprojects:
        tests_path = fs.join(introspection.pts_subproject_dir(subproject), "tests")
        if fs.is_directory(tests_path) and not fs.is_empty(tests_path): subprojects.append(subproject)
    return subprojects

# -----------------------------------------------------------------

def tests_for_subproject(subproject):

    """
    This function ...
    :param subproject:
    :return:
    """

    # Get test path for subproject
    subproject_path = introspection.pts_subproject_dir(subproject)
    tests_path = fs.join(subproject_path, "tests")

    # If directory doesn't exist, return empty list for the tests
    if not fs.is_directory(tests_path): return []
    else: return fs.directories_in_path(tests_path, returns="name")

# -----------------------------------------------------------------

class PTSTestSuite(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):


        """
        This class ...
        :param config:
        """

        # Call the constructor of the base class
        super(PTSTestSuite, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Prompt for which test has to be executed
        self.prompt()

        # Test the import statements
        self.test_imports()

        # Load tests
        self.load_tests()

        # Run tests
        self.run_tests()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PTSTestSuite, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition()

        # Loop over the specified subprojects
        for subproject in self.config.subprojects:
            definition.add_required(subproject + "_tests", "string_list", "test to perform from the " + subproject + " subproject", choices=tests_for_subproject(subproject))

        # Get config
        setter = InteractiveConfigurationSetter("subproject_tests", add_logging=False, add_cwd=False)
        config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def test_imports(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the validity of import statements ...")

    # -----------------------------------------------------------------

    def test_external_imports(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing external import statements ...")

    # -----------------------------------------------------------------

    def test_internal_imports(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing internal import statements ...")

    # -----------------------------------------------------------------

    def load_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the tests ...")

    # -----------------------------------------------------------------

    def run_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running the tests ...")

# -----------------------------------------------------------------
