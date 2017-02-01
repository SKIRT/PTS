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

# Import standard modules
import importlib
import imp

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter
from .imports import ImportsChecker

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

        # The import statements checker
        self.checker = None

        # The test names
        self.test_names = dict()

        # The tests
        self.tests = dict()

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
        if not self.test_names_for_all_subprojects: self.prompt()

        # 3. Check the import statements
        self.check_imports()

        # 4. Load tests
        self.load_tests()

        # 5. Run tests
        self.run_tests()

        # 6. Show
        if self.config.show: self.show()

        # 7. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    @property
    def test_names_for_all_subprojects(self):

        """
        This function ...
        :return:
        """

        for subproject in self.config.subprojects:

            if subproject not in self.test_names: return False
            elif len(self.test_names[subproject]) == 0: return False

        # No problems encountered
        return True

    # -----------------------------------------------------------------

    def has_test_names(self, subproject):

        """
        This function ...
        :return:
        """

        return subproject in self.test_names and len(self.test_names[subproject]) > 0

    # -----------------------------------------------------------------

    @property
    def subprojects_without_test_names(self):

        """
        This function ...
        :return:
        """

        subprojects = []
        for subproject in self.config.subprojects:
            if self.has_test_names(subproject): continue
            subprojects.append(subproject)
        return subprojects

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PTSTestSuite, self).setup(**kwargs)

        # Tests are specified
        if "tests" in kwargs:
            tests = kwargs.pop("tests")
            assert isinstance(tests, dict)
            self.test_names = tests

        # Check whether the specific tests that have been defined exist
        if self.config.tests is not None:

            # Check whether only one subproject has been specified
            if len(self.config.subprojects) > 1: raise ValueError("Can only specify tests when the number of specified subprojects is 1")

            # Get the test names
            subproject = self.config.subprojects[0]
            subproject_tests = tests_for_subproject(subproject)

            # Check
            for test in self.config.tests:
                if test not in subproject_tests: raise ValueError("Test '" + test + "' for the " + subproject + " subproject does not exist")

            # Set the test names
            self.test_names[subproject] = self.config.tests

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for test names ...")

        # Create definition
        definition = ConfigurationDefinition()

        # Loop over the specified subprojects
        for subproject in self.config.subprojects:
            if self.has_test_names(subproject): continue
            definition.add_required(subproject + "_tests", "string_list", "test to perform from the " + subproject + " subproject", choices=tests_for_subproject(subproject))

        # Get config
        setter = InteractiveConfigurationSetter("subproject_tests", add_logging=False, add_cwd=False)
        config = setter.run(definition, prompt_optional=False)

        # Set the tests
        for subproject in self.subprojects_without_test_names:

            # Get tests
            self.test_names[subproject] = config[subproject + "_tests"]

    # -----------------------------------------------------------------

    def check_imports(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the validity of import statements ...")

        # Create and run the imports checker
        self.checker = ImportsChecker()
        self.checker.config.show = False
        self.checker.config.write = False
        self.checker.run()

    # -----------------------------------------------------------------

    def load_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the tests ...")

        print(self.test_names)

        # Loop over the subprojects
        for subproject in self.test_names:

            # Tests path for subproject
            tests_path = fs.join(introspection.pts_subproject_dir(subproject), "tests")

            # Loop over the test names
            for name in self.test_names[subproject]:

                # Determine the test path
                path = fs.join(tests_path, name)

                # Find file with name test.py
                filepath = fs.join(path, "test.py")
                #pythonic_path = filepath.split("PTS/")[1].replace("/", ".")[:-3]

                # Load the python file
                #test_module = importlib.import_module(pythonic_path)

                # Load the test module
                print(filepath)
                test_module = imp.load_source(name, filepath)

                commands = test_module.commands
                settings = test_module.settings
                setup_function = test_module.setup
                test_function = test_module.test

    # -----------------------------------------------------------------

    def run_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running the tests ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing results ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write report
        self.write_report()

    # -----------------------------------------------------------------

    def write_report(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing report ...")

# -----------------------------------------------------------------
