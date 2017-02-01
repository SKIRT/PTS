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
import imp
import importlib
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter
from .imports import ImportsChecker
from .test import PTSTest
from ..tools import time

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

scripts = introspection.get_scripts()
tables = introspection.get_arguments_tables()

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
        self.tests = defaultdict(list)

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

        # Loop over the subprojects
        for subproject in self.test_names:

            # Tests path for subproject
            tests_path = fs.join(introspection.pts_subproject_dir(subproject), "tests")

            # Loop over the test names
            for name in self.test_names[subproject]:

                # Determine the test path
                test_path = fs.join(tests_path, name)

                # Determine an output path for the test
                temp_tests_path = fs.create_directory_in(introspection.pts_temp_dir, "tests")
                temp_path = fs.create_directory_in(temp_tests_path, time.unique_name(name))

                # Find file with name test.py
                filepath = fs.join(test_path, "test.py")

                # Load the test module
                test_module = imp.load_source(name, filepath)

                # Get properties of the test module
                description = test_module.description

                # Iterate over these:
                commands = test_module.commands
                settings = test_module.settings
                cwds = test_module.cwds

                setup_function = test_module.setup
                test_function = test_module.test

                # Create Test instance
                test = PTSTest(name, description, setup_function, test_function, temp_path)

                # Loop over the commands
                for command, setting_dict, cwd in zip(commands, settings, cwds):

                    # Find matches
                    matches = introspection.find_matches_scripts(command, scripts)
                    table_matches = introspection.find_matches_tables(command, tables)

                    # Find match
                    if len(matches) + len(table_matches) == 0: raise ValueError("Invalid PTS command: '" + command + "'")
                    elif len(matches) == 1 and len(table_matches) == 0: raise NotImplementedError("Not implemented yet")
                    elif len(table_matches) == 1 and len(matches) == 0:

                        subproject, index = table_matches[0]
                        command_name = tables[subproject]["Command"][index]
                        hidden = False
                        if command_name.startswith("*"):
                            hidden = True
                            command_name = command_name[1:]
                        command_description = tables[subproject]["Description"][index]
                        class_path_relative = tables[subproject]["Path"][index]
                        class_path = "pts." + subproject + "." + class_path_relative
                        module_path, class_name = class_path.rsplit('.', 1)

                        configuration_method_table = tables[subproject]["Configuration method"][index]
                        subproject_path = introspection.pts_subproject_dir(subproject)

                        # Get the class of the configurable of which an instance has to be created
                        module = importlib.import_module(module_path)
                        try: cls = getattr(module, class_name)
                        except AttributeError:
                            raise Exception("The class name for the '" + command_name + "' command is incorrectly specified in the 'commands.dat' file of the '" + subproject + "' subproject")

                        configuration_name = tables[subproject]["Configuration"][index]
                        if configuration_name == "--": configuration_name = command_name
                        configuration_module_path = "pts." + subproject + ".config." + configuration_name

                        # print(configuration_module_path)

                        # try:
                        configuration_module = importlib.import_module(configuration_module_path)
                        # has_configuration = True
                        definition = getattr(configuration_module, "definition")

                        # Parse the configuration
                        setter = DictConfigurationSetter(setting_dict, name, description=None)
                        config = setter.run(definition)

                        # Set working directory (output directory)
                        output_path = fs.absolute_path(fs.join(temp_path, cwd))
                        config.path = output_path

                        # Create the class instance, configure it with the configuration settings
                        inst = cls(config)

                        # Add component to the test
                        test.add_component(inst)

                    # Ambigious command
                    else: raise ValueError("The command '" + command + "' is ambigious")

                # Add the test to the suite
                self.tests[subproject].append(test)

    # -----------------------------------------------------------------

    def run_tests(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running the tests ...")

        # Loop over the subprojects
        for subproject in self.tests:

            # Inform the user
            log.info("Running tests for " + subproject + " subproject ...")

            # Loop over the tests and perform them
            counter = 0
            for test in self.tests[subproject]:

                # Debugging
                log.debug("Performing test " + str(counter+1) + " of " + str(len(self.tests[subproject])))

                # Perform
                test.perform()

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
