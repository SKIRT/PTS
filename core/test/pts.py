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
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from .imports import ImportsChecker
from .test import PTSTest
from ..tools import time
from ..tools import stringify
from ..remote.utils import DetachedCalculation
from .table import load_tests_table
from ..basics.configuration import create_configuration_flexible
from ..tools import formatting as fmt

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

def tests_and_descriptions_for_subproject(subproject):

    """
    This function ...
    :param subproject:
    :return:
    """

    # Get test path for subproject
    subproject_path = introspection.pts_subproject_dir(subproject)
    tests_path = fs.join(subproject_path, "tests")

    # Return empty dictionary if there are no tests
    if not fs.is_directory(tests_path): return dict()

    tests = dict()

    # Loop over the tests
    for path, name in fs.directories_in_path(tests_path, returns=["path", "name"]):

        # Determine path of the test.py file
        filepath = fs.join(path, "test.py")

        # Check if the file is present
        if not fs.is_file(filepath):
            log.warning("The test definition for '" + name + "' is not complete")
            continue

        # Load the test module
        test_module = imp.load_source(name, filepath)

        # Get the description
        try: description = test_module.description
        except AttributeError:
            log.warning("Description for test '" + name + "' is not given")
            description = "no description"

        # Add to the tests dictionary
        tests[name] = description

    # Return the tests dictionary
    return tests

# -----------------------------------------------------------------

def path_for_test(subproject, name):

    """
    This function ...
    :param subproject:
    :param name:
    :return:
    """

    # Determine the path
    path = fs.join(introspection.pts_subproject_dir(subproject), "tests", name)

    # Check if exists
    if not fs.is_directory(path): raise ValueError("'" + name + "' is not a test in the '" + subproject + "' subproject")

    # Return the path
    return path

# -----------------------------------------------------------------

tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

class PTSTestSuite(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):


        """
        This class ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PTSTestSuite, self).__init__(*args, **kwargs)

        # The import statements checker
        self.checker = None

        # The test names
        self.test_names = dict()

        # The tests
        self.tests = defaultdict(list)

        # The tests table
        self.table = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Prompt for which test has to be executed
        if not self.test_names_for_all_subprojects and not self.config.only_checks: self.prompt()

        # 3. Check the import statements
        if self.config.check_imports and not self.config.only_tests: self.check_imports()

        # Check the commands
        if self.config.check_commands and not self.config.only_tests: self.check_commands()

        # Check configurations
        if self.config.check_commands and not self.config.only_tests: self.check_configurations()

        # Check package definitions
        if self.config.check_packages and not self.config.only_tests: self.check_packages()

        # 4. Load tests
        if not self.config.only_checks: self.load_tests()

        # 5. Run tests
        if not self.config.only_checks: self.run_tests()

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

        # Load the tests table
        self.table = load_tests_table()

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

        # If the 'all' flag has been enabled
        if self.config.all:
            for subproject in self.config.subprojects:
                self.test_names[subproject] = tests_for_subproject(subproject)

        # If the 'open_output' flag is enabled, also enable the 'keep' flag
        if self.config.open_output: self.config.keep = True

        # Check
        if self.config.settings is not None and len(self.test_names) > 1: raise ValueError("Cannot specifiy test settings when more than one test is getting launched")

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for test names ...")

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Loop over the specified subprojects
        multiple_subprojects = len(self.config.subprojects) > 1
        for subproject in self.config.subprojects:
            if self.has_test_names(subproject): continue
            if multiple_subprojects: definition.add_optional(subproject + "_tests", "string_list", "test to perform from the " + subproject + " subproject", choices=tests_and_descriptions_for_subproject(subproject))
            else: definition.add_required(subproject + "_tests", "string_list", "test to perform from the " + subproject + " subproject", choices=tests_and_descriptions_for_subproject(subproject))

        # Get config
        setter = InteractiveConfigurationSetter("subproject_tests", add_logging=False, add_cwd=False)
        config = setter.run(definition, prompt_optional=multiple_subprojects)

        # Set the tests
        for subproject in self.subprojects_without_test_names:

            # No test must be done for this subproject
            if config[subproject + "_tests"] is None: continue

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

    def check_commands(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Check the command tables ...")

        # Get the argument tables
        tables = introspection.get_arguments_tables()

        # Loop over the subprojects
        for subproject in tables:

            table = tables[subproject]

            # Loop over the commands
            for index in range(len(table["Command"])):

                command_name = table["Command"][index]
                description = table["Description"][index]
                class_path_relative = table["Path"][index]
                class_path = "pts." + subproject + "." + class_path_relative
                module_path, class_name = class_path.rsplit('.', 1)

                #configuration_method = table["Configuration method"][index]

                # Determine the configuration module path
                configuration_name = table["Configuration"][index]
                if configuration_name == "--": configuration_name = command_name
                configuration_module_path = "pts." + subproject + ".config." + configuration_name

                # Find definition
                #try:
                #    configuration_module = importlib.import_module(configuration_module_path)
                #    definition = getattr(configuration_module, "definition")
                #except ImportError:
                #    log.warning("No configuration definition found for the " + class_name + " class")
                #    definition = ConfigurationDefinition()  # Create new configuration definition

                # Find configuration module
                configuration_module_file_path = fs.join(introspection.pts_subproject_dir(subproject), "config", configuration_name + ".py")
                if not fs.is_file(configuration_module_file_path): log.warning("The configuration module cannot be found for the '" + class_name + "' class")

                # Find the class
                try: cls = introspection.get_class(module_path, class_name)
                except (ValueError, ImportError): log.warning("The class '" + class_name + "' could not be found in module '" + module_path + "'")

    # -----------------------------------------------------------------

    def check_configurations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the configuration modules with their corresponding classes ...")

    # -----------------------------------------------------------------

    def check_packages(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the package definitions ...")

        # Check these kinds of definitions at the top of each module:

        # ## \package pts.magic.misc.cortese Contains the GalametzTIRCalibration class.

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

            # Inform the user
            log.info("Loading tests for '" + subproject + "' subproject ...")

            # Loop over the test names
            for name in self.test_names[subproject]:

                # Determine the test path
                test_path = fs.join(tests_path, name)

                # Determine an output path for the test
                unique_name = time.unique_name(name)
                temp_path = fs.create_directory_in(introspection.pts_tests_dir, unique_name)

                # Debugging
                log.debug("Creating temporary directory '" + temp_path + "' for the test '" + name + "' ...")

                # Find file with name test.py
                filepath = fs.join(test_path, "test.py")

                # Find file with the name config.py
                config_path = fs.join(test_path, "config.py")

                # Load the configuration definition
                config_module = imp.load_source(name + "_config", config_path)
                definition = config_module.definition

                # Show test settings
                if self.config.settings is not None:
                    # Debugging
                    log.debug("Setting options for the '" + name + "' test from the following dictionary:")
                    if log.is_debug: fmt.print_dictionary(self.config.settings)

                # Create the test configuration
                config = create_configuration_flexible(name, definition, self.config.settings, self.config.default)

                ### REMOVE PREVIOUS

                # If remove_previous is enabled, remove previous output directories of the same test
                if self.config.remove_previous:

                    # Determine which directories to avoid
                    if "reference_test" in config: exact_not_name = [unique_name, config.reference_test]
                    else: exact_not_name = unique_name

                    # Find
                    previous_paths = fs.directories_in_path(introspection.pts_tests_dir, startswith=name, exact_not_name=exact_not_name)

                    # Debugging
                    if len(previous_paths) > 0: log.debug("Removing output directories of previous test runs: " + stringify.stringify(previous_paths)[1] + " ...")

                    # Remove
                    if "reference_path" in config and config.reference_path is not None:
                        log.debug("Except for when it contains the reference path [" + config.reference_path + "]")
                        fs.remove_directories_but_keep(previous_paths, config.reference_path)
                    else: fs.remove_directories(previous_paths)

                ###

                # Load the test module
                test_module = imp.load_source(name, filepath)

                # Get properties of the test module
                try: description = test_module.description
                except AttributeError: raise RuntimeError("Description not specified for the '" + name + "' test")

                # Get the test function
                #test_function = test_module.test
                test_function = None

                # Get the implementation class
                implementation_cls = introspection.classes_in_module(test_module)[0]

                # Create class instance with the configuration
                implementation = implementation_cls(config)

                # Create Test instance
                test = PTSTest(name, description, implementation, test_function, temp_path, self.config.keep, self.config.open_output)

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
                log.debug("Performing test " + str(counter+1) + " of " + str(len(self.tests[subproject])) + " ...")

                # Set log path
                log_path = fs.join(test.output_path, "log.txt")

                # Write log output to test log file
                with log.log_to_file(log_path, filter_level="DEBUG"):

                    # Start
                    log.start("Starting test '" + test.name + "' ...")

                    # Run the test
                    try: test.run()
                    except DetachedCalculation as detached:

                        # Give warning
                        log.warning("The test '" + test.name + "' of the '" + subproject + "' subproject is being detached: progress and retrieval information are being saved into the tests table ...")

                        # Save the test
                        test.save()

                        # Add an entry to the table
                        self.table.add_test()

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
