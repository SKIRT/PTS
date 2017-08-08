#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.test_coverage Determine the test coverage.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import imp

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.basics.configuration import find_command
from pts.core.tools import formatting as fmt
from pts.core.test.pts import tests_for_subproject, path_for_test
from pts.core.tools import filesystem as fs
from pts.core.tools import sequences
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Get arguments tables
tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

# Get all configurable classes
classes = list(introspection.all_concrete_configurable_classes())

# Initialize container
in_tests = [False] * len(classes)

# -----------------------------------------------------------------

log = setup_log(level="INFO")

# -----------------------------------------------------------------

# Loop over the subprojects
for subproject in introspection.subprojects:

    # Get the test names for this subproject
    for name in tests_for_subproject(subproject):

        # Determine the test path
        test_path = path_for_test(subproject, name)

        # Find file with name test.py
        filepath = fs.join(test_path, "test.py")

        # Load the test module
        test_module = imp.load_source(name, filepath)

        # Get properties of the test module
        try: description = test_module.description
        except AttributeError: log.warning("Description not specified for test '" + name + "'")

        # Iterate over these:
        commands = test_module.commands

        # Loop over the commands
        for command in commands:

            the_command = command.command

            # Find match
            match = introspection.resolve_command_tables(command, tables)

            # Get the class name
            cls = introspection.get_class(match.module_path, match.class_name)

            # Find index of class in list of classes
            index = sequences.find_exact_index(classes, cls)

            # Set to True
            in_tests[index] = True

# -----------------------------------------------------------------

# Report
for index in range(len(classes)):

    class_name = classes[index].__name__
    if in_tests[index]: print(fmt.green + class_name + ": in tests" + fmt.reset)
    else: print(fmt.red + class_name + ": untested" + fmt.reset)

# -----------------------------------------------------------------
