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
        description = test_module.description

        # Iterate over these:
        commands = test_module.commands
        input_dicts = test_module.input_dicts
        settings = test_module.settings
        cwds = test_module.cwds

        # Loop over the commands
        for command, input_dict, settings_dict, cwd in zip(commands, input_dicts, settings, cwds):

            # Find matches
            matches = introspection.find_matches_scripts(command, scripts)
            table_matches = introspection.find_matches_tables(command, tables)

# -----------------------------------------------------------------
