#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.test.pts import subprojects_with_tests

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("subprojects", "string_list", "update on a remote system", subprojects_with_tests())
definition.add_positional_optional("tests", "string_list", "test(s) to perform (when one subproject is specified)")

# Add flags
definition.add_flag("keep", "keep the output")
definition.add_flag("show", "show results")
definition.add_flag("write", "write results")
definition.add_flag("open_output", "open the output directory after each test for manual inspection")

# Advanced
definition.add_flag("all", "perform all tests for the given subproject")
definition.add_flag("default", "use all default options for the specific tests")

definition.add_flag("remove_previous", "remove the output of previous runs of the same tests that are executed now")

# -----------------------------------------------------------------
