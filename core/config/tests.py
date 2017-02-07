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

# -----------------------------------------------------------------
