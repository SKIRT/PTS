#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_status Show the status page.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.html.status import StatusPageGenerator
from pts.modeling.core.environment import verify_modeling_cwd
from pts.core.tools import browser

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_flag("generate", "first (re)generate the HTML", False)
config = parse_arguments("show_status", definition)

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Generate the HTML
if config.generate:

    # Generate the HTML
    generator = StatusPageGenerator()
    generator.config.path = modeling_path
    generator.run()

# -----------------------------------------------------------------

# Check whether the status page is present
if not fs.is_file(environment.html_status_path): raise ValueError("The status page is not present")
else: browser.open_path(environment.html_status_path)

# -----------------------------------------------------------------
