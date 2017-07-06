#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_attenuation Check the attenuation correction of the prepared images.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.preparation.preparer import load_statistics
from pts.core.filter.filter import parse_filter
from pts.magic.services.attenuation import GalacticAttenuation
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Create the configuration
config = parse_arguments("check_attenuation", definition)

# -----------------------------------------------------------------

# Modeling path
modeling_path = fs.cwd()

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------
