#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.list_parameters

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import formatting as fmt
from pts.modeling.config.parameters import possible_parameter_types, possible_parameter_types_descriptions, default_units
from pts.core.units.stringify import represent_unit as ru
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

definition = ConfigurationDefinition()

definition.add_positional_optional("model_name", "string", "model name")
config = parse_arguments(definition)

# -----------------------------------------------------------------
