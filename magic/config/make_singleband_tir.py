#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.maps.tir.single import possible_filters

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# The filters
definition.add_positional_optional("filters", "filter_list", "filters to make TIR maps from", default=possible_filters, choices=possible_filters, convert_default=True)

# -----------------------------------------------------------------
