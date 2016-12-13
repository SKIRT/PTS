#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.modeler import modeling_methods

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add required settings
definition.add_required("galaxy_name", "string", "name of the galaxy for which to initiate the radiative transfer modeling")
definition.add_required("method", "string", "method to use for the modeling", choices=modeling_methods)
definition.add_required("host_ids", "string_list", "remote hosts to use for heavy computations (in order of preference)", choices=find_host_ids(schedulers=False))

# Add optional arguments
definition.add_optional("fitting_host_ids", "string_list", "remote hosts to use for performing simulations as part of the fitting", choices=find_host_ids(), default=find_host_ids())

# -----------------------------------------------------------------
