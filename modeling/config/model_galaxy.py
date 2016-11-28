#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids
from pts.modeling.modeler import modeling_methods, modeling_methods_descriptions

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add required arguments
definition.add_required("galaxy_name", "string", "name of the galaxy for which to perform the radiative transfer modeling")
definition.add_required("method", "string", "method to use for the modeling", choices=modeling_methods, choice_descriptions=modeling_methods_descriptions)
definition.add_required("host_id", "string", "remote host to use for heavy computations", choices=find_host_ids())

# Add optional arguments
definition.add_optional("fitting_host_ids", "string_list", "remote hosts to use for performing simulations as part of the fitting", choices=find_host_ids(), default=find_host_ids())

# -----------------------------------------------------------------
