#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.core.config.analyse_basic import definition as basic_definition
from pts.core.config.analyse_batch import definition as batch_definition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_positional_optional("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "integer", "ID of the simulation")

# Add flags
definition.add_flag("ignore_missing_data", "ignore missing data when analysing the simulations", False)

# Analyse
definition.add_flag("do_basic", "do basic analysis", True)
definition.add_flag("do_batch", "do batch analysis", True)
definition.add_flag("do_scaling", "do scaling analysis", True)
definition.add_flag("do_extra", "do extra analysis", True)

# -----------------------------------------------------------------

# Import settings for basic, batch analysis
definition.import_section("basic", "basic analysis options", basic_definition)
definition.import_section("batch", "batch analysis options", batch_definition)

# Remove certain settings from the sections
definition.sections["basic"].remove_positional_optional("remote")
definition.sections["basic"].remove_positional_optional("id")
definition.sections["batch"].remove_positional_optional("remote")
definition.sections["batch"].remove_positional_optional("id")

# -----------------------------------------------------------------
