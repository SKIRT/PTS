#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import filesystem as fs
from pts.modeling.component.component import load_modeling_configuration, load_fitting_configuration, load_modeling_history
from pts.modeling.component.galaxy import load_preparation_statistics
from pts.magic.core.frame import Frame

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("model_info", "Show information about the radiative transfer model")
config = setter.run(definition)

# -----------------------------------------------------------------



# -----------------------------------------------------------------
