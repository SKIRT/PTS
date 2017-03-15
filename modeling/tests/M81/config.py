#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Optional settings
definition.add_optional("nwavelengths", "positive_integer", "number of wavelengths for the reference simulation", 50)
definition.add_optional("npackages", "positive_integer", "number of photon packages per wavelength for the reference simulation", int(1e5))
definition.add_optional("dust_grid_relative_scale", "real", "smallest scale of the dust grid relative to the pixelscale of the input maps", 2.)
definition.add_optional("dust_grid_min_level", "positive_integer", "level for the dust grid", 2)
definition.add_optional("dust_grid_max_mass_fraction", "positive_real", "max mass fraction for the dust grid", 1e-6)

# Flags
definition.add_flag("transient_heating", "enable transient heating", False)
definition.add_flag("selfabsorption", "enable dust selfabsorption", False)

# For remote execution of reference simulation
definition.add_optional("host_ids", "string_list", "remote hosts to use for heavy computations (in order of preference)", choices=find_host_ids(schedulers=False))

# -----------------------------------------------------------------
