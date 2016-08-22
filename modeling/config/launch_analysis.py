#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Optional settings
definition.add_optional("remotes", "string_list", "the list of remote hosts on which to launch the simulations", ["nancy"], choices=find_host_ids())

# Simulation options
definition.add_optional("packages", "real", "the number of photon packages per wavelength", 1e7)

# Settings for the wavelength grid
definition.add_optional("nwavelengths", "integer", "the number of wavelengths to simulate the best model", 450)

# Settings for the dust grid
definition.add_section("dg", "options for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 1.)
definition.sections["dg"].add_optional("min_level", "integer", "the minimum division level for the tree", 8)
definition.sections["dg"].add_optional("max_mass_fraction", "the maximum mass fraction per cell", 1e-6)

# Parallelization options
definition.add_optional("nodes", "integer", "the number of nodes to use for the simulations", 4)

# -----------------------------------------------------------------
