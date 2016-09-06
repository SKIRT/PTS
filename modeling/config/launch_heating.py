#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids
from pts.core.tools import filesystem as fs
from pts.modeling.analysis.component import get_analysis_run_names, get_last_run_name

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Positional option
definition.add_positional_optional("run", "string", "name of the analysis run for which to launch the heating simulations", get_last_run_name(fs.cwd()), get_analysis_run_names(fs.cwd()))

# Optional settings
definition.add_optional("remote", "string", "remote host on which to launch the simulations", "nancy", choices=find_host_ids())

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e7)

# Settings for the wavelength grid
definition.add_section("wg", "options for the wavelength grid")
definition.sections["wg"].add_optional("range", "quantity_range", "the wavelength range", "0.1 micron>10micron", convert_default=True)
definition.sections["wg"].add_optional("npoints", "integer", "the number of wavelength points", 25)

# Settings for the dust grid
definition.add_section("dg", "options for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 1.)
definition.sections["dg"].add_optional("min_level", "integer", "the minimum division level for the tree", 8)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "the maximum mass fraction per cell", 1e-6)

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations", 4)

# -----------------------------------------------------------------
