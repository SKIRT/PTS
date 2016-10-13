#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids
from pts.modeling.fitting.component import get_finished_generations, get_last_finished_generation
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Set the default option for the generation name
last_generation_name = get_last_finished_generation(fs.cwd())
if last_generation_name is None: last_generation_name = "first_guess"

# Set the choices for the generationn name
generation_names = ["first_guess"] + get_finished_generations(fs.cwd())

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Required settings
definition.add_positional_optional("generation", "string", "the name of the (finished) generation for which to relaunch the best simulation", default=last_generation_name, choices=generation_names)

# Settings for the wavelength grid
definition.add_optional("nwavelengths", "integer", "the number of wavelengths to simulate the best model", 450)

# Settings for the dust grid
definition.add_section("dg", "options for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 1.)
definition.sections["dg"].add_optional("min_level", "integer", "the minimum division level for the tree", 8)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "the maximum mass fraction per cell", 1e-6)

# Add optional arguments
definition.add_optional("remote", "string", "the remote host on which to launch the simulations", "nancy", choices=find_host_ids())
definition.add_optional("images_remote", "string", "the remote host on which to make the observed images", "nancy", choices=find_host_ids())
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# Simulation options
definition.add_optional("npackages", "real", "number of photon packages per wavelength", 1e6)
definition.add_flag("selfabsorption", "dust self-absorption", True)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", True)

# Parallelization options
definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations (for scheduler)", 4)
definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 4)
definition.add_flag("data_parallel", "data parallelization mode", True)

# Special options
definition.add_flag("dry", "dry-run (don't actually launch simulations)")

# -----------------------------------------------------------------
