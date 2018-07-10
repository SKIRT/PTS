#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import verify_modeling_cwd
from pts.core.simulation.grids import cartesian, bintree, octtree
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
make_seds = "seds"
make_images = "images"
make_choices = [make_seds, make_images]
default_make = [make_seds]
default_make_contributions = False

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

origins = ["model", "fitting_run"]
dust_grid_types = [cartesian, bintree, octtree]
default_dust_grid_type = bintree

# # FITTING RUN
# if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
# elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
# else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# Settings for the wavelength grid
#definition.add_optional("nwavelengths", "integer", "the number of wavelengths to simulate the best model", 450)

# -----------------------------------------------------------------

definition.add_required("origin", "string", "origin of the analysis model", choices=origins)

# -----------------------------------------------------------------

definition.add_required("name", "string_no_spaces", "name for the run")

# -----------------------------------------------------------------

# Remote execution
definition.add_positional_optional("remote", "string", "remote host ID for running the simulation", choices=find_host_ids())
definition.add_flag("attached", "launch remote executions in attached mode", True)

# -----------------------------------------------------------------

# Type of output
definition.add_optional("make", "string_list", "let these forms of output be made", default_make, choices=make_choices)
definition.add_flag("contributions", "make output for the contributions to the flux of the various sources (transparent, dust ...)", default_make_contributions)
definition.add_flag("faceon", "include face-on maps and SEDs", False)
definition.add_flag("edgeon", "include edge-on maps and SEDs", False)

# Full instrument properties
definition.add_flag("scattering_levels", "record scattering levels", 0)
definition.add_flag("counts", "record photon counts for the earth instrument", False)

# To create the output
definition.add_flag("spectral_convolution", "use spectral convolution to calculate observed fluxes", False)

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e5)
definition.add_flag("selfabsorption", "dust self-absorption", False)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", False)

# -----------------------------------------------------------------

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints", "positive_integer", "range of the wavelength grid size", 100)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)
definition.sections["wg"].add_optional("range", "quantity_range", "range of wavelengths", "0.1 micron > 2000 micron", convert_default=True)

# -----------------------------------------------------------------

# # Settings for the dust grid
# definition.add_section("dg", "options for the dust grid")
# definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
# definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 1.)
# definition.sections["dg"].add_optional("min_level", "integer", "the minimum division level for the tree", 8)
# definition.sections["dg"].add_optional("max_mass_fraction", "real", "the maximum mass fraction per cell", 1e-6)

# Dust grid properties
definition.add_section("dg", "settings for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale", "real", "number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 5.) # was 0.5
definition.sections["dg"].add_optional("bintree_min_level", "integer", "minimum depth level for binary trees", 9)
definition.sections["dg"].add_optional("octtree_min_level", "integer", "minimum depth level for octrees", 3)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "maximum mass fraction in each cell", 5e-6)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 15.)

# -----------------------------------------------------------------

# Regenerate
definition.add_flag("regenerate_wavelength_grid", "regenerate the wavelength grid", True)
definition.add_flag("regenerate_dust_grid", "regenerate the dust grid", True)

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", False)

# -----------------------------------------------------------------

# The number of parallel processes for local execution
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")
definition.add_flag("data_parallel_local", "use data-parallelization", False)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", None)

# -----------------------------------------------------------------

# Add optional arguments
#definition.add_optional("remote", "string", "the remote host on which to launch the simulations", "nancy", choices=find_host_ids())
#definition.add_optional("images_remote", "string", "the remote host on which to make the observed images", "nancy", choices=find_host_ids())
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# # Simulation options
# definition.add_optional("npackages", "real", "number of photon packages per wavelength", 1e6)
# definition.add_flag("selfabsorption", "dust self-absorption", True)
# definition.add_flag("transient_heating", "transient (non-LTE) dust heating", True)

# # Parallelization options
# definition.add_optional("nnodes", "integer", "the number of nodes to use for the simulations (for scheduler)", 4)
# definition.add_optional("cores_per_process", "integer", "number of cores per process (for non-scheduler)", 4)
# definition.add_flag("data_parallel", "data parallelization mode", True)

# Special options
#definition.add_flag("dry", "dry-run (don't actually launch simulations)")

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 15)
definition.add_flag("from_projection", "create the projections for the other orientations from the earth projection (instead of from the deprojection model)", True)
definition.add_optional("radial_factor", "real", "factor with which to multiply the radial extent of the projections", 1.5)

# -----------------------------------------------------------------
