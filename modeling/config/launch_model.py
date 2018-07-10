#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.core.remote.host import find_host_ids
from pts.core.simulation.grids import cartesian, bintree, octtree
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
make_seds = "seds"
make_images = "images"
make_choices = [make_seds, make_images]
#default_make = make_choices
default_make = [make_seds]
default_make_contributions = False

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

origins = ["model", "fitting_run"]
dust_grid_types = [cartesian, bintree, octtree]
default_dust_grid_type = bintree

# -----------------------------------------------------------------

default_nwavelengths = 100
default_wavelength_range = "0.02 micron > 2000 micron"

# -----------------------------------------------------------------

definition.add_required("origin", "string", "origin of the analysis model", choices=origins)

# -----------------------------------------------------------------

definition.add_required("name", "string_no_spaces", "name for the simulation")
definition.add_flag("clear", "clear launch directory when existing")

# -----------------------------------------------------------------

# Remote execution
definition.add_positional_optional("remote", "string", "remote host ID for running the simulation", choices=find_host_ids())
definition.add_flag("attached", "launch remote executions in attached mode", True)

# Keep?
definition.add_flag("keep_remote_input", "keep the remote input directory after the simulation is retrieved", False)
definition.add_flag("keep_remote_input_and_output", "keep the remote input and output after the simulation is retrieved", False)

# -----------------------------------------------------------------

# Type of output
definition.add_optional("make", "string_list", "let these forms of output be made", default_make, choices=make_choices)
definition.add_flag("contributions", "make output for the contributions to the flux of the various sources (transparent, dust ...)", default_make_contributions)
definition.add_flag("faceon", "include face-on maps and SEDs", True)
definition.add_flag("edgeon", "include edge-on maps and SEDs", True)

# USE IMAGES TO MAKE SEDs?
definition.add_flag("make_image_seds", "make SEDs based on the (masked) output images", False)

# Full instrument properties
definition.add_optional("scattering_levels", "positive_integer" ,"record scattering levels", 0)
definition.add_flag("counts", "record photon counts for the earth instrument", False)

# To create the output
definition.add_flag("spectral_convolution", "use spectral convolution to calculate observed fluxes", False)

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 1e5)
definition.add_flag("selfabsorption", "dust self-absorption", False)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", False)

# Adjust the number of packages?
definition.add_flag("adjust_npackages", "adjust the number of photon packages to the dust grid and instruments", True)

# -----------------------------------------------------------------

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints", "positive_integer", "range of the wavelength grid size", default_nwavelengths)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", False)
definition.sections["wg"].add_optional("range", "quantity_range", "range of wavelengths", default_wavelength_range, convert_default=True)

# -----------------------------------------------------------------

# Dust grid properties
definition.add_section("dg", "settings for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale", "real", "number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 2.) # was 0.5
definition.sections["dg"].add_optional("bintree_min_level", "integer", "minimum depth level for binary trees", 9)
definition.sections["dg"].add_optional("octtree_min_level", "integer", "minimum depth level for octrees", 3)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "maximum mass fraction in each cell", 1e-6)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 15.)

# -----------------------------------------------------------------

# Regenerate
definition.add_flag("regenerate_wavelength_grid", "regenerate the wavelength grid", True)
definition.add_flag("regenerate_dust_grid", "regenerate the dust grid", True)

# -----------------------------------------------------------------

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", False)

# -----------------------------------------------------------------

# The number of parallel processes
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes for local execution", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes for remote execution")
definition.add_flag("data_parallel_local", "use data-parallelization", False)
definition.add_flag("data_parallel_remote", "use data-parallelization for remote execution", None)

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 15)
definition.add_flag("from_projection", "create the projections for the other orientations from the earth projection (instead of from the deprojection model)", True)
definition.add_optional("radial_factor", "real", "factor with which to multiply the radial extent of the projections", 1.5)

# -----------------------------------------------------------------
