#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.tests.base import possible_free_parameters, default_free_parameters
from pts.evolve.solve.extremizer import genetic_definition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Optional settings
definition.add_optional("nwavelengths", "positive_integer", "number of wavelengths for the reference simulation", 100)
definition.add_optional("npackages", "positive_integer", "number of photon packages per wavelength for the reference simulation", int(1e4))
definition.add_optional("dust_grid_relative_scale", "real", "smallest scale of the dust grid relative to the pixelscale of the input maps", 10.)
definition.add_optional("dust_grid_min_level", "positive_integer", "level for the dust grid", 0)
definition.add_optional("dust_grid_max_mass_fraction", "positive_real", "max mass fraction for the dust grid", 1e-5)

# Flags
definition.add_flag("transient_heating", "enable transient heating", False)
definition.add_flag("selfabsorption", "enable dust selfabsorption", False)

# For remote execution
definition.add_optional("host_ids", "string_list", "remote hosts to use for heavy computations and simulations", choices=find_host_ids())
definition.add_flag("attached", "launch remote executions in attached mode", False)

# Fitting
definition.add_optional("ngenerations", "positive_integer", "number of generations", 2)
definition.add_optional("nsimulations", "even_positive_integer", "number of simulations per generation", 2)
definition.add_optional("npackages_fitting", "positive_integer", "number of photon packages for each fitting simulation", int(1e3))

# Free parameters
definition.add_optional("free_parameters", "string_list", "free parameter labels", choices=possible_free_parameters, default=default_free_parameters)
definition.add_optional("relative_range_initial", "real_range", "relative range for generating the initial parameter values", default="0.3>3", convert_default=True)
definition.add_optional("relative_range_fitting", "real_range", "relative range of the free parameter values for the fitting", default="0.1>20", convert_default=True)

# Wavelength grid
definition.add_optional("wavelength_range", "quantity_range", "range of wavelengths", "0.1 micron > 2000 micron", convert_default=True)

# Dust grid
definition.add_optional("physical_domain_disk_ellipse_factor", "positive_real", "factor of the disk ellipse region to take as the physical domain of the model", 0.82)

# Genetic section
definition.import_section("genetic", "genetic algorithm options", genetic_definition)
#definition.sections["genetic"].optional["crossover_method"].default = ""
#definition.sections["genetic"].optional["nelite_individuals"].default = 2
definition.sections["genetic"].optional["mutation_rate"].default = 0.07
definition.sections["genetic"].optional["genome_type"].default = "binary_string"
#definition.sections["genetic"].flags["gray_code"].default = True # IS ALREADY TRUE BY DEFAULT

# Flags
definition.add_flag("spectral_convolution", "use spectral convolution to calculate observed fluxes", False)
definition.add_flag("cheat", "cheat by putting the real parameter values into the initial population of models", False)

default_fitting_method = "genetic"
fitting_methods = ["genetic", "grid"]
definition.add_optional("fitting_method", "string", "fitting method", default_fitting_method, fitting_methods)

definition.add_optional("reference_path", "directory_path", "use the simulation in this directory as the reference simulation and infer the real parameter values from it")
definition.add_optional("reference_test", "string", "use the reference simulation of this previous test")

# Scale
default_scale = "logarithmic"
scales = ["logarithmic", "linear"]
definition.add_optional("scale", "string", "scale to use for the generation of grid/genetic individuals", default_scale, scales)

# Recurrence
# Check recurrence
definition.add_flag("check_recurrence", "check for recurrence of models that have been simulated previously", True)
definition.add_optional("recurrence_rtol", "positive_real", "relative tolerance for recurrence checking", 1e-5)
definition.add_optional("recurrence_atol", "positive_real", "absolute tolerance for recurrence checking", 1e-8)

# -----------------------------------------------------------------
