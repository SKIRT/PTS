#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.tests.base import default_free_parameters, possible_free_parameters
from pts.evolve.solve.extremizer import genetic_definition

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

# For remote execution
definition.add_optional("host_ids", "string_list", "remote hosts to use for heavy computations and simulations", choices=find_host_ids())
definition.add_flag("attached", "launch remote executions in attached mode", False)

# Fitting
definition.add_optional("ngenerations", "positive_integer", "number of generations", 5)
definition.add_optional("nsimulations", "even_positive_integer", "number of simulations per generation", 30)
definition.add_optional("npackages_fitting", "positive_integer", "number of photon packages for each fitting simulation", int(1e4))

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

# Flags
definition.add_flag("spectral_convolution", "use spectral convolution to calculate observed fluxes", False)

# -----------------------------------------------------------------
