#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition
from pts.modeling.config.calculate_weights import definition as calculate_weights_definition

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# The fitting run for which to adapt the configuration
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# -----------------------------------------------------------------

# Name for the refitting
definition.add_positional_optional("name", "string", "name for the refitting (required when not refitting in-place, or when refitting as new fitting run)")

# Generations?
definition.add_optional("generations", "string_list", "generation names")

# Simulations?
definition.add_optional("simulations", "string_list_or_string_column", "simulation names (in the case of one generation) of filepath with simulation names")
#definition.add_optional("simulation_file", "file_path", "path of file with simulation names") # REPLACED BY 'STRING_LIST_OR_STRING_COLUMN' ABOVE
#definition.add_optional("simulation_file_column", "positive_integer", "index of the column containing the simulation names")

# -----------------------------------------------------------------

# Refit as a new fitting run (copy current fitting run)
definition.add_optional("as_run", "string", "create new fitting run", forbidden=runs.names)
definition.add_optional("in_place", "string", "change everything in place, but make a backup with the specified name")

# -----------------------------------------------------------------

# New fitting filters
definition.add_optional("filters", "filter_list", "filters to use for the evaluation (None means default fitting filters)")
definition.add_optional("not_filters", "lazy_filter_list", "filters to ignore for the fitting")

# -----------------------------------------------------------------

# Flags
definition.add_flag("show", "show", True)
definition.add_flag("plot", "make plots", True)

# -----------------------------------------------------------------

# Redo
definition.add_flag("reweigh", "recalculate weights", None)
definition.add_flag("reflux", "recalculate mock observed fluxes", None)
definition.add_flag("rediff", "recalculate differences", None)

# -----------------------------------------------------------------

# Flux calculation options
definition.add_flag("spectral_convolution", "enable spectral convolution for calculating fluxes (None means default for generation)", None)
definition.add_flag("fluxes_from_images", "use images to calculate observed fluxes (None means default for generation)", None)

# -----------------------------------------------------------------

# For weights
definition.import_section("weighing", "wavelength filter weighing", calculate_weights_definition, pos_optional=False)

# -----------------------------------------------------------------

definition.add_optional("nbest", "positive_integer", "number of best simulations to show", 10)
definition.add_flag("show_best_sed", "plot the SED of the best simulation", False)
definition.add_flag("show_counts", "show the counts of different parameter values in the best simulations", True)

# -----------------------------------------------------------------

# Plotting flags
definition.add_flag("plot_counts", "plot the counts", True)
definition.add_flag("plot_chi_squared", "plot the chi squared distributions", True)
definition.add_flag("plot_probabilities", "plot probabilities", True)
definition.add_flag("plot_distributions", "plot distributions", True)

# -----------------------------------------------------------------

# Additional relative error
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")

# -----------------------------------------------------------------
