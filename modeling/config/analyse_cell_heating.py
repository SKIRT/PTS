#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import warnings

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.config.component import definition
from pts.modeling.analysis.heating.cell import disk_components, ionizing

# -----------------------------------------------------------------

default_disk_component = ionizing

# -----------------------------------------------------------------

default_spacing_measure = "min"
spacing_measures = ["min", "max", "mean", "median"]

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# The analysis run
if runs.empty: warnings.warn("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.names)

# -----------------------------------------------------------------

# The number of bins
definition.add_optional("nbins", "positive_integer", "number of bins", 20)
definition.add_optional("nradial_bins", "positive_integer", "number of radial bins", 200)

# For the midplane heating
definition.add_optional("midplane_component", "string", "disk component of which to use the scaleheight as the reference for defining the midplane height", default_disk_component, choices=disk_components)
definition.add_optional("midplane_factor", "positive_real", "factor to be multiplied with the component scaleheight to define the height of the midplane", 0.5)

# For the scale of the heating map
definition.add_optional("map_spacing_measure", "string", "measure to be used to determined the spacing of the heating map w.r.t. to the spacing of the cell coordinates", default_spacing_measure, choices=spacing_measures)
definition.add_optional("map_spacing_factor", "positive_real", "factor to be multiplied with the spacing determined from the average of the x and y spacing measure", 5.)

# For the scale of the midplane heating map
definition.add_optional("midplane_spacing_measure", "string", "measure to be used to determine the spacing of the midplane heating map w.r.t. the spacing of the cell coordinates in the midplane", default_spacing_measure, choices=spacing_measures)
definition.add_optional("midplane_spacing_factor", "positive_real", "factor to be multiplied with the spacing determined from the average of the x and y spacing measure", 2.)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "do plotting", True)
definition.add_flag("plot_distribution", "plot distribution", True)
definition.add_flag("plot_radial_distribution", "plot radial distribution", True)
definition.add_flag("plot_map", "plot map", True)
definition.add_flag("plot_map_midplane", "plot map of the midplane", True)

# -----------------------------------------------------------------

# Recreate table
definition.add_flag("recreate_table", "recreate the absorption table")

# Recalculate fractions, distributions
definition.add_flag("recalculate_fractions", "recalculate the heating fractions")
definition.add_flag("recalculate_distributions", "recreate the heating distributions")
definition.add_flag("recalculate_distribution", "recalculate")
definition.add_flag("recalculate_distribution_diffuse", "recalculate")
definition.add_flag("recalculate_radial_distribution", "recalculate")

# Recreate maps
definition.add_flag("recreate_maps", "recreate the map and the map in the midplane")
definition.add_flag("recreate_map", "recreate the map")
definition.add_flag("recreate_map_midplane", "recreate the map in the midplane")

# Re-interpolate maps
definition.add_flag("reinterpolate_maps", "recreate the interpolated maps")
definition.add_flag("reinterpolate_map", "recreate the interpolated map")
definition.add_flag("reinterpolate_map_midplane", "recreate the interpolated map in the midplane")
definition.add_optional("not_nans_dilation_radius", "positive_real", "radius for dilating not-nans", 3)

# Replot
definition.add_flag("replot", "replot")
definition.add_flag("replot_distribution", "replot the distribution")
definition.add_flag("replot_radial_distribution", "replot the radial distribution")
definition.add_flag("replot_map", "replot the map")
definition.add_flag("replot_map_midplane", "replot the map in the midplane")
definition.add_flag("replot_maps", "replot the map and the map in the midplane")

# -----------------------------------------------------------------

# Keep consistent
definition.add_flag("consistency", "assure heating fraction maps are consistent with ncells and stddev maps", True)

# -----------------------------------------------------------------

default_format = "pdf"
formats = ["png", "pdf"]

# -----------------------------------------------------------------

# Plotting
# M81: radii = [2.22953938405, 3.34430907608, 5.90827936773, 8.9181575362]  # kpc
definition.add_optional("radial_distribution_radii", "quantity_list", "radii at which to draw vertical lines on the radial distribution plot")
definition.add_optional("plotting_format", "string", "plotting format", default_format, choices=formats)
definition.add_optional("nlevels", "positive_integer", "number of levels for contour plots of the heating fraction maps", 5)
definition.add_flag("contours", "add contours to the plots")

# -----------------------------------------------------------------

# For creating interpolated maps
definition.add_optional("min_ncells", "positive_integer", "minimum number of cells for eachn pixel of the interpolated map of the heating fraction", 10)
definition.add_optional("min_ncells_midplane", "positive_integer", "minimum number of cells for each pixel of the interpolated map of the heating fraction in the midplane", 5)

# -----------------------------------------------------------------
