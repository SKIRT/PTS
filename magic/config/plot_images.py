# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.magic.config.plot_imagegrid import definition
from pts.magic.plot.imagegrid import default_residual_cmap, default_absolute_residual_cmap, default_direction, directions, observation_name, observation_or_model
from pts.core.basics.plot import diverging_colormaps, normal_colormaps

# ------------------------------------------------------------------------------

columns_name = "columns"
rows_name = "rows"
columns_or_rows = [columns_name, rows_name]

# ------------------------------------------------------------------------------

definition = definition.copy()

# ------------------------------------------------------------------------------

# Dimension of the grid
definition.add_optional("ncolumns", "positive_integer", "number of columns", 6)
definition.add_optional("nrows", "positive_integer", "number of rows", 10)
definition.add_optional("fixed", "string", "which dimension is fixed", columns_name, choices=columns_or_rows)

# ------------------------------------------------------------------------------

definition.add_flag("write_images", "write images", False)
definition.add_flag("write_frames", "write frames", False)
definition.add_flag("write_masks", "write masks", False)
definition.add_flag("write_regions", "write regions", False)

# ------------------------------------------------------------------------------
