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

definition = definition.copy()

# ------------------------------------------------------------------------------

# For residuals
definition.add_optional("residual_cmap", "string", "colormap for residuals", default_residual_cmap, choices=diverging_colormaps)
definition.add_optional("absolute_residual_cmap", "string", "colormap for absolute residuals", default_absolute_residual_cmap, choices=normal_colormaps)
definition.add_optional("residual_amplitude", "percentage", "amplitude of the residual plots", 1.)
definition.add_optional("residual_interval", "string", "interval for the residual plots", "pts")

# ------------------------------------------------------------------------------

# Writing settings
definition.add_flag("write_observations", "write observation frames", False)
definition.add_flag("write_models", "write model frames", False)
definition.add_flag("write_residuals", "write residual frames", False)
definition.add_flag("write_models", "write model frames", False)
definition.add_flag("write_distributions", "write distributions", False)
definition.add_flag("write_settings", "write plotting settings", False)

# ------------------------------------------------------------------------------

# Sigma-clipping
definition.add_flag("sigma_clip_distributions", "use sigma-clipping on the residual values before creating distributions", True)
definition.add_optional("sigma_clip_level", "positive_real", "sigma level for sigma clipping", 3.)

# ------------------------------------------------------------------------------

# Direction
definition.add_optional("direction", "string", "plotting direction", default_direction, choices=directions)

# ------------------------------------------------------------------------------

# Plot distributions
definition.add_flag("distributions", "plot distributions", False)

# ------------------------------------------------------------------------------

# Method for residuals
definition.add_flag("weighed", "use weighed residuals", False)
definition.add_optional("weighing_reference", "string", "reference for weighed residuals", observation_name, choices=observation_or_model)
definition.add_flag("relative", "use relative residuals", True)
definition.add_flag("absolute", "use the absolute values of the residuals", False)

# ------------------------------------------------------------------------------
