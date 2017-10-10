#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.basics.models import models_3D

# -----------------------------------------------------------------

default_npackages = 1e7
default_parallelization = "2:1:2" # 2 cores, 1 process, 2 threads per core

# -----------------------------------------------------------------

default_deprojection_method = "skirt"
deprojection_methods = ["pts", "skirt"]

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

definition.add_required("name", "string", "name for the model")
definition.add_required("model_type", "string", "type of model", choices=models_3D)

# -----------------------------------------------------------------

definition.add_optional("map", "file_path", "path of the input map (for deprojection model)")

# -----------------------------------------------------------------

definition.add_flag("project", "make projections", False)
definition.add_flag("deproject", "make deprojections", False)
definition.add_flag("view", "create 3D view", False)

# Deprojection method
definition.add_optional("deprojection_method", "string", "deprojection method", default=default_deprojection_method, choices=deprojection_methods)

# Show the view
definition.add_flag("show", "show the 3D view at the end", False)

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("scale_heights", "real", "number of times to take the scale height as the vertical radius of the model", 15.)

# -----------------------------------------------------------------

# SKIRT options
definition.add_optional("npackages", "positive_integer", "number of photon packages", default_npackages)
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations", default_parallelization, convert_default=True)

# -----------------------------------------------------------------

# Output
definition.add_optional("output", "directory_path", "output directory")

# -----------------------------------------------------------------

# Additional properties specified through the command line (not necessary)
definition.add_optional("center", "sky_or_pixel_coordinate", "coordinate of the galaxy center")
definition.add_optional("distance", "length_quantity", "galaxy distance")
definition.add_optional("position_angle", "angle", "position angle of the galaxy")
definition.add_optional("inclination", "angle", "inclination angle of the galaxy")
definition.add_optional("scale_height", "length_quantity", "vertical scale height")

# -----------------------------------------------------------------
