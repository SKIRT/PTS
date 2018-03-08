#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

definition.add_flag("plot", "do plotting", False)

definition.add_optional("fwhm", "positive_real", "FWHM in pixels", 5.)
definition.add_optional("shape", "integer_pair", "shape", (300, 500))

definition.add_optional("nsources", "positive_integer", "number of sources", 100)
definition.add_optional("flux_range", "real_range", "range of flux of sources", "500>1000", convert_default=True)

definition.add_optional("noise_stddev", "real", "stddev of noise", 2.)

definition.add_optional("constant_sky", "real", "constant sky value", 5.)

definition.add_flag("rotate", "rotate", True)
definition.add_optional("rotation_angle", "angle", "rotation angle", "20 deg", convert_default=True)

definition.add_optional("aperture_fwhm_factor", "positive_real", "aperture FWHM factor", 3.0)

definition.add_optional("galaxy_central_flux", "real", "central flux of galaxy", 100.)
definition.add_optional("galaxy_sersic_index", "real", "sersic index", 1.5)
definition.add_optional("galaxy_position", "pixelcoordinate", "galaxy position", "250,150", convert_default=True)
definition.add_optional("galaxy_effective_radius", "positive_real", "galaxy effective radius in pixels", 30.)
definition.add_optional("galaxy_angle", "angle", "galaxy angle", "52 deg", convert_default=True)
definition.add_optional("galaxy_axial_ratio", "positive_real", "axial ratio", 2.0)
definition.add_optional("galaxy_relative_asymptotic_radius", "positive_real", "radius of the galaxy mask relative to the effective radius", 4.)

definition.add_optional("polynomial_degree", "positive_integer", "degree of polynomial", 2)

definition.add_section("plotting", "plotting options")
definition.sections["plotting"].add_flag("meshes", "plot meshes")

# -----------------------------------------------------------------
