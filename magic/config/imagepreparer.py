#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add options
config.add_optional("error_frame_names", "string_list", "the names of the error planes to be included in the final error map", [])
config.add_flag("write_steps", "write the results of intermediate steps")
config.add_flag("write_sky_annuli", "write the sky annuli")
config.add_optional("sky_annuli_path", str, "the path to the sky annuli directory", None)
config.add_optional("calculate_poisson_noise", bool, "calculate poisson noise", True)
config.add_optional("calculate_calibration_uncertainties", bool, "calculate calibration uncertanties", True)
config.add_optional("extract_sources", bool, "extract sources", True)
#config.add_topic("extraction")
#config.topics["extraction"].add_...

# Star extraction
extract_sources: True
extraction:
{
}

# Galactic extinction
correct_for_extinction: True

# The galactic attenuation
attenuation: None

# Unit conversion
convert_unit: True
unit_conversion:
{
  to_unit: "MJy/sr"
}

# Convolution
convolve: True
convolution:
{
  kernel_path: None
  kernel_fwhm: None
  remote: None
}

# Rebinning
rebin: True
rebinning:
{
  reference_path: None
}

# Sky subtraction
subtract_sky: True
sky_subtraction:
{
}

# Uncertainties
set_uncertainties: True
uncertainties:
{
  calibration_error: None
}

# Cropping
crop: True
cropping:
{
  limits: None # e.g. [350, 725, 300, 825]
}
