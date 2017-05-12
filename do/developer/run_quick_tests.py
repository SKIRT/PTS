#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.run_quick_tests Run quick versions of some PTS tests.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.launch.pts import launch_local

# -----------------------------------------------------------------

pts_command = "tests"

# -----------------------------------------------------------------

sources_test_settings = dict()
sources_test_settings["rotate"] = False
sources_test_settings["nrandom_sources"] = 200
sources_test_settings["vary_fwhm"] = False
sources_test_settings["add_catalogued_sources"] = False
sources_test_settings["nprocesses"] = 1
#sources_test_settings["remote"] = "nancy"

# -----------------------------------------------------------------

settings = dict()
settings["subprojects"] = ["magic"]
settings["tests"] = ["Sources"]
settings["debug"] = True
settings["keep"] = True
settings["remove_previous"] = True
settings["settings"] = sources_test_settings

# -----------------------------------------------------------------

# Launch
launch_local(pts_command, settings)

# -----------------------------------------------------------------

# Other:
# pts tests evolve Stepwise --debug --keep --default --settings "'plot':True" --open_output --remove_previous
# pts tests magic Sources --debug --default --keep --remove_previous
# pts tests modeling M81_sed --settings "'free_parameters':['dust_mass']" --debug --keep

# pts tests modeling M81_sed --debug --settings "'free_parameters':['dust_mass'], 'reference_test':'M81_sed_2017-05-06--21-01-58-593', 'ngenerations': 10, 'nsimulations': 10" --keep --only_tests

# pts tests modeling M81_sed --debug --settings "'free_parameters':['dust_mass'], 'reference_test':'M81_sed_2017-05-06--21-01-58-593', 'ngenerations': 3, 'nsimulations': 5" --keep --only_tests

# -----------------------------------------------------------------
