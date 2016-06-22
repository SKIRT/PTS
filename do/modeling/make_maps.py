#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.makemaps Make the maps of dust and stars for the SKIRT radiative transfer modeling procedure.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.maps.mapmaking import MapMaker
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add optional settings
config.add_optional("map", str, "the map to be made (dust, old, NIY, IY)")

config.add_section("cutoff")
config.sections["cutoff"].add_optional("reference_path", str, "...", None)
config.sections["cutoff"].add_optional("level", float, "cutoff when signal < level * uncertainty (ilse: 5)", 3.0)
config.sections["cutoff"].add_optional("remove_holes", bool, "remove holes from the cutoff mask", True)

config.add_section("dust")
config.sections["dust"].add_section("ssfr")
config.sections["dust"].sections["ssfr"].add_optional("mask_low_fuv_snr", bool, "...", True)
config.sections["dust"].sections["ssfr"].add_optional("fuv_snr_level", float, "cut-off when signal(FUV) < fuv_snr_level * uncertainty(FUV)  (Ilse: 10.0)", 0.0)

config.add_section("old_stars")
config.sections["old_stars"].add_optional("irac_snr_level", float, "cut-off when signal(IRAC) < irac_snr_level * uncertainty(IRAC) (Ilse: 10.0)", 0.0)

config.add_section("ionizing_stars")
config.sections["ionizing_stars"].add_section("mips_young_stars")
config.sections["ionizing_stars"].sections["mips_young_stars"].add_optional("mips_snr_level", float, "cut-off when signal(MIPS) < mips_snr_level * uncertainty(MIPS)  (Ilse: 10.0)", 0.0)
config.sections["ionizing_stars"].add_optional("ha_snr_level", float, "cut-off ionizing stars map when signal(Ha) < ha_snr_level * uncertainty(Ha) (Ilse: 10.0)", 0.0)
config.sections["ionizing_stars"].add_optional("mips_snr_level", float, "cut-off ionizing stars map when signal(24 micron) < mips_snr_level * uncertainty(24 micron) (Ilse: 10.0)", 0.0)

config.add_section("non_ionizing_stars")
config.sections["non_ionizing_stars"].add_optional("fuv_snr_level", float, "cut-off when signal(FUV) < fuv_snr_level * uncertainty(FUV) (Ilse: 10.0)", 0.0)

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting make_maps ...")

# -----------------------------------------------------------------

# Create a MapMaker object
maker = MapMaker(config.get_settings())

# Run the map maker
maker.run()

# -----------------------------------------------------------------
