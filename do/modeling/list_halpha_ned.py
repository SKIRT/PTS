#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.list_halpha_ned This ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astroquery.ned import Ned

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.filter import Filter

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
definition.add_required("galaxy", "string", "galaxy name")

# Get configuration
setter = ArgumentConfigurationSetter("check_simulated_images")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting list_halpha_ned ...")

# -----------------------------------------------------------------

# Get the list
urls = Ned.get_image_list(config.galaxy)

images = []

# Print the list
for url in urls:

    # Get the name
    name = fs.strip_extension(fs.strip_extension(fs.name(url))) # strip both the .gz as the .fits extension

    # Get the bibcode
    try: bibcode = url.split("img/")[1].split("/")[0]
    except IndexError: bibcode = None

    if ":" in name:

        splitted = name.split(":")

        if splitted[0].startswith("NGC_"):
            band = splitted[0].split("NGC_")[1][5:]
            try:
                filter = Filter.from_string(band)
                splitted = [config.galaxy, None, band, splitted[1]]
            except: pass

        if len(splitted) == 3:

            splitted = [config.galaxy, None, splitted[1], splitted[2]]

        elif len(splitted) == 2:

            info_and_band = splitted[0].split("NGC_")[1][5:]
            splitted = [config.galaxy, None, info_and_band, splitted[1]]

        galaxy_name = splitted[0]
        unknown = splitted[1]
        band = splitted[2]
        source = splitted[3]

        try:
            year = int(source[-4:])
            if year < 1985: continue
        except ValueError: year = None

        images.append((band, year, bibcode, url))

    elif "_" in name:

        splitted = name.split("_")

        band = splitted[-1]

        images.append((band, None, bibcode, url))

    elif "." in name:

        splitted = name.split(".")

        galaxy_name = splitted[0]

        images.append((None, None, bibcode, url))

# Print
for band, year, bibcode, url in images:

    if band is not None and ("Ha" in band or "H-alpha" in band or "H_alph" in band):

        print(band, year, bibcode, url)

# -----------------------------------------------------------------
