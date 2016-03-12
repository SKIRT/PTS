#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.cache Cache ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import tables, inspection, filesystem
from pts.magic.tools import catalogs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("galaxies", type=str, help="the name of the file containing the galaxy catalog", nargs='?', default="galaxies.dat")
parser.add_argument("stars", type=str, help="the name of the file containing the stellar catalog", nargs='?', default="stars.dat")
parser.add_argument("--debug", type=str, help="debug mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

directory_path = os.getcwd()

# Determine the full paths to the galaxy and star catalog
galactic_catalog_path = os.path.join(directory_path, arguments.galaxies)
stellar_catalog_path = os.path.join(directory_path, arguments.stars)

# Open the galactic catalog (to get the name of the principal galaxy)
galactic_catalog = tables.from_file(galactic_catalog_path)

galaxy_name = None

# Loop over the entries of the galactic catalog
for i in range(len(galactic_catalog)):

    if galactic_catalog["Principal"][i]:

        galaxy_name = galactic_catalog["Name"][i]
        break

# If the galaxy name is still None, something is wrong with the galaxy catalog (principal not defined)
if galaxy_name is None: raise RuntimeError("The galactic catalog is invalid: principal galaxy not defined")

# Determine the path to the user catalogs directory
catalogs_user_path = os.path.join(inspection.pts_user_dir, "magic", "catalogs")

# Determint the path to the directory to contain the catalogs for this galaxy
galaxy_user_path = os.path.join(catalogs_user_path, galaxy_name)

# Cache the galaxy and stellar catalog
if os.path.isdir(galaxy_user_path):

    old_galactic_catalog_path = os.path.join(galaxy_user_path, "galaxies.dat")
    old_stellar_catalog_path = os.path.join(galaxy_user_path, "stars.dat")

    if os.path.isfile(old_galactic_catalog_path):

        # Open the 'old' galaxy catalog
        old_galaxy_catalog = tables.from_file(old_galactic_catalog_path)

        # Create merged galaxy catalog
        galaxy_catalog = catalogs.merge_galactic_catalogs(galactic_catalog, old_galaxy_catalog)

        # Save the merged catalog
        path = os.path.join(galaxy_user_path, "galaxies.dat")
        tables.write(galaxy_catalog, path)

    # If a galactic catalog file does not exist yet
    else: filesystem.copy_file(galactic_catalog_path, galaxy_user_path, "galaxies.dat")

    # Check whether there is a stellar catalog file in the galaxy's directory
    if os.path.isfile(old_stellar_catalog_path):

        # Open the new stellar catalog
        stellar_catalog = tables.from_file(stellar_catalog_path)

        # Open the 'old' stellar catalog
        old_stellar_catalog = tables.from_file(old_stellar_catalog_path)

        # Create merged stellar catalog
        stellar_catalog = catalogs.merge_stellar_catalogs(stellar_catalog, old_stellar_catalog)

        # Save the merged catalog
        path = os.path.join(galaxy_user_path, "stars.dat")
        tables.write(stellar_catalog, path)

    # If a stellar catalog file does not exist yet
    else: filesystem.copy_file(stellar_catalog_path, galaxy_user_path, "stars.dat")

else:

    # Create the directory to contain the catalogs for this galaxy
    filesystem.create_directory(galaxy_user_path)

    # Copy the galaxy and stellar catalog files into the new directory
    filesystem.copy_file(galactic_catalog_path, galaxy_user_path, "galaxies.dat")
    filesystem.copy_file(stellar_catalog_path, galaxy_user_path, "stars.dat")

# -----------------------------------------------------------------
