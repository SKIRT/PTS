#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.modelgalaxy Model a galaxy with Astromagic and SKIRT
#

# *****************************************************************

# Import standard modules
import os.path
import argparse
from datetime import datetime

# Import relevant PTS modules
from modeling.galaxymodeler import GalaxyModeler

# *****************************************************************

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument('--stage', type=str, help='the preparation stage')
parser.add_argument('--image', type=str, help='provide this argument to only prepare one specific image')
parser.add_argument('--report', action='store_true', help='write a report file')
parser.add_argument('--plot', action='store_true', help='plot the result of intermediate steps')

# Parse the command line arguments
args = parser.parse_args()

# *****************************************************************

# Get the path to the current working directory
working_directory = os.getcwd()

# *****************************************************************

# Create a GalaxyModeler object
modeler = GalaxyModeler(working_directory, args.config)

# Set configuration options passed as command line arguments
modeler.config.preparation.filter_name = args.image
if args.report:

    # Determine a unique report path and set the appropriate configuration entry
    timestamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
    report_path = os.path.join(working_directory, "report_" + timestamp + ".txt")
    modeler.config.logging.path = report_path

# *****************************************************************

# Run the modeling (or a specific stage)
if args.stage is None: modeler.run()
elif args.stage == "preparation": modeler.prepare_images()
elif args.stage == "decomposition": modeler.decompose()
elif args.stage == "mapmaking": modeler.make_maps()
elif args.stage == "fitting": modeler.fit_sed()
else: raise ValueError("Unkown stage (choose 'preparation', 'decomposition', 'mapmaking' or 'fitting')")

# *****************************************************************
