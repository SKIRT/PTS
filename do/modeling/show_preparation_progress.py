#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_preparation_progress Show the progress of the data preparation step.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.preparation.preparer import sort_image, status_to_steps, steps
from pts.core.filter.filter import parse_filter
from pts.core.tools import formatting as fmt
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Get configuration
config = parse_arguments("show_preparation_progress", definition)

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

prep_path = fs.join(modeling_path, "prep")

# -----------------------------------------------------------------

progress = defaultdict(dict)

# -----------------------------------------------------------------

# Loop over all images of the initial dataset
#for path in fs.files_in_path(prep_path, recursive=True, exact_name="initialized"):
for directory_path in fs.directories_in_path(prep_path):

    # Directory path
    #directory_path = fs.directory_of(path)
    image_name = fs.name(directory_path)

    # Determine filter
    fltr = parse_filter(image_name)

    # Sort
    label, filepath = sort_image(image_name, directory_path, read_only=True)

    category = fltr.instrument if fltr.instrument is not None else "other"

    # Add to progress data
    progress[category][fltr.band] = label

# -----------------------------------------------------------------

print("")

# Show
for instrument in progress:

    print(fmt.underlined + instrument + fmt.reset + ":")
    print("")

    for band in progress[instrument]:

        print(" " + fmt.blue + band + fmt.reset + ":")
        print(" status = " + progress[instrument][band])
        print("")

        done_steps = status_to_steps(progress[instrument][band])

        for step in steps:

            if step in done_steps: print("   - " + fmt.green + step + ": YES" + fmt.reset) #log.success("  - " + step + ": YES")
            else: print("   - " + fmt.red + step + ": NO" + fmt.reset) #log.error("  - " + step + ": NO")

        print("")

# -----------------------------------------------------------------
