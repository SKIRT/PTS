#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.collect Collect ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

# Import the relevant PTS classes and modules
from pts.core.basics.loggable import Loggable
from pts.magic.core import Image, Source
from pts.core.tools import filesystem, inspection

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("--debug", type=str, help="debug mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the path to the magic/saturation directory
collector_user_path = os.path.join(inspection.pts_user_dir, "magic", "collector")

saturation_path = os.path.join(collector_user_path, "saturation")
no_saturation_path = os.path.join(collector_user_path, "no saturation")

filesystem.create_directories([collector_user_path, saturation_path, no_saturation_path])

# -----------------------------------------------------------------

class Collector(Loggable):

    """
    This class
    """

    def __init__(self):

        """
        This function ..
        :return:
        """

        # Call the constructor of the base class
        super(Collector, self).__init__()

        # The source
        self.source = None

        # Current index of saturation sources
        self.current_index_yes = -1

        for path in filesystem.files_in_path(saturation_path, extension="fits"):
            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_yes: self.current_index_yes = index

        # Current index of sources without saturation
        self.current_index_no = -1

        for path in filesystem.files_in_path(no_saturation_path, extension="fits"):
            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_no: self.current_index_no = index

        # Setup
        super(Collector, self).setup()

    # -----------------------------------------------------------------

    def save_yes(self, event):

        """
        This function ...
        :return:
        """

        self.current_index_yes += 1

        path = os.path.join(saturation_path, str(self.current_index_yes) + ".fits")

        self.log.info("Saving the saturated star to " + path)

        #source.save(path)

        plt.close()

    # -----------------------------------------------------------------

    def save_no(self, event):

        """
        This function ...
        :return:
        """

        self.current_index_no += 1

        path = os.path.join(no_saturation_path, str(self.current_index_no) + ".fits")

        self.log.info("Saving the star without saturation to " + path)

        #source.save(path)

        plt.close()

# -----------------------------------------------------------------

collector = Collector()

# Loop over all FITS files found in the current directory
for file_path in filesystem.files_in_path(os.getcwd(), extension="fits", contains="saturation"):

    # Get information
    name = os.path.basename(file_path).split(".fits")[0]
    object_type, level, counter = name.split("_")

    # Open the image, select all frames
    image = Image(file_path, always_call_first_primary=False)
    image.frames.select_all()

    # Create a source
    source = Source.from_image(image)

    # Relative peak
    if source.peak is not None:
        rel_peak = source.cutout.rel_position(source.peak)
        peak_coordinates = [[rel_peak.x], [rel_peak.y]]
    else: peak_coordinates = None

    source.plot(title="Is this a saturated star?", show=False, scale="log")

    collector.source = source
    axyes = plt.axes([0.7, 0.05, 0.1, 0.075])
    axno = plt.axes([0.81, 0.05, 0.1, 0.075])
    yes_button = Button(axyes, 'Yes')
    yes_button.on_clicked(collector.save_yes)
    no_button = Button(axno, 'No')
    no_button.on_clicked(collector.save_no)

    # Show the plot
    plt.show()

# -----------------------------------------------------------------
