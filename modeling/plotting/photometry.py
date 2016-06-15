#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.photometry Contains the PhotometryPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ..core.sed import ObservedSED
from ...core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

class PhotometryPlotter(PlottingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PlottingComponent, self).__init__(config)

        # -- Attributes --

        # The SEDs
        self.seds = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the SEDs
        self.load_seds()

        # 2. Load the aperture data
        #self.load_apertures()

        # 3. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs ...")

        # Load the 'DustPedia' observed SED
        dustpedia_sed = ObservedSED.from_file(self.observed_sed_dustpedia_path)

        # Load the PTS observed SED
        pts_sed = ObservedSED.from_file(self.observed_sed_path)

        # Add the SEDs
        self.seds["DustPedia"] = dustpedia_sed
        self.seds["PTS"] = pts_sed

    # -----------------------------------------------------------------

    def load_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the apertures ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the SEDs
        self.plot_seds()

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs ...")

        # Create the SED plotter
        plotter = SEDPlotter()

        # Add the SEDs
        for label in self.seds: plotter.add_observed_sed(self.seds[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_photometry_path, "seds.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
