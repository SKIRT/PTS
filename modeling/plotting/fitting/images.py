#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.chisquared Contains the ChiSquaredPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....magic.plot.imagegrid import ResidualImageGridPlotter
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class ImagesPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(ImagesPlotter, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

        # The simulated images
        self.simulated_images = dict()

        # The observed imags
        self.observed_images = dict()

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Load simulated images
        self.load_simulated_images()

        # Load observed images
        self.load_observed_images()

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Determine the path to the fit/best directory
        best_path = self.fitting_run.best_path

        # Determine the path to the fit/best/images directory
        out_path = fs.join(best_path, "images")

        # Loop over all FITS files found in the fit/best/images directory
        for path, name in fs.files_in_path(out_path, extension="fits", returns=["path", "name"], contains="__"):

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # Get the filter name
            #filter_name = name.split("__")[1]

            # Open the image
            frame = Frame.from_file(path)

            # Get the filter name
            filter_name = str(frame.filter)

            # Add the image frame to the dictionary
            self.simulated_images[filter_name] = frame

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over all FITS files found in the 'truncated' directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Check whether a simulated image exists for this band
            if name not in self.simulated_images:
                log.warning("The simulated version of the " + name + " image could not be found, skipping " + name + " data ...")
                continue

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # The filter name is the image name
            filter_name = name

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.observed_images[filter_name] = frame

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual images ...")

        # Get the ellipse
        path = fs.join(self.truncation_path, "ellipse.reg")
        region = SkyRegion.from_file(path)
        ellipse = region[0]

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Image residuals")

        # Create list of filter names sorted by increasing wavelength
        sorted_filter_names = sorted(self.observed_images.keys(),
                                     key=lambda key: self.observed_images[key].filter.pivotwavelength())

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in sorted_filter_names:
            observed = self.observed_images[filter_name]
            simulated = self.simulated_images[filter_name]

            plotter.add_row(observed, simulated, filter_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(ellipse.bounding_box)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "images.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
