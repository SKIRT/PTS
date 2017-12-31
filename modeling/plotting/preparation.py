#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.preparation Contains the PreparationPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ..preparation.component import PreparationComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...magic.core.fits import get_frame_names
from ...magic.basics.mask import Mask, get_mask_names
from ...magic.region.list import PixelRegionList
from ...magic.plot.imagegrid import StandardImageGridPlotter
from ...core.plot.distribution import DistributionGridPlotter, DistributionPlotter
from ...core.basics.distribution import Distribution
from ...magic.plot.error import ErrorPlotter
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class PreparationPlotter(PlottingComponent, PreparationComponent):
    
    """
    This class...
    """

    # The load functions
    load_functions = dict()

    # The plot functions
    plot_functions = dict()

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructors of the base classes
        PlottingComponent.__init__(self, no_config=True)
        PreparationComponent.__init__(self, *args, **kwargs)

        # -- Attributes --

        # Features to plot
        self.features = None

        # The paths to the resulting FITS files
        self.result_paths = dict()

        # The paths to the sky directories
        self.sky_paths = dict()

        # The dictionary of prepared image frames
        self.images = dict()

        # The dictionary of error frames
        self.errors = dict()
        self.poisson_errors = dict()
        self.calibration_errors = dict()
        self.sky_errors = dict()

        # The dictionary of sources masks
        self.sources_masks = dict()

        # The dictionary of sky masks
        self.sky_masks = dict()

        # The dictionary of sky values
        self.sky_values = dict()

        # The dictionary of sky annuli
        self.annuli = dict()

        # The dictionary of sky apertures
        self.apertures = dict()

    # -----------------------------------------------------------------

    def run(self, features=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(features)

        # 2. Load the prepared images
        self.load_images()

        # 3. Load the error frame
        self.load_errors()

        # 4. Load the source and sky masks
        self.load_masks()

        # 5. Load the sky values
        self.load_sky()

        # 6. Load the galaxy and sky annuli
        self.load_annuli()

        # 7. Load the sky apertures
        self.load_apertures()

        # 8. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, features=None):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PreparationPlotter, self).setup()

        # Set features to plot
        self.features = features

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Add the image path to the dictionary
            self.result_paths[directory_name] = image_path

            # Look for the 'sky' directory
            sky_path = fs.join(directory_path, "sky")
            if not fs.is_directory(sky_path):
                log.warning("Sky directory is not present for " + directory_name)
                continue

            # Add the sky directory path to the dictionary
            self.sky_paths[directory_name] = sky_path

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prepared images ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the prepared image frame
            frame = Frame.from_file(self.result_paths[label])

            # Set the image name
            frame.name = label

            # Add the image to the dictionary
            self.images[label] = frame

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Load sources masks
        self.load_sources_masks()

        # Load sky masks
        self.load_sky_masks()

    # -----------------------------------------------------------------

    def load_sources_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sources masks ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Check whether the sources mask is present in the FITS file
            if not "sources" in get_mask_names(self.result_paths[label]):
                log.warning("The sources mask is not present in the " + label + " image")

            # Open the sources mask
            mask = Mask.from_file(self.result_paths[label], plane="sources")

            # Add the mask to the dictionary
            self.sources_masks[label] = mask

    # -----------------------------------------------------------------

    def load_sky_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky masks ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Check whether the sky mask is present in the FITS file
            if not "sky" in get_mask_names(self.result_paths[label]):
                log.warning("The sky mask is not present in the " + label + " image")
                continue

            # Open the sky mask
            mask = Mask.from_file(self.result_paths[label], plane="sky")

            # Add the sky mask to the dictionary
            self.sky_masks[label] = mask

    # -----------------------------------------------------------------

    def load_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky values ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the sky frame
            sky = Frame.from_file(self.result_paths[label], plane="sky")

            # Get the sky value (assuming the sky frame is constant)
            value = sky[0,0]

            # Add the sky value to the dictionary
            self.sky_values[label] = value

    # -----------------------------------------------------------------

    def load_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy and sky annuli ...")

        # Loop over the sky paths
        for label in self.sky_paths:

            # Look for the annulus region file
            region_path = fs.join(self.sky_paths[label], "annulus.reg")
            if not fs.is_file(region_path):
                log.warning("The annulus region could not be found for " + label)
                continue

            # Open the annulus region
            region = PixelRegionList.from_file(region_path).homogenized()

            # Add the region to the dictionary
            self.annuli[label] = region

    # -----------------------------------------------------------------

    def load_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky apertures ...")

        # Loop over the sky paths
        for label in self.sky_paths:

            # Look for the apertures FITS file
            apertures_path = fs.join(self.sky_paths[label], "apertures.fits")
            if not fs.is_file(apertures_path):
                log.warning("The apertures image could not be found for " + label)
                continue

            # Open the apertures image
            apertures = Frame.from_file(apertures_path)

            # Add the apertures image to the dictionary
            self.apertures[label] = apertures

    # -----------------------------------------------------------------

    def load_errors(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Loading the error frames ...")

        # Load the total errors
        self.load_total_errors()

        # Load the poisson errors
        self.load_poisson_errors()

        # Load the sky errors
        self.load_sky_errors()

        # Load the calibration errors
        self.load_calibration_errors()

    # -----------------------------------------------------------------

    def load_total_errors(self):

        """
        This function ...
        :return:
        """

        # Loop over the image paths
        for label in self.result_paths:

            # Open the errors frame
            frame = Frame.from_file(self.result_paths[label], plane="errors")

            # Set the image name
            frame.name = label

            # Add the error frame to the dictionary
            self.errors[label] = frame

    # -----------------------------------------------------------------

    def load_poisson_errors(self):

        """
        This function ...
        :return:
        """

        # Loop over the image paths
        for label in self.result_paths:

            # Check if the poisson_errors frame is present in the FITS file
            if not "poisson_errors" in get_frame_names(self.result_paths[label]): continue

            # Open the poisson errors frame
            errors = Frame.from_file(self.result_paths[label], plane="poisson_errors")

            # Add the error frame to the dictionary
            self.poisson_errors[label] = errors

    # -----------------------------------------------------------------

    def load_sky_errors(self):

        """
        This function ...
        :return:
        """

        # Loop over the image paths
        for label in self.result_paths:

            # Check if the sky_errors frame is present in the FITS file
            if not "sky_errors" in get_frame_names(self.result_paths[label]):
                log.warning("The sky_errors frame is not present in the " + label + " image")
                continue

            # Open the sky error frame
            errors = Frame.from_file(self.result_paths[label], plane="sky_errors")

            # Add the error frame to the dictionary
            self.sky_errors[label] = errors

    # -----------------------------------------------------------------

    def load_calibration_errors(self):

        """
        This function ...
        :return:
        """

        # Loop over the image paths
        for label in self.result_paths:

            # Check if the calibration_errors frame is present in the FITS file
            if not "calibration_errors" in get_frame_names(self.result_paths[label]):
                log.warning("The calibration_errors frame is not present in the " + label + " image")
                continue

            # Open the calibration error frame
            errors = Frame.from_file(self.result_paths[label], plane="calibration_errors")

            # Add the error frame to the dictionary
            self.calibration_errors[label] = errors

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot a grid of the prepared images
        if self.features is None or "images" in self.features: self.plot_images()

        # Plot the grid of images with the sources masks and sky annuli overlayed
        if self.features is None or "masks_annuli" in self.features: self.plot_masks_and_annuli()

        # Plot a grid of the apertures
        if self.features is None or "apertures" in self.features: self.plot_apertures()

        # Plot the sky values
        if self.features is None or "sky" in self.features: self.plot_sky()

        # Plot the distributions of the relative errors
        if self.features is None or "errors" in self.features: self.plot_errors()

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in self.sorted_labels: plotter.add_image(self.images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "preparation.pdf")

        # plotter.colormap = "hot"

        plotter.vmin = 0.0

        plotter.set_title("Prepared images")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sky values ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        sky_path = fs.join(self.plot_preparation_path, "sky")
        if not fs.is_directory(sky_path): fs.create_directory(sky_path)

        # Loop over the different images
        for label in self.sorted_labels:

            not_nan = Mask.is_nan(self.images[label]).inverse()

            # Create the distribution from the image pixel values
            distribution = Distribution.from_values("Pixel value", self.images[label][not_nan].flatten() + self.sky_values[label])

            # Create an array of all the pixels used for estimating the sky
            #notnan = np.logical_not(np.isnan(self.apertures))
            #print(self.apertures.dtype)
            notnan = Mask.is_nan(self.apertures[label]).inverse()
            sky_values = self.apertures[label][notnan]

            # Create the distribution of pixel values used for the sky estimation
            sky_distribution = Distribution.from_values("Pixel value", sky_values)

            # Add the distributions
            plotter.add_distribution(distribution, label)
            plotter.add_distribution(sky_distribution, label)

            # Plot seperately
            distr_plotter = DistributionPlotter()
            distr_plotter.add_distribution(distribution, "image")
            distr_plotter.add_distribution(sky_distribution, "sky")
            distr_plotter.run(fs.join(sky_path, label + ".pdf"))

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "sky_distribution.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_errors(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the errors ...")

        # Plot histograms of the absolute error values
        #self.plot_error_histograms_absolute()

        # Plot histograms of the relative error values
        #self.plot_error_histograms_relative()

        # Plot the relative errors of each pixel
        self.plot_errors_pixels()

    # -----------------------------------------------------------------

    def plot_error_histograms_absolute(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absolute error values in a histogram for each prepared image compared to the histogram of the actual image values ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        absolute_errors_path = fs.join(self.plot_preparation_path, "absolute_errors")
        if not fs.is_directory(absolute_errors_path): fs.create_directory(absolute_errors_path)

        # Loop over the different images
        for label in self.sorted_labels:

            not_nan = Mask.is_nan(self.images[label]).inverse()

            # Create the distribution from the image pixel values
            distribution = Distribution.from_values("Pixel value", self.images[label][not_nan].flatten())

            not_nan = Mask.is_nan(self.errors[label]).inverse()

            # Create the distribution from the error values
            error_distribution = Distribution.from_values("Error", self.errors[label][not_nan].flatten())

            # Add an entry to the distribution grid plotter
            plotter.add_distribution(distribution, label)
            plotter.add_distribution(error_distribution, label)

            # Plot seperately
            distr_plotter = DistributionPlotter()
            distr_plotter.add_distribution(distribution, "image")
            distr_plotter.add_distribution(error_distribution, "absolute errors")
            distr_plotter.run(fs.join(absolute_errors_path, label + ".pdf"))

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "absolute_errors.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_error_histograms_relative(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the relative error values in a histogram for each prepared image ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        relative_errors_path = fs.join(self.plot_preparation_path, "relative_errors")
        if not fs.is_directory(relative_errors_path): fs.create_directory(relative_errors_path)

        # Loop over the different images
        for label in self.sorted_labels:

            # Calculate the relative errors
            rel_errors = self.errors[label] / self.images[label]

            # Create a distribution from the relative errors
            rel_error_distribution = Distribution.from_values("Relative error", rel_errors)

            # Add the distribution to the plotter
            plotter.add_distribution(rel_error_distribution, label)

            # Plot seperately
            rel_error_distribution.plot(title="relative errors", path=fs.join(relative_errors_path, label + ".pdf"))

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "relative_errors.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_errors_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the relative error for each pixel in each prepared image ...")

        # Create the ErrorPlotter instance
        plotter = ErrorPlotter()

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "errors_pixels.png")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_masks_and_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images with the sources masks and sky annuli overlayed ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in self.sorted_labels: plotter.add_image(self.images[label], label, mask=self.sources_masks[label], region=self.annuli[label])

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "preparation_masks_annuli.pdf")

        plotter.vmin = 0.0

        plotter.set_title("Prepared images with sources masks and sky annuli")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the aperture frames with the sky annuli overlayed ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in self.sorted_labels: plotter.add_image(self.apertures[label], label, region=self.annuli[label])

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "preparation_apertures.pdf")

        plotter.vmin = 0.0

        plotter.set_title("Aperture frames with sky annuli")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_labels(self):

        """
        This function ...
        :return:
        """

        sorted_labels = sorted(self.images.keys(), key=lambda key: self.images[key].filter.pivotwavelength())
        return sorted_labels

# -----------------------------------------------------------------
