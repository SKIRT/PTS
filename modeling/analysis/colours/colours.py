#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.colours.colours Contains the ColourAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import ColourAnalysisComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame
from ....magic.plot.imagegrid import ResidualImageGridPlotter
from ....core.basics.distribution import Distribution
from ....magic.core.mask import Mask
from ....core.tools import sequences

# -----------------------------------------------------------------

# Names that identify the interesting wavelengths/filters
colour_names = ["70/100", "100/160", "160/250", "250/350", "350/500"]
keys = ["70", "100", "160", "250", "350", "500"]
ids = {"Pacs blue": "70", "Pacs green": "100", "Pacs red": "160", "SPIRE PSW": "250", "SPIRE PMW": "350", "SPIRE PLW": "500"}
filter_names = {"70": "Pacs blue", "100": "Pacs green", "160": "Pacs red", "250": "SPIRE PSW", "350": "SPIRE PMW", "500": "SPIRE PLW"}

# -----------------------------------------------------------------

class ColourAnalyser(ColourAnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ColourAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The dictionaries that contain the observed and simulated far-infrared images
        self.observed = dict()
        self.errors = dict()
        self.simulated = dict()

        # The dictionaries that contain the observed and simulated colour maps
        self.observed_colours = dict()
        self.simulated_colours = dict()

        # The residual colour maps
        self.residual_colours = dict()

        # The residual pixel value distributions
        self.residual_distributions = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Load the images
        self.load_images()

        # 2. Rebin the images to the same pixel grid
        self.rebin()

        # 3. Calculate the colour maps
        self.calculate_colours()

        # 4. Calculate the residual colour maps
        self.calculate_residuals()

        # Create distributions
        self.create_distributions()

        # 5. Writing
        self.write()

        # 6. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        Thisn fuction ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColourAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # 1. Load the observed images
        self.load_observed_images()

        # 2. Load the simulated images
        self.load_simulated_images()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over the appropriate observed images
        filter_names_keys = [filter_names[key] for key in keys]
        image_names_for_filters = self.dataset.get_names_for_filters(filter_names_keys)
        for key, image_name in zip(keys, image_names_for_filters):

            # Get the corresponding filter
            #fltr = BroadBandFilter(filter_names[key])

            # Check if not None (in database)
            if image_name is None:
                log.warning("No observed " + key + " image was found")
                continue

            # Get the observed image for this filter
            frame = self.dataset.get_frame(image_name)

            # Add the frame to the dictionary
            self.observed[key] = frame

            # Get the error map
            errors = self.dataset.get_errormap(image_name)

            # Add to the dictionary
            self.errors[key] = errors

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over the appropriate simulated images
        for key in keys:

            # Check if in observed
            if key not in self.observed:
                log.warning("An observed " + key + " image was not found: not loading the simulated version")
                continue

            # Determine the path to the image
            filter_name = filter_names[key]
            image_name = "earth__" + filter_name
            path = fs.join(self.analysis_run.misc_path, image_name + ".fits")

            # Check whether the image is present
            if not fs.is_file(path):
                log.warning("Simulated " + key + " micron image is not present")
                continue

            # Load the image frame
            frame = Frame.from_file(path)

            # Add the frame to the dictionary
            self.simulated[key] = frame

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning all images to the pixel grid of the observed image with the lowest resolution ...")

        # Set the target coordinate system to None initially
        target_wcs = None

        # Get the WCS of the observed image with the largest pixelscale
        for key in self.observed:
            if target_wcs is None or self.observed[key].wcs.average_pixelscale > target_wcs.average_pixelscale: target_wcs = self.observed[key].wcs

        # Debugging
        log.debug("The target coordinate system has a pixelscale of " + str(target_wcs.average_pixelscale) + " and a shape of (nx=" + str(target_wcs.naxis1) + ", ny=" + str(target_wcs.naxis2) + ")")

        # Rebin all observed images to the target coordinate system
        for key in self.observed:

            # Debugging
            log.debug("Checking observed images ...")

            # Check whether the coordinate systems match
            if self.observed[key].wcs == target_wcs:

                # Debugging
                log.debug("The coordinate system of the observed " + key + " micron image matches the target coordinate system")

            else:

                # Debugging
                log.debug("The coordinate system of the observed " + key + " micron image does not match the target coordinate system: rebinning ...")

                # Rebin the image to the specific coordinate system
                self.observed[key] = self.observed[key].rebinned(target_wcs)

                # Rebin the error map
                self.errors[key] = self.errors[key].rebinned(target_wcs)

        # Rebin all simulated images to the target coordinate system
        for key in self.simulated:

            # Check whether the coordinate systems match
            if self.simulated[key].wcs == target_wcs:

                # Debugging
                log.debug("The coordinate system of the simulated " + key + " micron image matches the target coordinate system")

            else:

                # Debugging
                log.debug("The coordinate system of the simulated " + key + " micron image does not match the target coordinate system: rebinning ...")

                # Rebin the image to the standardized coordinate system
                self.simulated[key] = self.simulated[key].rebinned(target_wcs)

    # -----------------------------------------------------------------

    def calculate_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the colour maps ...")

        # Loop over the different colours
        for colour_name in colour_names:

            # Get the wavelengths to create the colour map
            wavelength1, wavelength2 = colour_name.split("/")

            # Check whether the corresponding observed and simulated images exist
            if not wavelength1 in self.observed:
                log.warning("Observed '" + filter_names[wavelength1] + "' image not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength1 in self.simulated:
                log.warning("Simulated '" + filter_names[wavelength1] + "' image not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.observed:
                log.warning("Observed '" + filter_names[wavelength2] + "' image not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.simulated:
                log.warning("Simulated '" + filter_names[wavelength2] + "' image not present, " + colour_name + " can not be calculated")
                continue

            # Debugging
            log.debug("Calculating the observed and simulated " + colour_name + " colour maps ...")

            # Calculate the colour maps
            observed_colour = Frame(np.log10(self.observed[wavelength1]/self.observed[wavelength2]), wcs=self.observed[wavelength1].wcs)
            simulated_colour = Frame(np.log10(self.simulated[wavelength1]/self.simulated[wavelength2]), wcs=self.observed[wavelength1].wcs)

            # Replace infs
            observed_colour.replace_infs(0.0)
            simulated_colour.replace_infs(0.0)

            observed = self.observed[wavelength1]

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            #significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)
            significance_mask = self.create_clip_mask(wavelength1, wavelength2)

            # MASK TRUNCATION ELLIPSE
            observed_colour[truncation_mask] = 0.0
            simulated_colour[truncation_mask] = 0.0

            # MASK SIGNIFICANE
            observed_colour[significance_mask] = 0.0
            simulated_colour[significance_mask] = 0.0

            # Add the colour maps to the appropriate dictionary
            self.observed_colours[colour_name] = observed_colour
            self.simulated_colours[colour_name] = simulated_colour

    # -----------------------------------------------------------------

    def create_clip_mask(self, key1, key2, fuzzy=False):

        """
        This function ...
        :param key1:
        :param key2:
        :param fuzzy:
        :return:
        """

        # Debugging
        log.debug("Making clip mask for " + key1 + " and " + key2 + " images ...")

        frame1 = self.observed[key1]
        errors1 = self.errors[key1]
        filter_name1 = filter_names[key1]
        level1 = self.get_significance_level(filter_name1)

        frame2 = self.observed[key2]
        errors2 = self.errors[key2]
        filter_name2 = filter_names[key2]
        level2 = self.get_significance_level(filter_name2)

        # Create significance maps
        significance1 = frame1 / errors1
        significance2 = frame2 / errors2

        # Create masks
        if fuzzy:
            mask1 = create_fuzzy_mask_for_level(significance1, level1, fuzziness=fuzziness, offset=fuzziness_offset)
            mask2 = create_fuzzy_mask_for_level(significance2, level2, fuzziness=fuzziness, offset=fuzziness_offset)
        else:
            mask1 = create_mask_for_level(significance1, level1)
            mask2 = create_mask_for_level(significance2, level2)

        from ....magic.core.mask import intersection
        from ....magic.core.alpha import product

        # Make list
        masks = [mask1, mask2]

        # Create intersection mask
        if fuzzy: mask = product(*masks)
        else: mask = intersection(*masks)

        # Only keep largest patch
        mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

        # Fill holes
        mask.fill_holes()

        # Invert FOR NORMAL MASKS: WE HAVE TO SET PIXELS TO ZERO THAT ARE NOT ON THE MASK
        if not fuzzy: mask.invert()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating residuals between the observed and simulated colour maps ...")

        # Loop over the different colours
        for colour_name in self.observed_colours:

            # Debugging
            log.debug("Calculating residuals between the observed and simulated " + colour_name + " colour maps ...")

            # Calculate the residual
            residual = (self.observed_colours[colour_name] - self.simulated_colours[colour_name])/self.observed_colours[colour_name]

            # Replace infs and nans
            residual.replace_infs(0.0)
            residual.replace_nans(0.0)

            # Add the residual map to the dictionary
            self.residual_colours[colour_name] = residual

    # -----------------------------------------------------------------

    # def calculate_residual_distributions(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Calculating distributions of residual pixel values ...")
    #
    #     # Loop over the different colours
    #     for colour_name in self.observed_colours:
    #
    #         # Debugging
    #         log.debug("Calculating the distribution for the pixels of the " + colour_name + " residual map ...")
    #
    #         # Get an 1D array of the valid pixel values
    #         pixel_values = None
    #
    #         # Create the distribution
    #         distribution = Distribution.from_values(pixel_values)
    #
    #         # Debugging
    #         #log.debug("Median " + colour_name + " residual: " + str(np.nanmedian(np.abs(residual))))
    #         #log.debug("Standard deviation of " + colour_name + " residual: " + str(np.nanstd(residual)))
    #
    #         # Add the distribution to the dictionary
    #         self.residual_distributions[colour_name] = distribution

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the distributions ...")

        # Loop over the residual colour maps
        for colour_name in self.residual_colours:

            # Get the map
            residuals = self.residual_colours[colour_name]

            # Get the values within the truncation ellipse
            values = residuals.values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
            indices = np.argwhere(values < -20)
            values = np.delete(values, indices)
            indices = np.argwhere(values > 20)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create residuals distribution for the '" + colour_name + "' colour")
                continue

            # Create distribution
            distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

            # Add the distribution
            self.residual_distributions[colour_name] = distribution

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the observed colour maps
        self.write_observed_colours()

        # Write the simulated colour maps
        self.write_simulated_colours()

        # Write the residual colour maps
        self.write_residual_colours()

        # Write the distributions of the residual pixel values
        self.write_residual_distributions()

    # -----------------------------------------------------------------

    def write_observed_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the observed colour maps ...")

        # Loop over the observed colour maps
        for colour_name in self.observed_colours:

            # Determine the path
            path = fs.join(self.colours_observed_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.observed_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_simulated_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulated colour maps ...")

        # Loop over the simulated colour maps
        for colour_name in self.simulated_colours:

            # Determine the path
            path = fs.join(self.colours_simulated_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.simulated_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_residual_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual colour maps ...")

        # Loop over the residual colour maps
        for colour_name in self.residual_colours:

            # Determine the path
            path = fs.join(self.colours_residuals_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.residual_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_residual_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual colour pixel distributions ...")

        # Loop over the residual pixel distributions
        for colour_name in self.residual_distributions:

            # Determine the path
            path = fs.join(self.colours_distributions_path, colour_name.replace("/", "-") + ".dat")

            # Save the distribution data
            self.residual_distributions[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot a grid with the observed, simulated and residual colour maps
        #self.plot_image_grid()

        # Plot the residual maps
        #self.plot_residuals()

        # Plot the distributions of the pixel values of the residual maps
        self.plot_residual_distributions()

    # -----------------------------------------------------------------

    def plot_image_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual colour maps ...")

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Colours")

        # Add the rows
        for colour_name in self.observed_colours:

            observed_colour = self.observed_colours[colour_name]
            simulated_colour = self.simulated_colours[colour_name]

            plotter.add_row(observed_colour, simulated_colour, colour_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(self.truncation_box)

        # Determine the path to the plot file
        path = fs.join(self.analysis_run.colours_path, "colours.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residual maps ...")

        # Loop over the colours
        for colour_name in self.residual_colours:

            # Determine the path
            path = fs.join(self.colours_residuals_path, colour_name.replace("/", "-") + ".pdf")

            # Save as PNG
            residuals = self.residual_colours[colour_name]
            vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
                                              interval=self.config.interval,
                                              scale=self.config.scale, alpha=self.config.alpha_method,
                                              peak_alpha=self.config.peak_alpha)

    # -----------------------------------------------------------------

    def plot_residual_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the pixel value distributions for the residual maps ...")

        # Loop over the residual pixel distributions
        for colour_name in self.residual_distributions:

            # Determine the path
            path = fs.join(self.colours_distributions_path, colour_name.replace("/", "-") + ".pdf")

            # Save the distribution data
            self.residual_distributions[colour_name].plot(title="Distribution of residual pixel values for " + colour_name + " colour", path=path)

# -----------------------------------------------------------------

def create_mask_for_level(significance, level):

    """
    This function ...
    :param significance:
    :param level:
    :return:
    """

    # Create the mask
    mask = Mask(significance > level, wcs=significance.wcs)

    # Only keep largest patch
    #mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

    # Fill holes
    #mask.fill_holes()

    # Return the mask
    return mask

# -----------------------------------------------------------------

def create_fuzzy_mask_for_level(significance, level, fuzziness, offset=1.0):

    """
    This function ...
    :param significance:
    :param level:
    :param fuzziness:
    :param offset:
    :return:
    """

    # Determine the maximum significance
    max_significance = np.nanmax(significance)

    # Debugging
    log.debug("Maximal significance: " + str(max_significance))

    # Construct value range
    lower_relative = 1. - fuzziness  # example: 1. - 0.1
    upper_relative = 1. + fuzziness

    # ADAPT IF THE RANGE IS TOO BROAD
    if level * upper_relative > max_significance - offset:
        log.warning("Changing the upper relative sigma level for the fuzzy edge from " + str(upper_relative) + " to 1")
        upper_relative = 1.

    value_range = RealRange.around(level, lower_relative, upper_relative)

    # Debugging
    log.debug("Sigma level range for fuzzy edge: " + str(value_range))

    # Check maximum of the range
    if value_range.max + offset > max_significance: raise ValueError("This should not happen")

    # Create the mask
    mask = AlphaMask.between(significance, value_range)

    # Only keep largest patch
    #mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

    # Fill holes
    #mask.fill_holes()

    # Return the mask
    return mask

# -----------------------------------------------------------------
