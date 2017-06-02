#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.truncation Contains the Truncator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...magic.region.list import SkyRegionList, PixelRegionList
from .component import TruncationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.dist_ellipse import distance_ellipse
from ...core.basics.range import RealRange
from ...core.basics.map import Map
from ...magic.core.mask import intersection
from ...core.remote.remote import Remote

# -----------------------------------------------------------------

class Truncator(TruncationComponent):
    
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
        super(Truncator, self).__init__(*args, **kwargs)

        # --- Attributes ---

        # The statistics for each image
        self.statistics = dict()

        # The frames and error maps
        self.frames = None
        self.errormaps = None
        self.masks = None

        # The sky ellipses
        self.sky_ellipses = dict()

        # Truncation ellipse
        self.ellipses = defaultdict(dict)

        # Paths
        self.paths = dict()

        # The remote host (if needed)
        self.remote = None

        # The remote cache path
        self.remote_truncation_path = None

        # Cache paths
        self.cache_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the data
        self.load_data()

        # 3. Create directories
        self.create_directories()

        # 4. Find the best radius for the truncation
        self.calculate_statistics()

        # 5. Create the ellipses
        self.create_ellipses()

        # 6. Writing
        self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Truncator, self).setup(**kwargs)

        # Setup the remote
        self.remote = Remote(host_id=self.environment.cache_host_id)

        # Create the cache directory
        self.remote_truncation_path = fs.join(self.remote.home_directory, self.galaxy_name + "_truncation")
        if self.config.cache:
            if not self.remote.is_directory(self.remote_truncation_path): self.remote.create_directory(self.remote_truncation_path)

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return len(self.frames)

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the data ...")

        # Get the frames
        self.frames = self.dataset.get_framelist()

        # Get the error maps
        self.errormaps = self.dataset.get_errormaplist()

        # Loop over all prepared images, get the images
        self.masks = dict()
        for name in self.dataset.names:

            # Get the mask
            mask_names = ["padded", "bad"]
            mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Set the mask
            if mask is None: continue
            self.masks[name] = mask

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating a directory for each image ...")

        # Loop over the image
        for name in self.frames.names:

            # Create directory
            path = fs.create_directory_in(self.truncation_path, name)

            # Set path
            self.paths[name] = path

            # Create remote directory
            if self.config.cache:

                # Create directory
                remote_path = self.remote.create_directory_in(self.remote_truncation_path, name)

                # Set path
                self.cache_paths[name] = remote_path

    # -----------------------------------------------------------------

    def calculate_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the statistics as a function of radius ...")

        # Get the angle
        center = self.disk_ellipse.center  # in sky coordinates
        semimajor = self.disk_ellipse.semimajor
        semiminor = self.disk_ellipse.semiminor
        angle = self.disk_ellipse.angle

        # Determine the ratio of semimajor and semiminor
        ratio = semiminor / semimajor

        # Loop over all prepared images
        for name in self.frames.names:

            # Get the image
            frame = self.dataset.get_frame(name)

            # Get the mask
            mask_names = ["padded", "bad"]
            mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Convert center to pixel coordinates
            center_pix = center.to_pixel(frame.wcs)

            # Create distance ellipse frame
            distance_frame = distance_ellipse(frame.shape, center_pix, ratio, angle)

            radius_list = []
            signal_to_noise_list = []
            nmasked_list = []

            # Loop over the radii
            min_distance = np.min(distance_frame)
            max_distance = np.max(distance_frame)
            step = (max_distance - min_distance) / float(self.config.nbins)

            # Set the first range
            radius_range = RealRange(min_distance, min_distance + step)

            # Loop, shifting ranges of radius
            while True:

                # Check the range
                if radius_range.min > max_distance: break

                # Get the average radius
                radius_center = radius_range.center

                above_min_mask = distance_frame >= radius_range.min
                below_max_mask = distance_frame < radius_range.max

                # Make a mask of the pixels corresponding to the current radius range
                #range_mask = radius_range.min <= distance_frame < radius_range.max

                range_mask = intersection(above_min_mask, below_max_mask)

                # Calculate the mean signal to noise in the pixels
                signal_to_noises = self.frames[name][range_mask] / self.errormaps[name][range_mask]

                # Calcalute the mean signal to noise
                signal_to_noise = np.mean(signal_to_noises)

                # Make a mask of all the pixels below the center radius
                below_mask = distance_frame < radius_center

                # Calculate the number of masked pixels
                nmasked = np.sum(mask[below_mask])
                ntotal = np.sum(below_mask)
                rel_nmasked = nmasked / ntotal

                # Add point
                radius_list.append(radius_center)
                signal_to_noise_list.append(signal_to_noise)
                nmasked_list.append(rel_nmasked)

                # Shift the range
                radius_range += step

            # Set the statistics for this image
            statistics = Map()
            statistics.radii = radius_list
            statistics.snr = signal_to_noise_list
            statistics.nmasked = nmasked_list
            self.statistics[name] = statistics

    # -----------------------------------------------------------------

    @lazyproperty
    def factors(self):

        """
        This function ..
        :return: 
        """

        return self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

    # -----------------------------------------------------------------

    def create_ellipses(self):

        """
        This function ....
        :return: 
        """

        # Inform the user
        log.info("Creating ellipses ...")

        # Loop over the different scale factors
        for factor in self.factors:

            # Get the scaled ellipse
            sky_ellipse = self.disk_ellipse * factor

            # Add the sky ellipse
            self.sky_ellipses[factor] = sky_ellipse

            # Loop over the frames
            for name in self.frames.names:

                # Convert to pixel ellipse
                pixel_ellipse = sky_ellipse.to_pixel(self.frames[name].wcs)

                # Add the ellipse
                self.ellipses[name][factor] = pixel_ellipse

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the truncation ellipse
        self.write_ellipses()

        # Write the truncated images
        self.write_images()

    # -----------------------------------------------------------------

    def write_ellipses(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ellipses ...")

        # Write sky ellipses
        self.write_sky_ellipses()

        # Write image ellipses
        self.write_image_ellipses()

    # -----------------------------------------------------------------

    def write_sky_ellipses(self):

        """
        Thisf ucntion ...
        :return: 
        """

        # Inform the user
        log.info("Writing the truncation ellipse region ...")

        # Determine the path to the region file
        path = fs.join(self.truncation_path, "ellipses.reg")

        # Create the region list
        regions = SkyRegionList()

        # Loop over the ellipses
        for factor in self.sky_ellipses:

            # Add ellipse
            ellipse = self.sky_ellipses[factor]
            ellipse.meta["text"] = str(factor)
            regions.append(ellipse)

        # Write
        regions.saveto(path)

    # -----------------------------------------------------------------

    def write_image_ellipses(self):

        """
        Thisf ucntion ...
        :return: 
        """

        # Inform the user
        log.info("Writing the image ellipses ...")

        # Loop over the images
        for name in self.ellipses:

            # Get the path
            path = fs.join(self.paths[name], "ellipses.reg")

            # Create region list
            regions = PixelRegionList()

            # Loop over the ellipses
            for factor in self.ellipses[name]:

                # Add ellipse
                ellipse = self.ellipses[name][factor]
                ellipse.meta["text"] = str(factor)
                regions.append(ellipse)

            # Write
            regions.saveto(path)

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the truncated images ...")

        # Loop over the the images
        index = 0
        for name in self.frames.names:

            # Debugging
            index += 1
            progress = float(index) / float(self.nframes)
            log.debug("Writing truncated images for the " + name + " image (" + str(index) + " of " + str(self.nframes) + ") ...")

            # Loop over the factors
            for factor in self.factors:

                # Determine the local path
                filename = str(factor) + ".fits"
                path = fs.join(self.paths[name], filename)

                # Determine the remote path
                remote_path = fs.join(self.cache_paths[name], filename)

                # Already existing
                if fs.is_file(path):

                    # Debugging
                    log.debug("Truncated " + name + " image with factor " + str(factor) + "is already present: not creating it again")

                    # Cache if requested
                    if self.config.cache:

                        # Upload
                        self.remote.upload_file_to(path, self.cache_paths[name], remove=True)

                        # Debugging
                        log.debug("Truncated " + name + " image with factor " + str(factor) + " has been cached to '" + remote_path + "'")

                # Already present remotely
                elif self.remote.is_file(remote_path):

                    # Debugging
                    log.debug("Truncated " + name + " image with factor " + str(factor) + " is already present and cached on remote host '" + self.remote.host_id + "'")

                # Not yet present, create truncated image (and cache)
                else:

                    # Debugging
                    log.debug("Creating the truncated " + name + " image with factor " + str(factor) + "...")

                    # Get the pixel ellipse
                    ellipse = self.ellipses[name][factor]

                    # Convert into mask
                    mask = ellipse.to_mask(self.frames[name].xsize, self.frames[name].ysize)

                    # Truncate the frame
                    frame = self.frames[name]
                    frame[mask] = 0.0

                    # Save
                    frame.saveto(path)

                    # Cache
                    if self.config.cache:

                        # Upload, and remove local file
                        remote_path = self.remote.upload_file_to(path, self.cache_paths[name], remove=True)

                        # Debugging
                        log.debug("Truncated " + name + " image with factor " + str(factor) + " has been cached to '" + remote_path + "'")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the curves
        self.plot_snr()

        # Plot nmasked pixels
        self.plot_nmasked()

    # -----------------------------------------------------------------

    def plot_snr(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting the snr curves ...")

        # Loop over the frame names
        for name in self.statistics:

            # Get x and y
            radii = self.statistics[name].radii
            snr = self.statistics[name].snr

            # Create plot
            plt.figure()
            plt.plot(radii, snr)

            # Add vertical lines
            for factor in self.ellipses[name]:
                radius = self.ellipses[name][factor].major
                plt.axvline(x=radius)

            # Determine the path
            path = fs.join(self.paths[name], "snr.pdf")

            # Save the figure
            plt.savefig(path)
            plt.close()

    # -----------------------------------------------------------------

    def plot_nmasked(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting the nmasked pixels curves ...")

        # Loop over the frame nems
        for name in self.statistics:

            # Get x and y
            radii = self.statistics[name].radii
            nmasked = self.statistics[name].nmasked

            # Create plot
            plt.figure()
            plt.plot(radii, nmasked)

            # Add vertical lines
            for factor in self.ellipses[name]:
                radius = self.ellipses[name][factor].major
                plt.axvline(x=radius)

            # Determine the path
            path = fs.join(self.paths[name], "nmasked.pdf")

            # Save the figure
            plt.savefig(path)
            plt.close()

# -----------------------------------------------------------------
