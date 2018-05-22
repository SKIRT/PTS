#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.inspector Contains the DataInspector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

# Import astronomical modules
from astropy.io import fits

# Import the relevant PTS classes and modules
from .component import DataComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.dataset import DataSet
from ...core.basics.plot import pretty_colors
from ...core.tools import archive
from ...magic.core.fits import contains_pc_and_cd, contains_pc_and_cdelt, remove_cd_keywords, remove_cdelt_keywords
from ...core.tools import time
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

essential_origins = ["Herschel", "GALEX", "Spitzer"]

# -----------------------------------------------------------------

class DataInspector(DataComponent):
    
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
        super(DataInspector, self).__init__(*args, **kwargs)

        # Paths
        self.paths = None
        self.error_paths = None

        # The initial dataset
        self.set = DataSet()

        # The frame list
        self.frames = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inspect the directories
        self.inspect_directories()

        # 2. Get the image paths
        self.get_paths()

        # Create the dataset
        self.create_dataset()

        # Check headers
        self.check_headers()

        # Get the frames
        self.get_frames()

        # 3. Writing
        self.write()

        # Show
        self.show()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DataInspector, self).setup(**kwargs)

        # Set the output path
        directory_name = time.unique_name(self.command_name())
        self.config.output = fs.create_directory_in(self.inspect_path, directory_name)

    # -----------------------------------------------------------------

    def inspect_directories(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting directories ...")

        # Inspect origins
        self.inspect_origins()

        # Inspect contents
        self.inspect_contents()

    # -----------------------------------------------------------------

    def inspect_origins(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the origins ...")

        # Get the origins
        origin_names = fs.directories_in_path(self.data_images_path, returns="name")

        # Loop over the essential origins
        for essential_origin in essential_origins:
            if essential_origin not in origin_names: raise RuntimeError("No " + essential_origin + " images are present. Run get_images again to fix this.")

    # -----------------------------------------------------------------

    def inspect_contents(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the contents ...")

        # Loop over the directories in data/image
        for origin_path, origin_name in fs.directories_in_path(self.data_images_path, returns=["path", "name"]):

            # Check whether not empty
            if fs.is_empty(origin_path):
                log.warning("No images for '" + origin_name + "'")
                if origin_name in essential_origins:
                    raise RuntimeError("No " + origin_name + " images are present. Run get_images again to fix this.")

            # Check whether compressed files present
            compressed_paths = fs.compressed_files_in_path(origin_path)
            if len(compressed_paths) > 0:
                log.warning("Compressed files present for '" + origin_name + "': decompressing ...")
                archive.decompress_files(compressed_paths, remove=True)

    # -----------------------------------------------------------------

    def get_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for images and error frames ...")

        # Get paths
        self.paths, self.error_paths = self.get_data_image_and_error_paths()

    # -----------------------------------------------------------------

    def create_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the initial dataset ...")

        # Loop over the image paths
        for prep_name in self.paths:

            # Add entry to the dataset
            self.set.add_path(prep_name, self.paths[prep_name])

    # -----------------------------------------------------------------

    def check_headers(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the image headers ...")

        # Loop over the image paths
        for name in self.paths:

            # Load the header
            header = self.set.get_header(name)

            cleaned = False

            # Check keywords
            if contains_pc_and_cd(header):
                log.warning("Header of '" + name + "' image contains PC and CD keywords: removing CD information ...")
                remove_cd_keywords(header)
                cleaned = True

            # Replace the header
            if contains_pc_and_cdelt(header):
                log.warning("Header of '" + name + "' image contains PC and CDELT keywords: removing CDELT information ...")
                remove_cdelt_keywords(header)
                cleaned = True

            # Set the new header
            if cleaned:

                # Open
                hdulist = fits.open(self.set.paths[name])

                # Get header
                hdu = hdulist[0]

                hdu.header = header

                # Remove original file
                fs.remove_file(self.set.paths[name])

                # Replace
                hdu.writeto(self.set.paths[name])

    # -----------------------------------------------------------------

    @lazyproperty
    def names_per_origin(self):

        """
        This function ...
        :return: 
        """

        origins = defaultdict(list)

        for name in self.set.names:

            path = self.set.get_frame_path(name)
            origin = fs.name(fs.directory_of(path))

            # Add the name
            origins[origin].append(name)

        # Return the origins
        return origins

    # -----------------------------------------------------------------

    def get_frames(self):

        """
        This function ...
        :return: 
        """

        self.frames = self.set.get_framelist()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the uer
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This funciton ...
        :return: 
        """

        # Inform the user
        log.info("Showing ...")

        # Show the coordinate ranges
        self.show_coordinate_ranges()

        # Show the pixelscales
        self.show_pixelscales()

    # -----------------------------------------------------------------

    def show_coordinate_ranges(self):

        """
        This function ...
        :return: 
        """

        print("")

        for origin in self.names_per_origin:

            print(origin)
            print("")

            for name in self.names_per_origin[origin]:

                print(" - " + name + ":")

                ra_range = self.frames[name].wcs.ra_range
                dec_range = self.frames[name].wcs.dec_range

                center, ra_span, dec_span = self.frames[name].wcs.coordinate_range

                #bounding_box = self.frames[name].wcs.bounding_box

                pixelscale = self.frames[name].wcs.pixelscale

                min_ra = ra_range.min
                max_ra = ra_range.max
                min_dec = dec_range.min
                max_dec = dec_range.max

                #print("  * RA range: " + str(ra_range))
                #print("  * DEC range: " + str(dec_range))
                #print("  * coordinate range: " + str(coordinate_range))
                #print("  * pixelscale: " + str(pixelscale))

                # THIS WILL NOT BE EQUAL: RA SPAN CANNOT BE ADDED OR SUBTRACTED TO RA COORDINATES (PROJECTION)
                min_ra_range = center.ra - 0.5 * ra_span
                max_ra_range = center.ra + 0.5 * ra_span
                min_dec_range = center.dec - 0.5 * dec_span
                max_dec_range = center.dec + 0.5 * dec_span

                min_ra_deg = min_ra.to("deg").value
                max_ra_deg = max_ra.to("deg").value
                min_dec_deg = min_dec.to("deg").value
                max_dec_deg = max_dec.to("deg").value

                min_ra_range_deg = min_ra_range.to("deg").value
                max_ra_range_deg = max_ra_range.to("deg").value
                min_dec_range_deg = min_dec_range.to("deg").value
                max_dec_range_deg = max_dec_range.to("deg").value

                print(np.isclose(min_ra_deg, min_ra_range_deg))
                print(np.isclose(max_ra_deg, max_ra_range_deg))
                print(np.isclose(min_dec_deg, min_dec_range_deg))
                print(np.isclose(max_dec_deg, max_dec_range_deg))

                print("")

            print("")

    # -----------------------------------------------------------------

    def show_pixelscales(self):

        """
        This function ...
        :return: 
        """

        print("")

        for origin in self.names_per_origin:

            print(origin)
            print("")

            for name in self.names_per_origin[origin]:

                print(" - " + name + ": " + str(self.frames[name].average_pixelscale))

            print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the coordinate boxes
        self.plot_coordinate_boxes()

    # -----------------------------------------------------------------

    def plot_coordinate_boxes(self):

        """
        This function ...
        :return: 
        """

        plt.figure()

        min_ra = None
        max_ra = None
        min_dec = None
        max_dec = None

        ax = plt.gca()

        itercolors = iter(pretty_colors + pretty_colors)

        #for_legend = OrderedDict()

        # Loop over the coordinate boxes
        for name in self.set.names:

            # Debugging
            log.debug("Plotting coordinate box of " + name + " image ...")

            # Get coordinate box
            box = self.set.get_coordinate_box(name)

            x_center = box.center.ra.to("deg").value
            y_center = box.center.dec.to("deg").value
            width = 2. * box.radius.ra.to("deg").value
            height = 2. * box.radius.dec.to("deg").value

            #print(x_center, y_center)
            #print(width, height)

            color = itercolors.next()

            x_lower_left = x_center - 0.5 * width
            y_lower_left = y_center - 0.5 * height

            # x, y is lower left
            ph = ax.add_patch(Rectangle((x_lower_left, y_lower_left), width, height, fill=None, alpha=1, edgecolor=color, lw=1, label=name))

            if min_ra is None or x_center - 0.5 * width < min_ra: min_ra = x_center - 0.5 * width
            if max_ra is None or x_center + 0.5 * width > max_ra: max_ra = x_center + 0.5 * width
            if min_dec is None or y_center - 0.5 * height < min_dec: min_dec = y_center - 0.5 * height
            if max_dec is None or y_center + 0.5 * height > max_dec: max_dec = y_center + 0.5 * height

            # Add patch handle
            #for_legend[name] = ph

        #for name in for_legend:
        #    plt.legend()

        #plt.legend()

        # Put a legend below current axis
        legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=6, prop={'size': 12})

        #print(min_ra, max_ra)
        #print(min_dec, max_dec)

        # Add overlapping box
        overlap_box, data = self.set.get_overlap_data()

        print("Min RA frame: " + data.ra.min)
        print("Max RA frame: " + data.ra.max)
        print("Min DEC frame: " + data.dec.min)
        print("Max DEC frame: " + data.dec.max)

        x_center = overlap_box.center.ra.to("deg").value
        y_center = overlap_box.center.dec.to("deg").value
        width = 2. * overlap_box.radius.ra.to("deg").value
        height = 2. * overlap_box.radius.dec.to("deg").value

        x_lower_left = x_center - 0.5 * width
        y_lower_left = y_center - 0.5 * height

        # x, y is lower left
        ax.add_patch(Rectangle((x_lower_left, y_lower_left), width, height, fill=None, alpha=1))

        # Add galaxy ellipse
        ellipse = self.galaxy_ellipse
        galaxy_x_center = ellipse.center.ra.to("deg").value
        galaxy_y_center = ellipse.center.dec.to("deg").value
        galaxy_width = 2.0 * ellipse.radius.ra.to("deg").value
        galaxy_height = 2.0 * ellipse.radius.dec.to("deg").value
        angle = ellipse.angle.to("deg").value

        # xy, width, height, angle=0.0
        # angle in degrees anti-clockwise
        ax.add_patch(Ellipse((galaxy_x_center, galaxy_y_center), galaxy_width, galaxy_height, angle, fill=None, alpha=1))

        ra_width = max_ra - min_ra
        dec_width = max_dec - min_dec
        ra_center = 0.5 * (min_ra + max_ra)
        dec_center = 0.5 * (min_dec + max_dec)

        # Add 10 % on each side
        min_ra = ra_center - 1.1 * 0.5 * ra_width
        max_ra = ra_center + 1.1 * 0.5 * ra_width
        min_dec = dec_center - 1.1 * 0.5 * dec_width
        max_dec = dec_center + 1.1 * 0.5 * dec_width

        # Set limits
        plt.xlim((min_ra, max_ra))
        plt.ylim((min_dec, max_dec))

        #plt.show()

        # Determine the path
        path = self.output_path_file("coordinate_boxes.pdf")

        # Save the figure and close
        plt.savefig(path)
        plt.close()

# -----------------------------------------------------------------
