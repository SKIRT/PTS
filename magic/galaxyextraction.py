#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
import matplotlib.pylab as plt

# Import astronomical modules
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# Import the relevant AstroMagic classes and modules
from .basics import Mask, Region, Position, Extent
from .core import Source
from .sky import Galaxy
from .tools import catalogs, regions

# Import the relevant PTS classes and modules
from ..core.basics import Configurable

# -----------------------------------------------------------------

class GalaxyExtractor(Configurable):

    """
    This class
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(GalaxyExtractor, self).__init__(config)

        ## Attributes

        # Initialize an empty list for the galaxies
        self.galaxies = []

        # Initialize an empty list to contain the manual sources
        self.manual_sources = []

        # Set the frame to None
        self.frame = None

        # Set the mask to None
        self.mask = None

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup(frame)

        # 2. Find and remove the galaxies
        self.find_and_remove_galaxies()

        # 3. If a manual region was specified, remove the corresponding galaxies
        if self.config.manual_region is not None: self.set_and_remove_manual()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame):

        """
        This function ...
        """

        # Call the setup function of the base class
        super(GalaxyExtractor, self).setup()

        # Make a local reference to the passed frame
        self.frame = frame

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask(np.zeros_like(self.frame))

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the galaxy extractor")

        # Clear the list of galaxies
        self.galaxies = []

        # Clear the list of manual sources
        self.manual_sources = []

        # Clear the frame and the mask
        self.frame = None
        self.mask = None

    # -----------------------------------------------------------------

    def find_and_remove_galaxies(self):

        """
        This function ...
        :return:
        """

        # Get list of galaxies
        self.fetch_galaxies()

        # Set special galaxies
        if self.config.special_region is not None: self.set_special()

        # Find the sources
        self.find_sources()

        # Find apertures
        if self.config.find_apertures: self.find_apertures()

        # If requested, remove
        if self.config.remove: self.remove_galaxies()

    # -----------------------------------------------------------------

    def set_and_remove_manual(self):

        """
        This function ...
        :return:
        """

        # Set manual galaxies
        self.set_manual()

        # If requested, remove the manually specified galaxies
        if self.config.remove_manual: self.remove_manual()

    # -----------------------------------------------------------------

    @property
    def positions(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the galaxies
        for galaxy in self.galaxies:

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(galaxy.pixel_position(self.frame.wcs))

        # Return the list
        return positions

    # -----------------------------------------------------------------

    @property
    def principal(self):

        """
        This function ...
        :return:
        """

        # Loop over the list of galaxies
        for galaxy in self.galaxies:

            # Check if it is the principal galaxy; if so, return it
            if galaxy.principal: return galaxy

    # -----------------------------------------------------------------

    @property
    def companions(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the companion galaxies
        companions = []

        # Loop over the list of galaxies
        for galaxy in self.galaxies:

            # Check if it is a companion galaxy; if so, add it to the list
            if galaxy.companion: companions.append(galaxy)

        # Return the list of companion galaxies
        return companions

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Looking for sources near the galaxy positions")

        # Loop over all galaxies in the list
        for galaxy in self.galaxies:

            # If this sky object should be ignored, skip it
            if galaxy.ignore: continue

            # If requested, use the galaxy extents obtained from the catalog to create the source
            if self.config.detection.use_d25 and galaxy.has_extent:

                outer_factor = self.config.detection.background_outer_factor
                expansion_factor = self.config.detection.d25_expansion_factor
                galaxy.source_from_parameters(self.frame, outer_factor, expansion_factor)

            else:

                # Find a source
                try: galaxy.find_source(self.frame, self.config.detection)
                except Exception:
                    #import traceback
                    self.log.error("Error when finding source")
                    #print(type(e))
                    #print(e)
                    #traceback.print_exc()
                    if self.config.plot_track_record_if_exception:
                        if galaxy.has_track_record: galaxy.track_record.plot()
                        else: self.log.warning("Track record is not enabled")

            # If a source was not found for the principal or companion galaxies, force it
            outer_factor = self.config.detection.background_outer_factor
            if galaxy.principal and not galaxy.has_source: galaxy.source_from_parameters(self.frame, outer_factor)
            elif galaxy.companion and not galaxy.has_source and galaxy.has_extent: galaxy.source_from_parameters(self.frame, outer_factor)

        # Inform the user
        self.log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.galaxies), self.have_source/len(self.galaxies)*100.0))

    # -----------------------------------------------------------------

    def fetch_galaxies(self):

        """
        This function ...
        :param image:
        :param config:
        :return:
        """

        # Inform the user
        self.log.info("Fetching galaxy positions from an online catalog")

        # Get the range of right ascension and declination of the image
        center, ra_span, dec_span = self.frame.coordinate_range()

        # Find galaxies in the box defined by the center and RA/DEC ranges
        for galaxy_name in catalogs.galaxies_in_box(center, ra_span, dec_span):

            # Create a Galaxy object and add it to the list
            self.galaxies.append(Galaxy(galaxy_name))

        # Define a function that returns the length of the major axis of the galaxy
        def major_axis(galaxy):

            if galaxy.major is None: return 0.0
            else: return galaxy.major

        # Indicate which galaxy is the principal galaxy
        principal_galaxy = max(self.galaxies, key=major_axis)
        principal_galaxy.principal = True

        # Loop over the galaxies, and enable track record if requested
        if self.config.track_record:

            for galaxy in self.galaxies: galaxy.enable_track_record()

        # Loop over the galaxies, check if they are companion galaxies
        for galaxy in self.galaxies:

            # Skip the principal galaxy
            if galaxy.principal: continue

            if galaxy.type == "HII": galaxy.companion = True
            if galaxy.name[:-1].lower() == principal_galaxy.name[:-1].lower(): galaxy.companion = True

            if galaxy.companion:
                principal_galaxy.companions.append(galaxy.name)
                galaxy.parent = principal_galaxy.name

        # Inform the user
        self.log.debug("Number of galaxies: " + str(len(self.galaxies)))
        self.log.debug(self.principal.name + " is the principal galaxy in the frame")
        self.log.debug("The following galaxies are its companions: " + str(self.principal.companions))

    # -----------------------------------------------------------------

    def find_apertures(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        self.log.info("Constructing elliptical apertures regions to encompass the detected galaxies")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If this galaxy should be ignored, skip it
            if galaxy.ignore: continue

            # If the galaxy does not have a source, continue
            if galaxy.has_source: galaxy.find_aperture(self.frame, self.config.apertures)

    # -----------------------------------------------------------------

    def remove_galaxies(self):

        """
        This function ...
        :param image:
        :param galaxies:
        :param config:
        :return:
        """

        # Inform the user
        self.log.info("Removing the galaxies from the frame (except for the principal galaxy and its companions)")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If this galaxy should be ignored, skip it
            if galaxy.ignore: continue

            # Remove the galaxy from the frame
            if not galaxy.principal and not galaxy.companion: galaxy.remove(self.frame, self.mask, self.config.removal)

        # Add the principal and companion galaxies to the mask
        self.mask += self.principal_mask
        self.mask += self.companion_mask

    # -----------------------------------------------------------------

    def set_special(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        self.log.info("Setting special region from " + self.config.special_region)

        # Load the region and create a mask from it
        region = Region.from_file(self.config.special_region, self.frame.wcs)
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Get the position of this object in pixel coordinates
            position = galaxy.pixel_position(self.frame.wcs)

            # Set special if position is covered by the mask
            if special_mask.masks(position): galaxy.special = True

    # -----------------------------------------------------------------

    def set_ignore(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        self.log.info("Setting region to ignore for subtraction from " + self.config.ignore_region)

        # Load the region and create a mask from it
        region = Region.from_file(self.config.ignore_region, self.frame.wcs)
        ignore_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Get the position of this object in pixel coordinates
            position = galaxy.pixel_position(self.frame.wcs)

            # Ignore if position is covered by the mask
            if ignore_mask.masks(position): galaxy.ignore = True

    # -----------------------------------------------------------------

    def set_manual(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Setting region for manual galaxy extraction from " + self.config.manual_region)

        # Load the region and create a mask from it
        region = Region.from_file(self.config.manual_region, self.frame.wcs)

        # Loop over the shapes in the region
        for shape in region:

            # Get the center and radius of the shape (can be a circle or an ellipse)
            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            # Create a source
            source = Source(self.frame, Position(x_center, y_center), Extent(x_radius, y_radius), Angle(0.0, u.deg), self.config.manual.background_outer_factor)

            # Add the source to the list of manual sources
            self.manual_sources.append(source)

    # -----------------------------------------------------------------

    def remove_manual(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Removing manually specified galaxies from the frame")

        # Loop over each item in the list of manual sources
        for source in self.manual_sources:

            # Estimate the background for the source
            source.estimate_background(self.config.manual.interpolation_method, self.config.manual.sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

    # -----------------------------------------------------------------

    @property
    def principal_mask(self):

        """
        This function ...
        :return:
        """

        # Create a new mask with the dimensions of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Add the principal galaxy's mask to the total mask
        mask[self.principal.source.cutout.y_slice, self.principal.source.cutout.x_slice] = self.principal.source.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @property
    def companion_mask(self):

        """
        This function ...
        :return:
        """

        # Create a new mask with the dimension of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Loop over all companion galaxies
        for galaxy in self.companions:

            # Check if the galaxy has a source and add its mask to the total mask
            if galaxy.has_source: mask[galaxy.source.cutout.y_slice, galaxy.source.cutout.x_slice] = galaxy.source.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @property
    def region(self):

        """
        This function ...
        :param image:
        :param stars:
        :param config:
        :return:
        """

        # TODO: improve this function

        ra_list = []
        dec_list = []
        height_list = []
        width_list = []
        angle_list = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            ra_list.append(galaxy.position.ra.value)
            dec_list.append(galaxy.position.dec.value)

            if galaxy.major is not None:

                width = galaxy.major.to("arcsec").value

                if galaxy.minor is not None: height = galaxy.minor.to("arcsec").value
                else: height = width

            else:

                width = self.config.region.default_radius * self.frame.pixelscale.value
                height = self.config.region.default_radius * self.frame.pixelscale.value

            angle = galaxy.pa.degree if galaxy.pa is not None else 0.0

            height_list.append(height)
            width_list.append(width)
            angle_list.append(angle)

        # Create a region and return it
        return regions.ellipses(ra_list, dec_list, height_list, width_list, angle_list)

    # -----------------------------------------------------------------

    def write_region(self):

        """
        This function ...
        :param frame:
        :return:
        """

        path = self.config.writing.region_path
        annotation = self.config.writing.region_annotation

        # Inform the user
        self.log.info("Writing galaxy region to " + path)

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Get the center in pixel coordinates
            x_center, y_center = galaxy.position.to_pixel(self.frame.wcs, origin=0)

            # Set the angle
            angle = galaxy.pa.degree if galaxy.pa is not None else 0.0

            if galaxy.major is None:

                color = "red"

                x_radius = self.config.region.default_radius
                y_radius = self.config.region.default_radius

            elif galaxy.minor is None or galaxy.pa is None:

                color = "green"

                x_radius = 0.5 * galaxy.major.to("arcsec").value / self.frame.pixelscale.value
                y_radius = x_radius

            else:

                color = "green"

                x_radius = 0.5 * galaxy.major.to("arcsec").value / self.frame.pixelscale.value
                y_radius = 0.5 * galaxy.minor.to("arcsec").value / self.frame.pixelscale.value

            # Annotation
            if annotation == "name": text = "text = {" + galaxy.name + "}"
            elif annotation == "redshift": text = "text = {" + str(galaxy.redshift) + "}"
            elif annotation == "type": text = "text = {" + str(galaxy.type) + "}"
            elif annotation is None: text = ""
            else: raise ValueError("Invalid option for annotation")

            color_suffix = " # color = " + color
            point_suffix = " # point = x " + text

            # Add point for the center
            print("image;point({},{})".format(x_center, y_center) + point_suffix, file=f)
            print("image;ellipse({},{},{},{},{})".format(x_center, y_center, x_radius, y_radius, angle) + color_suffix, file=f)

            # Add aperture
            if galaxy.has_aperture:

                ap_x_center, ap_y_center = galaxy.aperture.positions[0]
                major = galaxy.aperture.a
                minor = galaxy.aperture.b
                angle = galaxy.aperture.theta / math.pi * 180

                aperture_suffix = " # color = white"

                print("image;ellipse({},{},{},{},{})".format(ap_x_center, ap_y_center, major, minor, angle) + aperture_suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing table to " + self.config.writing.table_path)

        # Write the table to file
        self.table.write(self.config.writing.table_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Writing masked frame to " + self.config.writing.masked_frame_path)

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(self.config.writing.masked_frame_path)

    # -----------------------------------------------------------------

    def write_result(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing resulting frame to " + self.config.writing.result_path)

        # Write out the resulting frame
        self.frame.save(self.config.writing.result_path)

    # -----------------------------------------------------------------

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for galaxy in self.galaxies: count += galaxy.has_source
        return count

    # -----------------------------------------------------------------

    @property
    def table(self):

        """
        This function ...
        :return:
        """

        names = ["name", "type", "distance", ""]
        values = [[] for i in range(len(names))]

        #for galaxy in self.galaxies:

            #for index, name in enumerate(names):

                #value = getattr(galaxy, name)
                #values[index].append(value)

        names = []
        types = []
        distances = []
        ascensions = []
        declinations = []
        majors = []
        d25s = []
        position_angles = []
        inclinations = []
        principal_flags = []
        companion_flags = []
        sources = []

        # Loop over all stars
        for galaxy in self.galaxies:

            names.append(galaxy.name)
            types.append(galaxy.type)
            if galaxy.distance is not None: distances.append(galaxy.distance.to(u.Mpc).value)
            else: distances.append(None)
            ascensions.append(galaxy.position.ra.value)
            declinations.append(galaxy.position.dec.value)
            if galaxy.major is not None: majors.append(galaxy.major.to(u.arcmin).value)
            else: majors.append(None)
            if galaxy.d25 is not None: d25s.append(galaxy.d25.to(u.arcmin).value)
            else: d25s.append(None)
            if galaxy.pa is not None: position_angles.append(galaxy.pa.degree)
            else: position_angles.append(None)
            if galaxy.inclination is not None: inclinations.append(galaxy.inclination.degree)
            else: inclinations.append(None)
            principal_flags.append(galaxy.principal)
            companion_flags.append(galaxy.companion)
            sources.append(galaxy.has_source)

        # Create the table
        table = Table([names, types, distances, ascensions, declinations, majors, d25s, position_angles, inclinations, principal_flags, companion_flags, sources],
                     names=('Name', 'Type', 'Distance', 'RA', 'DEC', 'Major axis length', 'D25', 'Position angle', 'Inclination', 'Principal', 'Companion', 'Source'),
                     meta={'name': 'galaxies'})

        # Set units for columns
        table['Distance'].unit = u.Mpc
        table['D25'].unit = u.arcmin
        table['Major axis length'].unit = u.arcmin
        table['Position angle'].unit = u.deg
        table['Inclination'].unit = u.deg

        # Return the table
        return table

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        x_centers = []
        y_centers = []
        apertures = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            x_center, y_center = galaxy.position.to_pixel(self.frame.wcs)
            x_centers.append(x_center)
            y_centers.append(y_center)

            # If the galaxy does not have a source, continue
            if galaxy.has_aperture: apertures.append(galaxy.aperture)

        # Initialize the plot
        #norm = ImageNormalize(stretch=SqrtStretch())
        norm = ImageNormalize(stretch=LogStretch())
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Determine the maximum value in the box and the mimimum value for plotting
        vmax = np.max(self.frame)
        vmin = np.min(self.frame) if vmax <= 0 else 0.0

        # Plot the frame and the segments mask
        ax1.imshow(self.frame, origin='lower', interpolation='nearest', norm=norm, vmin=vmin, vmax=vmax)
        ax2.imshow(self.mask, origin='lower', cmap='jet')

        # Set axes limits
        plt.xlim(0, self.frame.xsize-1)
        plt.ylim(0, self.frame.ysize-1)

        # Plot the apertures
        for aperture in apertures:

            aperture.plot(color='white', lw=1.5, alpha=0.5, ax=ax1)
            aperture.plot(color='white', lw=1.5, alpha=1.0, ax=ax2)

        # Plot centers of the galaxies
        plt.plot(x_centers, y_centers, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)

        # Show the plot
        plt.show()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # If requested, write out a table with the galaxy properties
        if self.config.write_table: self.write_table()

        # If requested, write out the galaxy region
        if self.config.write_region: self.write_region()

        # If requested, write out the frame where the galaxies are masked
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the result
        if self.config.write_result: self.write_result()

# -----------------------------------------------------------------