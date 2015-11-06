#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import os.path
import numpy as np
import inspect
import matplotlib.pylab as plt

# Import Astromagic modules
from .objectextraction import ObjectExtractor
from ..tools import catalogs
from ..core import regions
from ..core.galaxy import Galaxy
from ..tools import configuration
from ..core.masks import Mask

# Import astronomical modules
from astropy import log
import astropy.units as u
from astropy.table import Table
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# *****************************************************************

class GalaxyExtractor(ObjectExtractor):

    """
    This class
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "galaxyextractor.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        # Call the constructor of the base class
        super(GalaxyExtractor, self).__init__()

    # *****************************************************************

    def run(self, frame):

        """
        This function ...
        """

        # Call the setup function
        self.setup(frame)

        # Get list of galaxies
        self.fetch_galaxies()

        # Set special galaxies
        if self.config.special_region is not None: self.set_special()

        # Find the sources
        self.find_sources()

        # If a source was not found for the principal galaxy, force it
        outer_factor = self.config.detection.background_outer_factor
        if not self.principal.has_source: self.principal.source_from_parameters(frame, outer_factor)

        # Find apertures
        if self.config.find_apertures: self.find_apertures()

        # If requested, remove
        if self.config.remove: self.remove_galaxies()

        # If a manual region was specified, remove the corresponding galaxies
        if self.config.manual_region is not None: self.remove_manual()

        # If requested, write a table with the galaxy properties
        if self.config.save_table: self.save_table()

        # If requested, save the galaxy region
        if self.config.save_region: self.save_region()

        # If requested, save the frame where the galaxies are masked
        if self.config.save_masked_frame: self.save_masked_frame()

        # If requested, save the result
        if self.config.save_result: self.save_result()

    # *****************************************************************

    @property
    def principal(self):

        """
        This function ...
        :return:
        """

        # Loop over the list of galaxies
        for galaxy in self.objects:

            # Check if it is the principal galaxy; if so, return it
            if galaxy.principal: return galaxy

    # *****************************************************************

    def fetch_galaxies(self):

        """
        This function ...
        :param image:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy positions from an online catalog")

        # Get the range of right ascension and declination of the image
        center, ra_span, dec_span = self.frame.coordinate_range()

        # Find galaxies in the box defined by the center and RA/DEC ranges
        for galaxy_name in catalogs.galaxies_in_box(center, ra_span, dec_span):

            # Create a Galaxy object and add it to the list
            self.objects.append(Galaxy(galaxy_name))

        # Define a function that returns the length of the major axis of the galaxy
        def major_axis(galaxy):

            if galaxy.major is None: return 0.0
            else: return galaxy.major

        # Indicate which galaxy is the principal galaxy
        principal_galaxy = max(self.objects, key=major_axis)
        principal_galaxy.principal = True

        # Loop over the galaxies, and enable track record if requested
        if self.config.track_record:

            for galaxy in self.objects: galaxy.enable_track_record()

        # Loop over the galaxies, check if they are companion galaxies
        for galaxy in self.objects:

            # Skip the principal galaxy
            if galaxy.principal: continue

            if galaxy.type == "HII": galaxy.companion = True
            if galaxy.name[:-1].lower() == principal_galaxy.name[:-1].lower(): galaxy.companion = True

            if galaxy.companion:
                principal_galaxy.companions.append(galaxy.name)
                galaxy.parent = principal_galaxy.name

        # Inform the user
        log.debug("Number of galaxies: " + str(len(self.objects)))
        log.debug(self.principal.name + " is the principal galaxy in the frame")
        log.debug("The following galaxies are its companions: " + str(self.principal.companions))

    # *****************************************************************

    def find_apertures(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Constructing elliptical apertures regions to encompass the detected galaxies")

        # Loop over all galaxies
        for galaxy in self.objects:

            # If this galaxy should be ignored, skip it
            if galaxy.ignore: continue

            # If the galaxy does not have a source, continue
            if galaxy.has_source: galaxy.find_aperture(sigma_level=self.config.apertures.sigma_level)

    # *****************************************************************

    def remove_galaxies(self):

        """
        This function ...
        :param image:
        :param galaxies:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Removing the galaxies from the frame (except for the principal galaxy and its companions)")

        # Loop over all galaxies
        for galaxy in self.objects:

            # If this galaxy should be ignored, skip it
            if galaxy.ignore: continue

            # Remove the galaxy from the frame
            if not galaxy.principal and not galaxy.companion: galaxy.remove(self.frame, self.config.removal)

    # *****************************************************************

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
        for galaxy in self.objects:

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

    # *****************************************************************

    def write_region(self, path, annotation="name"):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Writing galaxy region to " + path)

        # Create a file
        f = open(path,'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Loop over all galaxies
        for galaxy in self.objects:

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

    # *****************************************************************

    @property
    def mask(self):

        """
        This function ...

        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Loop over all galaxies
        for galaxy in self.objects:

            # If a source was found for this galaxy
            if galaxy.has_source:

                # Check whether the extent (as defined by the D25 isophotal radius) of the galaxy should be used for the mask
                if self.config.mask.use_d25 and galaxy.has_extent:

                    # Obtain the ellipse parameters (default radius should not be defined, since we only get here
                    # after we checked galaxy.has_extent
                    center, radius, angle = galaxy.ellipse_parameters(self.frame.wcs, self.frame.pixelscale, None)

                    # Create a mask from the ellipse parameters
                    d25_mask = Mask.from_ellipse(self.frame.xsize, self.frame.ysize, center, radius*self.config.mask.d25_expansion_factor, angle)

                    # Add the d25 mask to the total mask
                    mask += d25_mask

                # Else, check whether the aperture created from the center segment should be used
                elif self.config.mask.use_aperture and galaxy.has_aperture:

                    # Create a mask from the aperture of the object (expand if specified under self.config.aperture_mask)
                    object_mask_frame = Mask.from_aperture(self.frame.xsize, self.frame.ysize, galaxy.aperture, expansion_factor=self.config.aperture_removal.expansion_factor)

                    # Now, we don't limit setting the mask within the source's cutout, because we expanded the apertures to perhaps a size larger than this cutout,
                    # so just add the object_mask_frame to the total frame
                    mask += object_mask_frame

                # Else, use the source mask (the detected center segment)
                else:

                    # Add the source mask of this galaxy to the total mask
                    mask[galaxy.source.cutout.y_min:galaxy.source.cutout.y_max, galaxy.source.cutout.x_min:galaxy.source.cutout.x_max] += galaxy.source.mask

            # If a source was not detected for this galaxy but the inclusion of undetected galaxies in the mask is enabled
            elif self.config.mask.include_undetected:

                # If the galaxy has no extents (D25 is not defined), skip it
                if not galaxy.has_extent: continue

                # Obtain the ellipse parameters (default radius should not be defined, since we only get here
                # after we checked galaxy.has_extent
                center, radius, angle = galaxy.ellipse_parameters(self.frame.wcs, self.frame.pixelscale, None)

                # Create a mask from the ellipse parameters
                d25_mask = Mask.from_ellipse(self.frame.xsize, self.frame.ysize, center, radius, angle)

                # Add the d25 mask to the total mask
                mask += d25_mask

        # Return the mask
        return mask

    # *****************************************************************

    @property
    def table(self):

        """
        This function ...
        :return:
        """

        names = ["name", "type", "distance", ""]
        values = [[] for i in range(len(names))]

        #for galaxy in self.objects:

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
        for galaxy in self.objects:

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

    # *****************************************************************

    def plot(self):

        """
        This function ...
        :return:
        """

        x_centers = []
        y_centers = []
        apertures = []

        # Loop over all galaxies
        for galaxy in self.objects:

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

# *****************************************************************
