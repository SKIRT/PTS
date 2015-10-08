#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import numpy as np
from config import Config
import inspect
import matplotlib.pylab as plt

# Import Astromagic modules
from ..core import regions
from ..core.masks import Mask
from ..core.galaxy import Galaxy
from ..core.vector import Position

# Import astromagic modules
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.ned import Ned
import astroquery.exceptions
from astropy.coordinates import Angle

# *****************************************************************

class GalaxyExtractor(object):

    """
    This class
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "galaxyextractor.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

        # Initialize an empty list for the galaxies
        self.galaxies = []

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the list of galaxies
        self.galaxies = []

    # *****************************************************************

    def run(self, frame):

        """
        This function ...
        """

        # Get list of galaxies
        self.fetch_galaxies(frame)

        # Find the sources
        self.find_sources(frame)

        # Find apertures
        if self.config.find_apertures: self.find_apertures()

        # If requested, remove
        if self.config.remove: self.remove_galaxies(frame)

    # *****************************************************************

    def fetch_galaxies(self, frame):

        """
        This function ...
        :param image:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy positions from an online catalog")

        # Get the range of right ascension and declination of the image
        center, ra_span, dec_span = frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["galaxies", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=["VII/237"])
        table = result[0]

        # Loop over all galaxies in the table
        for entry in table:

            # Get the PGC number
            pgc_id = entry["PGC"]

            # Obtain more information about this galaxy
            try:
                ned_result = Ned.query_object("PGC " + str(pgc_id))
                ned_entry = ned_result[0]

                #notes = Ned.get_table("PGC " + str(pgc_id), table='objectnotes')

                # Get the NGC number
                name = ned_entry["Object Name"]

                # Get the redshift
                redshift = ned_entry["Redshift"]

                # Get the type (G=galaxy, HII ...)
                type = ned_entry["Type"]

            except astroquery.exceptions.RemoteServiceError:

                # Set attributes
                name = "PGC " + str(pgc_id)
                redshift = None
                type = None

            satellite = False
            if type == "HII": satellite = True

            # Get the right ascension and the declination
            center = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # Get the names given to this galaxy
            names = entry["ANames"].split() if entry["ANames"] else None

            # Get the size of the galaxy
            ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
            diameter = np.power(10.0, entry["logD25"])*0.1*u.arcmin if entry["logD25"] else None

            # Get the size of major and minor axes
            major = diameter
            minor = diameter / ratio if diameter is not None and ratio is not None else None

            # Get the position angle of the galaxy
            pa = Angle(entry["PA"]-90, u.deg) if entry["PA"] else None

            # Add a new galaxy to the list
            self.galaxies.append(Galaxy(pgc_id=pgc_id, name=name, type=type, redshift=redshift, center=center, names=names, major=major, minor=minor, pa=pa, satellite=satellite))

        # Define a function that returns the length of the major axis of the galaxy
        def major_axis(galaxy):

            if galaxy.major is None: return 0.0
            else: return galaxy.major

        # Indicate which galaxy is the principal galaxy
        galaxy = max(self.galaxies, key=major_axis)
        galaxy.principal = True

        # Inform the user
        log.debug("Number of galaxies: " + str(len(self.galaxies)))
        log.debug(galaxy.name + " is the principal galaxy in the frame")

    # *****************************************************************

    def find_sources(self, frame):

        """
        This function ...
        :param image:
        :return:
        """

        # Inform the user
        log.info("Looking for sources near the galaxy positions")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Find a source
            galaxy.find_source(frame, self.config.detection)

        # Inform the user
        log.debug("Success ratio: {0:.2f}%".format(self.have_source/len(self.galaxies)*100.0))

    # *****************************************************************

    def find_apertures(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If the galaxy does not have a source, continue
            if galaxy.has_source: galaxy.find_aperture(sigma_level=self.config.apertures.sigma_level)

    # *****************************************************************

    def remove_galaxies(self, frame):

        """
        This function ...
        :param image:
        :param galaxies:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Removing the galaxies from the frame (except for the principal galaxy and its satellites)")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Remove the galaxy from the frame
            if not galaxy.principal and not galaxy.satellite: galaxy.remove(frame, self.config.removal)

    # *****************************************************************

    def create_region(self, frame):

        """
        This function ...
        :param image:
        :param stars:
        :param config:
        :return:
        """

        ra_list = []
        dec_list = []
        height_list = []
        width_list = []
        angle_list = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            ra_list.append(galaxy.center.ra.value)
            dec_list.append(galaxy.center.dec.value)

            if galaxy.major is not None:

                width = galaxy.major.to("arcsec").value

                if galaxy.minor is not None: height = galaxy.minor.to("arcsec").value
                else: height = width

            else:

                width = self.config.region.default_radius * frame.pixelscale.value
                height = self.config.region.default_radius * frame.pixelscale.value

            angle = galaxy.pa.degrees if galaxy.pa is not None else 0.0

            height_list.append(height)
            width_list.append(width)
            angle_list.append(angle)

        # Create a region and return it
        return regions.ellipses(ra_list, dec_list, height_list, width_list, angle_list)

    # *****************************************************************

    def create_mask(self, frame):

        """
        This function ...
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(frame))

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If no source was found for the galaxy, skip it
            if not galaxy.has_source: continue

            # Add this galaxy to the mask
            if self.config.mask.use_aperture and galaxy.has_aperture:

                galaxy_mask_frame = Mask.from_aperture(frame.xsize, frame.ysize, galaxy.aperture)
                galaxy_mask = galaxy_mask_frame[galaxy.source.cutout.y_min:galaxy.source.cutout.y_max, galaxy.source.cutout.x_min:galaxy.source.cutout.x_max]

            else: galaxy_mask = galaxy.source.mask

            # Add this galaxy to the total mask
            mask[galaxy.source.cutout.y_min:galaxy.source.cutout.y_max, galaxy.source.cutout.x_min:galaxy.source.cutout.x_max] += galaxy_mask

        # Expand the mask
        #if self.config.mask.dilate: mask.dilate(self.config.mask.connectivity, self.config.mask.iterations)

        # Return the mask
        return mask

    # *****************************************************************

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for galaxy in self.galaxies: count += galaxy.has_source
        return count

    # *****************************************************************

    def position_list(self, frame):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the galaxy positions
        position_list = []

        # Loop over the galaxies
        for galaxy in self.galaxies:

            # Calculate the pixel coordinate in the frame and add it to the list
            x_center, y_center = galaxy.center.to_pixel(frame.wcs)
            position_list.append(Position(x=x_center, y=y_center))

        # Return the list
        return position_list

    # *****************************************************************

    def table(self):

        """
        This function ...
        :return:
        """

        names = []
        types = []
        redshifts = []
        ascensions = []
        declinations = []
        majors = []
        position_angles = []
        principal_flags = []
        satellite_flags = []
        sources = []
        #other_names = []

        # Loop over all stars
        for galaxy in self.galaxies:

            names.append(galaxy.name)
            types.append(galaxy.type)
            redshifts.append(galaxy.redshift)
            ascensions.append(galaxy.center.ra.value)
            declinations.append(galaxy.center.dec.value)
            if galaxy.major is not None: majors.append(galaxy.major)
            else: majors.append(None)
            position_angles.append(galaxy.pa)
            principal_flags.append(galaxy.principal)
            satellite_flags.append(galaxy.satellite)
            sources.append(galaxy.has_source)
            #other_names.append(galaxy.names)

        print(len(names), len(types), len(redshifts), len(ascensions), len(declinations), len(majors), len(position_angles), len(principal_flags), len(satellite_flags), len(sources))

        print(redshifts[3], type(redshifts[3]))

        # Create and return the table
        return Table([names, types, redshifts, ascensions, declinations, majors, position_angles, principal_flags, satellite_flags, sources],
                     names=('Name', 'Type', 'Redshift', 'RA', 'DEC', 'Major axis', 'Position angle', 'Principal', 'Satellite', 'Source'),
                     meta={'name': 'galaxies'})

    # *****************************************************************

    def plot(self, frame):

        """
        This function ...
        :return:
        """

        from astropy.visualization import SqrtStretch, LogStretch
        from astropy.visualization.mpl_normalize import ImageNormalize

        x_centers = []
        y_centers = []
        apertures = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            x_center, y_center = galaxy.center.to_pixel(frame.wcs)
            x_centers.append(x_center)
            y_centers.append(y_center)

            # If the galaxy does not have a source, continue
            if galaxy.has_aperture: apertures.append(galaxy.aperture)

        # Initialize the plot
        #norm = ImageNormalize(stretch=SqrtStretch())
        norm = ImageNormalize(stretch=LogStretch())
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Determine the maximum value in the box and the mimimum value for plotting
        vmax = np.max(frame)
        vmin = np.min(frame) if vmax <= 0 else 0.0

        # Plot the frame and the segments mask
        ax1.imshow(frame, origin='lower', interpolation='nearest', norm=norm, vmin=vmin, vmax=vmax)
        ax2.imshow(self.create_mask(frame), origin='lower', cmap='jet')

        # Set axes limits
        plt.xlim(0, frame.xsize-1)
        plt.ylim(0, frame.ysize-1)

        # Plot the apertures
        for aperture in apertures:

            aperture.plot(color='white', lw=1.5, alpha=0.5, ax=ax1)
            aperture.plot(color='white', lw=1.5, alpha=1.0, ax=ax2)

        # Plot centers of the galaxies
        plt.plot(x_centers, y_centers, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)

        # Show the plot
        plt.show()

# *****************************************************************
