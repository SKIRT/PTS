#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import math
import numpy as np
import os.path
import inspect
import matplotlib.pyplot as plt
from config import Config

# Import Astromagic modules
from .objectextraction import ObjectExtractor
from ..core.star import Star
from ..tools import statistics
from ..core.vector import Position
from ..core.regions import Region
from ..core.source import Source

# Import astronomical modules
import aplpy
import astropy.io.fits as pyfits
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from astroquery.vizier import Vizier

# *****************************************************************

class StarExtractor(ObjectExtractor):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # If no configuration is given, use the default configuration
        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "starextractor.cfg")
            config = Config(file(config_path))

        # Call the constructor of the base class
        super(StarExtractor, self).__init__(config)

    # *****************************************************************

    def run(self, frame, galaxyextractor=None, regionpath=None):

        """
        This function ...
        """

        # Create a list of stars based on online catalogs
        self.fetch_stars(frame, galaxyextractor)

        # Set special stars
        if self.config.special_region is not None: self.set_special(frame)

        # Set ignore stars
        if self.config.ignore_region is not None: self.set_ignore(frame)

        # For each star, find a corresponding source in the image
        self.find_sources(frame)

        # Fit analytical models to the stars
        self.fit_stars(frame)

        #self.write_region(frame, regionpath[:-4]+"_0.reg", annotation="has_background")

        # If requested, remove the stars
        if self.config.remove: self.remove_stars(frame, galaxyextractor)

        #self.write_region(frame, regionpath[:-4]+"_1.reg", annotation="has_background")

        # If requested, remove saturation in the image
        if self.config.remove_saturation: self.remove_saturation(frame, galaxyextractor)

        #self.write_region(frame, regionpath[:-4]+"_2.reg", annotation="has_background")

        # If requested, find apertures
        if self.config.find_apertures: self.find_apertures()

    # *****************************************************************

    def fetch_stars(self, frame, galaxyextractor=None):

        """
        This function ...
        """

        # Inform the user
        log.info("Fetching star positions from an online catalog")

        # Get the range of right ascension and declination of this image
        center, ra_span, dec_span = frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=self.config.fetching.catalogs)
        table = result[0]

        # If a galaxyextractor is passed, get the positions of the galaxies in the frame
        if galaxyextractor is not None: galaxies = galaxyextractor.position_list(frame)
        else: galaxies = []

        # Loop over all stars in the table
        norms = []
        for entry in table:

            # Get the PGC number
            ucac_id = entry["UCAC4"]

            # Get the right ascension and the declination
            position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # If this star does not lie within the frame, skip it
            if not frame.contains(position): continue

            # Get the mean error on the right ascension and declination
            position_error = entry["ePos"]*u.mas
            ra_error = entry["e_RAJ2000"]*u.mas
            dec_error = entry["e_DEJ2000"]*u.mas

            # Loop over all galaxy positions in the list
            for galaxy in galaxies:

                # Calculate the distance between the star's position and the galaxy's center
                x_center, y_center = position.to_pixel(frame.wcs)
                difference = galaxy - Position(x=x_center, y=y_center)

                norms.append(difference.norm)

                if difference.norm <= self.config.fetching.min_difference_from_galaxy:

                    # Remove the position of this galaxy from the list
                    galaxies.remove(galaxy)
                    break

            else:

                # Get the magnitude in different bands
                k_mag = entry["Kmag"]*u.mag
                b_mag = entry["Bmag"]*u.mag
                v_mag = entry["Vmag"]*u.mag
                r_mag = entry["rmag"]*u.mag
                i_mag = entry["imag"]*u.mag

                # Create a star object
                star = Star(ucac_id=ucac_id, position=position, position_error=position_error, ra_error=ra_error,
                            dec_error=dec_error, k_mag=k_mag, b_mag=b_mag, v_mag=v_mag, r_mag=r_mag, i_mag=i_mag)

                # If requested, enable track record
                if self.config.track_record: star.enable_track_record()

                # Add the star to the list of stars
                self.objects.append(star)

        # Inform the user
        if galaxyextractor is not None: log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(norm) for norm in sorted(norms)[:10]))

        # Inform the user
        log.debug("Number of stars: " + str(len(self.objects)))

    # *****************************************************************

    def fit_stars(self, frame):

        """
        This function ...
        """

        # Inform the user
        log.info("Fitting analytical profiles to the sources")

        # Loop over all stars in the list
        for star in self.objects:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Check if the star has a source (has been detected)
            if not star.has_source and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                center, radius, angle = star.ellipse_parameters(frame.wcs, frame.pixelscale, self.config.fitting.initial_radius)

                # Create a source object
                source = Source(frame, center, radius, angle, self.config.fitting.background_outer_factor)

            else: source = None

            # Find a source
            if star.has_source or source is not None: star.fit_model(self.config.fitting, source)

        # Inform the user
        log.debug("Found a model for {0} out of {1} stars with source ({2:.2f}%)".format(self.have_model, self.have_source, self.have_model/self.have_source*100.0))

    # *****************************************************************

    def remove_stars(self, frame, galaxyextractor=None):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing the stars from the frame")

        # Get a list of the fwhm of all the stars that were fitted
        fwhm_list = self.fwhm_list

        # Determine the fwhm for the stars that were not fitted
        if self.config.removal.no_source_fwhm == "max": default_fwhm = max(fwhm_list)
        elif self.config.removal.no_source_fwhm == "mean": default_fwhm = np.mean(fwhm_list)
        elif self.config.removal.no_source_fwhm == "median": default_fwhm = np.median(fwhm_list)
        else: raise ValueError("Unkown measure for determining the fwhm of stars without detected source")

        # Inform the user
        log.debug("Default FWHM used when star could not be fitted: {0:.2f} pixels".format(default_fwhm))

        # Loop over all stars in the list
        for star in self.objects:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Determine whether we want the background to be sigma-clipped when interpolating over the source
            if galaxyextractor.principal.contains(star.pixel_position(frame.wcs)) and self.config.removal.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = self.config.removal.sigma_clip

            # Remove the star in the frame
            star.remove(frame, self.config.removal, default_fwhm, sigma_clip)

    # *****************************************************************

    def remove_saturation(self, frame, galaxyextractor=None):

        """
        This function ...
        :param image:
        :return:
        """

        # Inform the user
        log.info("Removing saturation from the frame")

        # Get a list of the fwhm of all the stars that were fitted
        fwhm_list = self.fwhm_list

        # Determine the fwhm for the stars that were not fitted
        if self.config.saturation.no_source_fwhm == "max": default_fwhm = max(fwhm_list)
        elif self.config.saturation.no_source_fwhm == "mean": default_fwhm = np.mean(fwhm_list)
        elif self.config.saturation.no_source_fwhm == "median": default_fwhm = np.median(fwhm_list)
        else: raise ValueError("Unkown measure for determining the fwhm of stars without detected source")

        # Set the number of stars where saturation was removed to zero initially
        removed = 0

        # Detect and remove saturation for all stars
        if self.config.saturation.method == "all":

            # Inform the user on the number of stars that have a source
            log.debug("Number of stars with source = " + str(self.have_source))

            # Loop over all stars
            for star in self.objects:

                # If this star should be ignored, skip it
                if star.ignore: continue

                # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
                if galaxyextractor.principal.contains(star.pixel_position(frame.wcs)) and self.config.saturation.no_sigma_clip_on_galaxy: sigma_clip = False
                else: sigma_clip = self.config.saturation.sigma_clip

                # Find a saturation source and remove it from the frame
                success = star.remove_saturation(frame, self.config.saturation, default_fwhm, sigma_clip)
                if success: star.has_saturation = True
                removed += success

            # Inform the user
            log.debug("Removed saturation in " + str(removed) + " out of " + str(self.have_source) + " stars with source ({0:.2f}%)".format(removed/self.have_source*100.0))

        # Detect and remove saturation for the brightest stars
        elif self.config.saturation.method == "brightest":

            # TODO: allow "central_brightness" as config.saturation.criterion

            # Get a list of the fluxes of the stars
            flux_list = self.flux_list

            # Calculate the minimal flux/central brightness
            minimum = statistics.cutoff(flux_list, self.config.saturation.brightest_method, self.config.saturation.limit)

            # Inform the user
            quantity = "flux" if self.config.saturation.criterion == "flux" else "central brightness"
            log.debug("Minimum value of the " + quantity + " for saturation removal: {0:.2f}".format(minimum))

            # Inform the user
            eligible = len([flux for flux in flux_list if flux >= minimum])
            log.debug("Number of stars eligible for saturation removal: " + str(eligible) + " ({0:.2f}%)".format(eligible/len(self.objects)*100.0))

            # Inform the user on the number of stars that have a source
            log.debug("Number of stars with source = " + str(self.have_source))

            # Loop over all stars
            for star in self.objects:

                # If this star should be ignored, skip it
                if star.ignore: continue

                # If a source was not found for this star, skip it
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Calculate the value (flux or brightness) for this star
                try: value = star.flux
                except AttributeError: value = 0.0

                # Remove the saturation if the value is greater than the minimum value or the star has no source and 'remove_if_undetected' is enabled
                if value >= minimum or (self.config.saturation.remove_if_undetected and not star.has_source):

                    # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
                    if galaxyextractor.principal.contains(star.pixel_position(frame.wcs)) and self.config.saturation.no_sigma_clip_on_galaxy: sigma_clip = False
                    else: sigma_clip = self.config.saturation.sigma_clip

                    # Find a saturation source and remove it from the frame
                    success = star.remove_saturation(frame, self.config.saturation, default_fwhm, sigma_clip)
                    removed += success

            # Inform the user
            log.debug("Removed saturation in {0} out of {1} stars ({2:.2f}%)".format(removed, eligible, removed/eligible*100.0))

        # Unkown saturation
        else: raise ValueError("Unknown method (should be 'brightest' or 'all'")

    # *****************************************************************

    def find_apertures(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Loop over all stars
        for star in self.objects:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If the galaxy does not have a source, continue
            if star.has_saturation: star.find_aperture(sigma_level=self.config.apertures.sigma_level)

    # *****************************************************************

    def create_frame(self, frame):

        """
        This function ...
        """

        return None

    # *****************************************************************

    def create_region(self, frame, type="sky"):

        """
        This function ...
        :return:
        """

        # Initialize lists
        position_list = []
        radius_list = []
        color_list = []

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.objects:

            position_list.append(star.position)

            if star.has_model:

                fwhm = star.fwhm
                color = "green"

            else:

                fwhm = default_fwhm
                color = "red"

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * self.config.region.sigma_level

            # Add the radius (in arcseconds) and the color the appropriate list
            radius_list.append(radius * frame.pixelscale)
            color_list.append(color)

        # Create a region
        region = Region.circles(position_list, radius_list, color_list)

        if type == "sky": return region
        elif type == "image": return region.as_imagecoord(frame.wcs.to_header())
        else: raise ValueError("Type should be either 'sky' or 'image'")

    # *****************************************************************

    def write_region(self, frame, path, annotation=None):

        """
        This function ...
        :param frame:
        :return:
        """

        # Create a file
        f = open(path,'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.0", file=f)

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.objects:

            # Get the center in pixel coordinates
            x_center, y_center = star.position.to_pixel(frame.wcs, origin=0)

            if star.has_source:

                if star.has_model:

                    fwhm = star.fwhm
                    color = "blue"

                else:

                    fwhm = default_fwhm
                    color = "green"

            else:

                fwhm = default_fwhm
                color = "red"

            if annotation == "flux":

                if star.has_source and star.source.has_background:

                    text = "text = {" + str(int(round(star.flux))) + "}"

                else: text = ""

            elif annotation == "has_source":

                text = "text = {" + str(star.has_source) + "}"

            elif annotation == "has_background":

                if star.has_source: text = "text = {" + str(star.source.has_background) + "}"
                else: text = ""

            elif annotation is None: text = ""
            else: raise ValueError("Invalid option for annotation")

            color_suffix = " # color = " + color

            point_suffix = " # point = x " + text

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * self.config.region.sigma_level

            if star.has_source:

                if star.source.has_peak:

                    print("image;point({},{})".format(star.source.peak.x, star.source.peak.y) + point_suffix, file=f)

            # Show a circle for the star
            print("image;circle({},{},{})".format(x_center, y_center, radius) + color_suffix, file=f)

            # Aperture created from saturation mask
            if star.has_aperture:

                ap_x_center, ap_y_center = star.aperture.positions[0]
                major = star.aperture.a
                minor = star.aperture.b
                angle = star.aperture.theta * math.pi / 180

                aperture_suffix = " # color = white"

                print("image;ellipse({},{},{},{},{})".format(ap_x_center, ap_y_center, major, minor, angle) + aperture_suffix, file=f)

        # Close the file
        f.close()

    # *****************************************************************

    @property
    def have_model(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.objects: count += star.has_model
        return count

    # *****************************************************************

    @property
    def fwhm_list(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fwhm of the fitted stars
        fwhm_list = []

        # Loop over all stars
        for star in self.objects:

            # If the star contains a model, add the fwhm of that model to the list
            if star.has_model: fwhm_list.append(star.fwhm)

        # Return the list
        return fwhm_list

    # *****************************************************************

    @property
    def flux_list(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fluxes of the stars
        flux_list = []

        # Loop over all stars
        for star in self.objects:

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            if star.has_source and star.source.has_background:

                # Add the flux to the list
                flux_list.append(star.flux)

        # Return the list
        return flux_list

    # *****************************************************************

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fwhm of the fitted stars
        fwhm_list = []

        # Loop over all stars
        for star in self.objects:

            # If the star contains a model, add the fwhm of that model to the list
            if star.has_model: fwhm_list.append(star.fwhm)

        # Return the mean of the fwhm for all fitted stars
        return np.mean(fwhm_list)

    # *****************************************************************

    def table(self):

        """
        This function ...
        :return:
        """

        ids = []
        ascensions = []
        declinations = []
        sources = []
        models = []
        fwhms = []

        # Loop over all stars
        for star in self.objects:

            ids.append(star.ucac_id)
            ascensions.append(star.position.ra.value)
            declinations.append(star.position.dec.value)
            sources.append(star.has_source)
            models.append(star.has_model)
            if star.has_model: fwhms.append(star.fwhm)
            else: fwhms.append(None)

        # Create and return the table
        return Table([ids, ascensions, declinations, sources, models, fwhms], names=('UCAC-ID', 'RA', 'DEC', 'Source', 'Model', 'FWHM'), meta={'name': 'stars'})

    # *****************************************************************

    def plot(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Create a HDU from this frame with the image header
        hdu = pyfits.PrimaryHDU(frame, frame.wcs.to_header())

        # Create a figure canvas
        figure = plt.figure(figsize=(20, 20))

        # Create a figure from this frame
        plot = aplpy.FITSFigure(hdu, figure=figure)

        # Plot in color scale
        plot.show_colorscale()

        # Add a color bar if requested
        if self.config.plotting.show_colorbar: plot.add_colorbar()

        # Add these shapes to the plot
        plot.show_regions(self.create_region(frame, type="image"))

        # Show the plot
        plt.show()

# *****************************************************************
