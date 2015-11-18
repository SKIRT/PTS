#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
import matplotlib.pyplot as plt

# Import Astromagic modules
from .objectextraction import ObjectExtractor
from ..core.star import Star
from ..tools import statistics
from ..core.vector import Position, Extent
from ..core.regions import Region
from ..core.source import Source
from ..core.masks import Mask
from ..tools import configuration
from ..tools import fitting

# Import astronomical modules
import aplpy
import astropy.io.fits as pyfits
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

# *****************************************************************

class StarExtractor(ObjectExtractor):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        ## Configuration

        self.config = configuration.set("skyextractor", config)

        ## Base class

        # Call the constructor of the base class
        super(StarExtractor, self).__init__()

        ## Attributes

        # Initialize an empty list for NOMAD stars
        self.nomad_stars = []

    # *****************************************************************

    def run(self, frame, galaxyextractor=None):

        """
        This function ...
        """

        # Call the setup function
        self.setup(frame)

        # Create a list of stars based on online catalogs
        self.fetch_stars(galaxyextractor)

        # Set special stars
        if self.config.special_region is not None: self.set_special()

        # Set ignore stars
        if self.config.ignore_region is not None: self.set_ignore()

        # Set manual stars
        if self.config.manual_region is not None: self.set_manual()

        # For each star, find a corresponding source in the image
        self.find_sources()

        # Fit analytical models to the stars
        self.fit_stars()

        # If requested, remove the stars
        if self.config.remove: self.remove_stars(galaxyextractor)

        # If requested, remove saturation in the image
        if self.config.remove_saturation: self.remove_saturation(galaxyextractor)

        # If requested, find apertures
        if self.config.find_apertures: self.find_apertures()

        # If requested, remove apertures
        if self.config.remove_apertures: self.remove_apertures()

        # If a manual region was specified, remove the corresponding stars
        if self.config.manual_region is not None: self.remove_manual()

        # If requested, write a table with the galaxy properties
        if self.config.save_table: self.save_table()

        # If requested, save the star region
        if self.config.save_region: self.save_region()

        # If requested, save the frame where the stars are masked
        if self.config.save_masked_frame: self.save_masked_frame()

        # If requested, save the result
        if self.config.save_result: self.save_result()

    # *****************************************************************

    def fetch_stars(self, galaxyextractor=None):

        """
        This function ...
        """

        # Inform the user
        log.info("Fetching star positions from an online catalog")

        # Get the range of right ascension and declination of this image
        center, ra_span, dec_span = self.frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=self.config.fetching.catalog)
        table = result[0]

        # Loop over all stars in the table
        distances = []
        for entry in table:

            # Initialize an empty dictionary to contain the magnitudes and the magnitude errors
            magnitudes = {}
            mag_errors = {}

            # Get the ID of this star in the catalog
            if self.config.fetching.catalog == "UCAC4": star_id = entry["UCAC4"]
            elif self.config.fetching.catalog == "NOMAD": star_id = entry["NOMAD1"]
            elif self.config.fetching.catalog == "II/246": star_id = entry["_2MASS"]
            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # Get the position of the star as a SkyCoord object
            position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # If this star does not lie within the frame, skip it
            if not self.frame.contains(position): continue

            # Get the mean error on the right ascension and declination
            if self.config.fetching.catalog == "UCAC4" or self.config.fetching.catalog == "NOMAD":

                ra_error = entry["e_RAJ2000"]*u.mas
                dec_error = entry["e_DEJ2000"]*u.mas

            elif self.config.fetching.catalog == "II/246":

                error_maj = entry["errMaj"] * u.arcsec
                error_min = entry["errMin"] * u.arcsec
                error_theta = Angle(entry["errPA"], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = error_maj
                dec_error = error_maj

            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # Loop over all galaxies
            galaxies = galaxyextractor.objects if galaxyextractor is not None else []
            for galaxy in galaxies:

                # Calculate the pixel position of the galaxy
                galaxy_position = galaxy.pixel_position(self.frame.wcs)

                # Calculate the distance between the star's position and the galaxy's center
                x_center, y_center = position.to_pixel(self.frame.wcs)
                difference = galaxy_position - Position(x=x_center, y=y_center)

                # Add the star-galaxy distance to the list of distances
                distances.append(difference.norm)

                # The principal galaxy/galaxies
                if galaxy.principal:

                    # Check whether the star-galaxy distance is smaller than a certain threshold
                    if difference.norm <= self.config.fetching.min_distance_from_galaxy.principal:

                        # Remove the position of this galaxy from the list (one star is already identified with it)
                        #galaxies.remove(galaxy)
                        break

                # Companion galaxies
                elif galaxy.companion:

                    if difference.norm <= self.config.fetching.min_distance_from_galaxy.companion: break

                # All other galaxies in the frame
                else:

                    if difference.norm <= self.config.fetching.min_distance_from_galaxy.other: break

            # If a break is not encountered
            else:

                # Get the magnitude in different bands
                for name in entry.colnames:

                    # If this column name does not end with "mag", skip it
                    if not name.endswith("mag"): continue

                    # If the column name contains more than one character before "mag", skip it
                    if len(name.split("mag")[0]) > 1: continue

                    # Get the name of the band
                    band = name.split("mag")[0]

                    # Add the magnitude in this band to the dictionary
                    magnitudes[band] = entry[name] * u.mag

                    # Check whether an error on the magnitude is present
                    if "e_"+name in entry.colnames:

                        # If so, add it to the mag_errors dictionary
                        mag_errors[band] = entry["e_"+name] * u.mag

                # Create a star object
                star = Star(catalog=self.config.fetching.catalog, id=star_id, position=position, ra_error=ra_error,
                            dec_error=dec_error, magnitudes=magnitudes, magnitude_errors=mag_errors)

                # Check whether this star is on top of the galaxy, and label it so (by default, star.on_galaxy is False)
                if galaxyextractor is not None: star.on_galaxy = galaxyextractor.principal.contains(star.pixel_position(self.frame.wcs))

                # If requested, enable track record
                if self.config.track_record: star.enable_track_record()

                # Add the star to the list of stars
                self.objects.append(star)

        # Inform the user
        if galaxyextractor is not None: log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(distance) for distance in sorted(distances)[:10]))

        # Inform the user
        log.debug("Number of stars: " + str(len(self.objects)))

    # *****************************************************************

    def compare_nomad_stars(self, frame):

        """
        This function ...
        :return:
        """

        # Get the range of right ascension and declination of this image
        center, ra_span, dec_span = frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=["NOMAD"])
        table = result[0]

        for entry in table:

            # Get the right ascension and the declination
            position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            ra_error = entry["e_RAJ2000"]*u.mas
            dec_error = entry["e_DEJ2000"]*u.mas

            smallest_distance = float("inf")
            ra_error_pixels_ucac = None
            dec_error_pixels_ucac = None

            # Loop over the stars found in the other catalog
            for star in self.objects:

                # Create a pixel position for the NOMAD star
                x, y = position.to_pixel(frame.wcs, origin=0)
                pixel_position = Position(x, y)

                # Calculate the distance between the NOMAD star and the current other star
                difference = pixel_position - star.pixel_position(frame.wcs)
                distance = difference.norm

                # Check whether this is the smallest distance encountered so far
                if distance < smallest_distance:

                    smallest_distance = distance
                    ra_error_pixels_ucac = star.ra_error.to("arcsec") / frame.pixelscale
                    dec_error_pixels_ucac = star.dec_error.to("arcsec") / frame.pixelscale

            ra_error_pixels = ra_error.to("arcsec") / frame.pixelscale
            dec_error_pixels = dec_error.to("arcsec") / frame.pixelscale

            # Smallest distance
            print(smallest_distance, ra_error_pixels, dec_error_pixels, ra_error_pixels_ucac, dec_error_pixels_ucac)

    # *****************************************************************

    def fit_stars(self):

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
                center, radius, angle = star.ellipse_parameters(self.frame.wcs, self.frame.pixelscale, self.config.fitting.initial_radius)

                # Create a source object
                source = Source(self.frame, center, radius, angle, self.config.fitting.background_outer_factor)

            else: source = None

            # Find a source
            if star.has_source or source is not None: star.fit_model(self.config.fitting, source)

        # Inform the user
        log.debug("Found a model for {0} out of {1} stars with source ({2:.2f}%)".format(self.have_model, self.have_source, self.have_model/self.have_source*100.0))

    # *****************************************************************

    def remove_stars(self, galaxyextractor=None):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing the stars from the frame")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm

        # Inform the user
        log.debug("Default FWHM used when star could not be fitted: {0:.2f} pixels".format(default_fwhm))

        if self.config.removal.method == "model":

            # Calculate the relative differences between the ampltides of the fitted models and the corresponding sources
            differences = np.array(self.amplitude_differences) * 100.0

            print(np.mean(differences))
            print(np.median(differences))
            print(np.std(differences))
            print(differences)

        # Loop over all stars in the list
        for star in self.objects:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
            if not self.config.removal.remove_foreground and galaxyextractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

            # Remove the star in the frame
            star.remove(self.frame, self.config.removal, default_fwhm)

    # *****************************************************************

    def remove_saturation(self, galaxyextractor=None):

        """
        This function ...
        :param image:
        :return:
        """

        # Inform the user
        log.info("Removing saturation from the frame")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm

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

                # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
                if not self.config.saturation.remove_foreground and galaxyextractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

                # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Find a saturation source and remove it from the frame
                success = star.remove_saturation(self.frame, self.config.saturation, default_fwhm)
                if success: star.has_saturation = True
                removed += success

            # Inform the user
            log.debug("Removed saturation in " + str(removed) + " out of " + str(self.have_source) + " stars with source ({0:.2f}%)".format(removed/self.have_source*100.0))

        # Detect and remove saturation for the brightest stars
        elif self.config.saturation.method == "brightest":

            # TODO: allow "central_brightness" as config.saturation.criterion

            # Get a list of the fluxes of the stars
            flux_list = self.fluxes

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

                # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
                if not self.config.saturation.remove_foreground and galaxyextractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

                # If a source was not found for this star, skip it
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Calculate the value (flux or brightness) for this star
                try: value = star.flux
                except AttributeError: value = 0.0

                # Remove the saturation if the value is greater than the minimum value or the star has no source and 'remove_if_undetected' is enabled
                if value >= minimum or (self.config.saturation.remove_if_undetected and not star.has_source):

                    # Find a saturation source and remove it from the frame
                    success = star.remove_saturation(self.frame, self.config.saturation, default_fwhm)
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

        # Inform the user
        log.info("Constructing elliptical apertures regions to encompass saturated stars")

        # Loop over all stars
        for star in self.objects:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If the galaxy does not have a source, continue
            if star.has_saturation: star.find_aperture(sigma_level=self.config.apertures.sigma_level)

    # *****************************************************************

    def remove_apertures(self):

        """
        This function ...
        :param frame:
        :param factor:
        :return:
        """

        # Inform the user
        log.info("Replacing aperture regions by the estimated background")

        # Loop over all stars
        for star in self.objects:

            # If the object does not have an aperture, skip it
            if not star.has_aperture: continue

            # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
            if star.on_galaxy and self.config.saturation.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = self.config.saturation.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            if star.on_galaxy and self.config.saturation.polynomial_on_galaxy: interpolation_method = "polynomial"
            else: interpolation_method = self.config.saturation.interpolation_method

            # Expansion factor
            expansion_factor = self.config.aperture_removal.expansion_factor

            # Create a source object
            # Get the parameters of the elliptical aperture
            x_center, y_center = star.aperture.positions[0]
            center = Position(x=x_center, y=y_center)

            major = star.aperture.a * expansion_factor
            minor = star.aperture.b * expansion_factor

            radius = Extent(x=major, y=minor)

            # theta is in radians
            angle = Angle(star.aperture.theta, u.rad)

            # Create a source
            source = Source(self.frame, center, radius, angle, self.config.saturation.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background(interpolation_method, sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

    # *****************************************************************

    @property
    def region(self):

        """
        This function ...
        :return:
        """

        # TODO: improve this function

        type="sky"

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
            radius_list.append(radius * self.frame.pixelscale)
            color_list.append(color)

        # Create a region
        region = Region.circles(position_list, radius_list, color_list)

        if type == "sky": return region
        elif type == "image": return region.as_imagecoord(self.frame.wcs.to_header())
        else: raise ValueError("Type should be either 'sky' or 'image'")

    # *****************************************************************

    def write_region(self, path, annotation=None):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Writing stars region to " + path)

        # Create a file
        f = open(path,'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.objects:

            # Get the center in pixel coordinates
            x_center, y_center = star.position.to_pixel(self.frame.wcs, origin=0)

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

            # If the FWHM is defined, draw a circle for the star and draw a cross for its peak position (if defined)
            if fwhm is not None:

                # Calculate the radius in pixels
                radius = fwhm * statistics.fwhm_to_sigma * self.config.region.sigma_level

                # Draw a cross for the peak position
                if star.has_source and star.source.has_peak:

                    point_suffix = " # point = x " + text
                    print("image;point({},{})".format(star.source.peak.x, star.source.peak.y) + point_suffix, file=f)

                # Show a circle for the star
                color_suffix = " # color = " + color
                print("image;circle({},{},{})".format(x_center, y_center, radius) + color_suffix, file=f)

            # If the FWHM is undefined, simply draw a point for the star's position (e.g. when this function is called
            # after the fetch_stars method)
            else: print("image;point({},{})".format(x_center, y_center), file=f)

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
    def fwhms(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fwhm of the fitted stars
        fwhms = []

        # Loop over all stars
        for star in self.objects:

            # If the star contains a model, add the fwhm of that model to the list
            if star.has_model: fwhms.append(star.fwhm)

        # Return the list
        return fwhms

    # *****************************************************************

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fluxes of the stars
        fluxes = []

        # Loop over all stars
        for star in self.objects:

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            if star.has_source and star.source.has_background:

                # Add the flux to the list
                fluxes.append(star.flux)

        # Return the list
        return fluxes

    # *****************************************************************

    @property
    def amplitude_differences(self):

        """
        This function ...
        :return:
        """

        # Initialize
        differences = []

        # Loop over all stars
        for star in self.objects:

            # If the star was not fitted, skip it
            if not star.has_model: continue

            # Determine the amplitude and the position of the center of the model
            amplitude_model = star.model.amplitude
            center = star.source.cutout.rel_position(fitting.center(star.model))

            # Convert into integers
            x = int(round(center.x))
            y = int(round(center.y))

            # Calculate the value of the source at the model's center position
            amplitude_source = star.source.subtracted[y, x]

            # Calculate the difference of the amplitudes
            difference = abs(amplitude_model - amplitude_source)
            rel_difference = difference / amplitude_source

            # Add the relative difference to the list
            differences.append(rel_difference)

        # Return the list of differences
        return differences

    # *****************************************************************

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # If the list of FWHM values is empty (the stars were not fitted yet), return None
        if len(self.fwhms) == 0: return None

        # Determine the default FWHM and return it
        if self.config.fwhm.measure == "max": return max(self.fwhms)
        elif self.config.fwhm.measure == "mean": return np.mean(self.fwhms)
        elif self.config.fwhm.measure == "median": return np.median(self.fwhms)
        else: raise ValueError("Unkown measure for determining the default FWHM")

    # *****************************************************************

    @property
    def mask(self):

        """
        This function ...
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Loop over all stars
        for star in self.objects:

            # If a source was found for this star
            if star.has_source:

                # Check whether the aperture should be used for the mask
                if self.config.mask.use_aperture and star.has_aperture:

                    # Create a mask from the aperture of the object (expand if specified under self.config.aperture_removal)
                    object_mask_frame = Mask.from_aperture(self.frame.xsize, self.frame.ysize, star.aperture, expansion_factor=self.config.aperture_removal.expansion_factor)

                    # Now, we don't limit setting the mask within the source's cutout, because we expanded the apertures to perhaps a size larger than this cutout,
                    # so just add the object_mask_frame to the total frame
                    mask += object_mask_frame

                # Else, use the source mask (saturation or model-based source)
                else:

                    # Add this galaxy to the total mask
                    mask[star.source.cutout.y_min:star.source.cutout.y_max, star.source.cutout.x_min:star.source.cutout.x_max] += star.source.mask

            # If a source was not detected for this star but the inclusion of undetected stars in the mask is enabled
            elif self.config.mask.include_undetected:

                # Create a new source based on the default FWHM value and the specified sigma level and outer factor of the removal
                source = star.source_at_sigma_level(self.frame, self.fwhm, self.config.removal.sigma_level, self.config.removal.outer_factor)

                # Add this sky object to the total mask
                mask[source.cutout.y_min:source.cutout.y_max, source.cutout.x_min:source.cutout.x_max] += source.mask

        # Return the mask
        return mask

    # *****************************************************************

    @property
    def saturation_mask(self):

        """
        This function ...
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Loop over all sky objects
        for star in self.objects:

            # If no saturation was found for this star, skip it
            if not star.has_saturation: continue

            # Add this star to the mask
            if self.config.saturation_mask.use_aperture and star.has_aperture:

                object_mask_frame = Mask.from_aperture(self.frame.xsize, self.frame.ysize, star.aperture)
                object_mask = object_mask_frame[star.source.cutout.y_min:star.source.cutout.y_max, star.source.cutout.x_min:star.source.cutout.x_max]

            else: object_mask = star.source.mask

            # Add this galaxy to the total mask
            mask[star.source.cutout.y_min:star.source.cutout.y_max, star.source.cutout.x_min:star.source.cutout.x_max] += object_mask

        # Return the mask
        return mask

    # *****************************************************************

    @property
    def table(self):

        """
        This function ...
        :return:
        """

        # Initialize empty lists for the table columns
        catalogs = []
        ids = []
        ascensions = []
        declinations = []
        sources = []
        models = []
        fwhms = []

        # Loop over all stars
        for star in self.objects:

            catalogs.append(star.catalog)
            ids.append(star.id)
            ascensions.append(star.position.ra.value)
            declinations.append(star.position.dec.value)
            sources.append(star.has_source)
            models.append(star.has_model)
            if star.has_model: fwhms.append(star.fwhm)
            else: fwhms.append(None)

        # Create and return the table
        return Table([catalogs, ids, ascensions, declinations, sources, models, fwhms], names=('CATALOG', 'ID', 'RA', 'DEC', 'Source', 'Model', 'FWHM'), meta={'name': 'stars'})

    # *****************************************************************

    def plot(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Create a HDU from this frame with the image header
        hdu = pyfits.PrimaryHDU(self.frame, self.frame.wcs.to_header())

        # Create a figure canvas
        figure = plt.figure(figsize=(20, 20))

        # Create a figure from this frame
        plot = aplpy.FITSFigure(hdu, figure=figure)

        # Plot in color scale
        plot.show_colorscale()

        # Add a color bar if requested
        if self.config.plotting.show_colorbar: plot.add_colorbar()

        # Add these shapes to the plot
        plot.show_regions(self.region)

        # Show the plot
        plt.show()

# *****************************************************************
