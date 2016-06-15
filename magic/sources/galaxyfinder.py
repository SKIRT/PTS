#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.galaxyfinder Contains the GalaxyFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..basics.mask import Mask
from ..basics.region import Region
from ..basics.skyregion import SkyRegion
from ..basics.vector import Extent
from ..basics.geometry import Coordinate, Ellipse
from ..basics.skygeometry import SkyCoordinate
from ..core.frame import Frame
from ..object.galaxy import Galaxy
from ..basics.skygeometry import SkyEllipse
from ...core.basics.configurable import Configurable
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

class GalaxyFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(GalaxyFinder, self).__init__(config, "magic")

        # -- Attributes --

        # Initialize an empty list for the galaxies
        self.galaxies = []

        # The image frame
        self.frame = None

        # The mask covering objects that require special attentation (visual feedback)
        self.special_mask = None

        # The mask covering pixels that should be ignored
        self.ignore_mask = None

        # The mask of bad pixels
        self.bad_mask = None

        # The galactic catalog
        self.catalog = None

        # The galactic statistics
        self.statistics = None

        # The galaxy region
        self.region = None

        # The segmentation map
        self.segments = None

    # -----------------------------------------------------------------

    def run(self, frame, catalog, special=None, ignore=None, bad=None):

        """
        This function ...
        :param frame:
        :param catalog:
        :param special:
        :param ignore:
        :param bad:
        """

        # 1. Call the setup function
        self.setup(frame, catalog, special, ignore, bad)

        # 2. Find the galaxies
        self.find_galaxies()

        # 3. Set the statistics
        self.set_statistics()

        # 4. Create the region
        self.create_region()

        # 5. Create the segmentation map
        self.create_segments()

    # -----------------------------------------------------------------

    def setup(self, frame, catalog, special_mask=None, ignore_mask=None, bad_mask=None):

        """
        This function ...
        :param frame:
        :param catalog
        :param special_mask:
        :param ignore_mask:
        :param bad_mask:
        """

        # Call the setup function of the base class
        super(GalaxyFinder, self).setup()

        # Inform the user
        log.info("Setting up the galaxy extractor ...")

        # Make a local reference to the image frame
        self.frame = frame

        self.catalog = catalog

        # Masks
        self.special_mask = special_mask
        self.ignore_mask = ignore_mask
        self.bad_mask = bad_mask

        # Create an empty frame for the segments
        self.segments = Frame.zeros_like(self.frame)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the galaxy extractor ...")

        # Clear the list of galaxies
        self.galaxies = []

        # Clear the frame
        self.frame = None

    # -----------------------------------------------------------------

    def find_galaxies(self):

        """
        This function ...
        :return:
        """

        # Load the galaxies from the galactic catalog
        self.load_galaxies()

        # Find the sources
        self.find_sources()

        # Find apertures
        if self.config.find_apertures: self.find_contours()

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

        # If the principal galaxy is not determined, return None
        return None

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
        log.info("Looking for sources near the galaxy positions ...")

        # Loop over all galaxies in the list
        for galaxy in self.galaxies:

            # If this sky object should be ignored, skip it
            if galaxy.ignore: continue

            # If the galaxy is the principal galaxy and a region file is
            if galaxy.principal and self.config.principal_region is not None:

                # Load the principal galaxy region file
                region = SkyRegion.from_file(self.config.principal_region)
                shape = region[0].to_pixel(self.frame.wcs)

                # Create a source for the galaxy from the shape in the region file
                outer_factor = self.config.detection.background_outer_factor
                galaxy.source_from_shape(self.frame, shape, outer_factor)

            else:

                # If requested, use the galaxy extents obtained from the catalog to create the source
                if self.config.detection.use_d25 and galaxy.has_extent:

                    outer_factor = self.config.detection.background_outer_factor
                    expansion_factor = self.config.detection.d25_expansion_factor
                    galaxy.source_from_parameters(self.frame, outer_factor, expansion_factor)

                else:

                    # Find a source
                    try: galaxy.find_source(self.frame, self.config.detection)
                    except Exception as e:
                        #import traceback
                        log.error("Error when finding source")
                        print(type(e))
                        print(e)
                        #traceback.print_exc()
                        if self.config.plot_track_record_if_exception:
                            if galaxy.has_track_record: galaxy.track_record.plot()
                            else: log.warning("Track record is not enabled")

            # If a source was not found for the principal or companion galaxies, force it
            outer_factor = self.config.detection.background_outer_factor
            if galaxy.principal and not galaxy.has_source: galaxy.source_from_parameters(self.frame, outer_factor)
            elif galaxy.companion and not galaxy.has_source and galaxy.has_extent: galaxy.source_from_parameters(self.frame, outer_factor)

        # Inform the user
        log.info("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.galaxies), self.have_source/len(self.galaxies)*100.0))

    # -----------------------------------------------------------------

    def load_galaxies(self):

        """
        This function creates the galaxy list from the galaxy catalog.
        :return:
        """

        # Inform the user
        log.info("Loading the galaxies from the catalog ...")

        # Create the list of galaxies
        for i in range(len(self.catalog)):

            # Get the galaxy properties
            name = self.catalog["Name"][i]
            redshift = self.catalog["Redshift"][i] if not self.catalog["Redshift"].mask[i] else None
            galaxy_type = self.catalog["Type"][i] if not self.catalog["Type"].mask[i] else None
            distance = self.catalog["Distance"][i] * Unit("Mpc") if not self.catalog["Distance"].mask[i] else None
            inclination = Angle(self.catalog["Inclination"][i], Unit("deg")) if not self.catalog["Inclination"].mask[i] else None
            d25 = self.catalog["D25"][i] * Unit("arcmin") if not self.catalog["D25"].mask[i] else None
            major = self.catalog["Major axis length"][i] * Unit("arcmin") if not self.catalog["Major axis length"].mask[i] else None
            minor = self.catalog["Minor axis length"][i] * Unit("arcmin") if not self.catalog["Minor axis length"].mask[i] else None
            position_angle = Angle(self.catalog["Position angle"][i], Unit("deg")) if not self.catalog["Position angle"].mask[i] else None
            ra = self.catalog["Right ascension"][i]
            dec = self.catalog["Declination"][i]
            names = self.catalog["Alternative names"][i].split(", ") if not self.catalog["Alternative names"].mask[i] else []
            principal = self.catalog["Principal"][i]
            companions = self.catalog["Companion galaxies"][i].split(", ") if not self.catalog["Companion galaxies"].mask[i] else []
            parent = self.catalog["Parent galaxy"][i] if not self.catalog["Parent galaxy"].mask[i] else None

            # Create a SkyCoordinate for the galaxy center position
            position = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame="fk5")

            # If the galaxy falls outside of the frame, skip it
            if not self.frame.contains(position): continue

            # Create a new Galaxy instance
            galaxy = Galaxy(i, name, position, redshift, galaxy_type, names, distance, inclination, d25, major, minor, position_angle)

            # Calculate the pixel position of the galaxy in the frame
            pixel_position = galaxy.pixel_position(self.frame.wcs)

            # Set other attributes
            galaxy.principal = principal
            galaxy.companion = parent is not None
            galaxy.companions = companions
            galaxy.parent = parent

            # Enable track record if requested
            if self.config.track_record: galaxy.enable_track_record()

            # Set attributes based on masks (special and ignore)
            if self.special_mask is not None: galaxy.special = self.special_mask.masks(pixel_position)
            if self.ignore_mask is not None: galaxy.ignore = self.ignore_mask.masks(pixel_position)

            # If the input mask masks this galaxy's position, skip it (don't add it to the list of galaxies)
            if self.bad_mask is not None and self.bad_mask.masks(pixel_position) and not galaxy.principal: continue

            # Add the new galaxy to the list
            self.galaxies.append(galaxy)

        # Debug messages
        log.debug(self.principal.name + " is the principal galaxy in the frame")
        log.debug("The following galaxies are its companions: " + str(self.principal.companions))

    # -----------------------------------------------------------------

    def find_contours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing elliptical contours to encompass the detected galaxies ...")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If this galaxy should be ignored, skip it
            if galaxy.ignore: continue

            # If the galaxy does not have a source, continue
            if galaxy.has_source: galaxy.find_contour(self.frame, self.config.apertures)

    # -----------------------------------------------------------------

    @property
    def principal_shape(self):

        """
        This function ...
        :return:
        """

        return self.principal.shape

    # -----------------------------------------------------------------

    @property
    def principal_ellipse(self):

        """
        This function ...
        :return:
        """

        # Get the center in pixel coordinates
        center = self.principal.pixel_position(self.frame.wcs)

        # Get the angle
        angle = self.principal.pa_for_wcs(self.frame.wcs)

        x_radius = 0.5 * self.principal.major.to("arcsec").value / self.frame.xy_average_pixelscale.to("arcsec/pix").value
        y_radius = 0.5 * self.principal.minor.to("arcsec").value / self.frame.xy_average_pixelscale.to("arcsec/pix").value
        radius = Extent(x_radius, y_radius)

        # Create and return an ellipse
        return Ellipse(center, radius, angle)

    # -----------------------------------------------------------------

    @property
    def principal_sky_ellipse(self):

        """
        This function ...
        :return:
        """

        # Get the ellipse in image coordinates
        ellipse = self.principal_ellipse

        # Create a SkyEllipse
        sky_ellipse = SkyEllipse.from_pixel(ellipse, self.frame.wcs)

        # Return the sky ellipse
        return sky_ellipse

    # -----------------------------------------------------------------

    @property
    def principal_mask(self):

        """
        This function ...
        :return:
        """

        # Create a new mask with the dimensions of the frame
        mask = Mask.empty_like(self.frame)

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
        mask = Mask.empty_like(self.frame)

        # Loop over all companion galaxies
        for galaxy in self.companions:

            # Check if the galaxy has a source and add its mask to the total mask
            if galaxy.has_source: mask[galaxy.source.cutout.y_slice, galaxy.source.cutout.x_slice] = galaxy.source.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def set_statistics(self):

        """
        This function ...
        :return:
        """

        index_column = []
        have_source_column = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            index_column.append(galaxy.index)
            have_source_column.append(galaxy.has_source)

        # Create data structure and set column names
        data = [index_column, have_source_column]
        names = ["Galaxy index", "Detected"]

        # Create the statistics table
        self.statistics = tables.new(data, names)

    # -----------------------------------------------------------------

    def create_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating galaxy region ...")

        # Initialize the region
        self.region = Region()

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Get the center in pixel coordinates
            center = galaxy.pixel_position(self.frame.wcs)

            # Set the angle
            angle = galaxy.pa_for_wcs(self.frame.wcs).to("deg") if galaxy.pa is not None else 0.0

            if galaxy.major is None:

                color = "red"

                x_radius = self.config.region.default_radius
                y_radius = self.config.region.default_radius

            elif galaxy.minor is None or galaxy.pa is None:

                color = "green"

                x_radius = 0.5 * galaxy.major.to("arcsec").value / self.frame.xy_average_pixelscale.to("arcsec/pix").value
                y_radius = x_radius

            else:

                color = "green"

                x_radius = 0.5 * galaxy.major.to("arcsec").value / self.frame.xy_average_pixelscale.to("arcsec/pix").value
                y_radius = 0.5 * galaxy.minor.to("arcsec").value / self.frame.xy_average_pixelscale.to("arcsec/pix").value

            radius = Extent(x_radius, y_radius)

            # Add point for the center
            #print("image;point({},{})".format(center.x+1, center.y+1) + point_suffix, file=f)
            #print("image;ellipse({},{},{},{},{})".format(center.x+1, center.y+1, x_radius, y_radius, angle) + color_suffix, file=f)

            # Create a coordinate for the center and add it to the region
            meta = {"point": "x"}
            self.region.append(Coordinate(center.x, center.y, meta=meta))

            # Create an ellipse for the galaxy and add it to the region
            meta = {"text": galaxy.name, "color": color}
            shape = Ellipse(center, radius, angle, meta=meta)
            self.region.append(shape)

            # Add aperture
            #if galaxy.has_contour:

                #contour_center = galaxy.contour.center
                #major = galaxy.contour.major
                #minor = galaxy.contour.minor
                #angle = galaxy.contour.angle.degree

                #aperture_suffix = " # color = white"

                #print("image;ellipse({},{},{},{},{})".format(contour_center.x+1, contour_center.y+1, major, minor, angle) + aperture_suffix, file=f)

    # -----------------------------------------------------------------

    def create_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the segmentation map of galaxies ...")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Skip galaxies without source
            if not galaxy.has_source: continue

            # Determine the label for the galaxy
            if galaxy.principal: label = 1
            elif galaxy.companion: label = 2
            else: label = 3

            # Add the galaxy mask to the segmentation map
            self.segments[galaxy.source.y_slice, galaxy.source.x_slice][galaxy.source.mask] = label

    # -----------------------------------------------------------------

    def write_cutouts(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the cutouts directory
        directory_path = self.full_output_path(self.config.writing.cutouts_path)

        # Inform the user
        log.info("Writing cutout boxes to " + directory_path + " ...")

        # Keep track of the number of stars encountered
        principals = 0
        companions = 0
        with_source = 0

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Check if this is the principal galaxy
            if galaxy.principal:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_principal_" + str(principals) + ".fits")
                galaxy.source.save(path, origin=self.name)

                # Increment the counter of the number of principal galaxies (there should only be one, really...)
                principals += 1

            # Check if this is a companion galaxy
            elif galaxy.companion:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_companion_" + str(companions) + ".fits")
                galaxy.source.save(path, origin=self.name)

                # Increment the counter of the number of companion galaxies
                companions += 1

            # Check if this galaxy has a source
            elif galaxy.has_source:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_source_" + str(principals) + ".fits")
                galaxy.source.save(path, origin=self.name)

                # Increment the counter of the number of galaxies with a source
                with_source += 1

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
