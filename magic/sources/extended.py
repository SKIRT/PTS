#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.extended Contains the ExtendedSourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import traceback
from collections import OrderedDict

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..region.list import PixelRegionList, SkyRegionList
from ..basics.stretch import PixelStretch
from ..region.point import PixelPointRegion
from ..region.ellipse import PixelEllipseRegion, SkyEllipseRegion
from ..core.frame import Frame
from ...core.basics.configurable import Configurable
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.basics.table import SmartTable
from ..core.mask import Mask
from ..basics.coordinate import SkyCoordinate

# -----------------------------------------------------------------

class ExtendedSourceTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["RA"] = (float, "deg", "right ascension")
    _column_info["DEC"] = (float, "deg", "declination")
    _column_info["Detected"] = (bool, None, "Has source detected")
    _column_info["Flux"] = (float, "Jy", "flux for the point source")
    _column_info["Flux error"] = (float, "Jy", "error on the flux value")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExtendedSourceTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_source(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        if source is not None:

            # Inform the user
            log.info("Adding source " + str(source.index) + " to the table of extended sources ...")

            # Get extended source properties
            ra = source.position.ra
            dec = source.position.dec
            detected = source.has_detection
            flux = None
            flux_error = None

            # Construct the row
            values = [ra, dec, detected, flux, flux_error]

        else: values = [None, None, None, None, None]

        # Add a row
        self.add_row(values)

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index] * self["RA"].unit, dec=self["DEC"][index] * self["DEC"].unit, unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def is_detected(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Detected"][index]

    # -----------------------------------------------------------------

    def get_flux(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Flux", index)

    # -----------------------------------------------------------------

    def get_flux_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Flux error", index)

# -----------------------------------------------------------------

class ExtendedSourceFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExtendedSourceFinder, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Initialize the sources list
        self.sources = []

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

        # The region list
        self.regions = None

        # The segmentation map
        self.segments = None

        # The extended source table
        self.table = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        """

        # 2. Load the sources from the catalog
        self.load_sources()

        # 3. Detect the sources
        if not self.config.weak: self.detect_sources()
        else: self.set_detections()

        # Find apertures
        #if self.config.find_apertures: self.find_contours()

        # 4. Create the region list
        self.create_regions()

        # 5. Create the segmentation map
        self.create_segments()

        # 10. Create the table
        self.create_table()

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the setup function of the base class
        super(ExtendedSourceFinder, self).setup()

        # Inform the user
        log.info("Setting up the extended source finder ...")

        # Make a local reference to the image frame and catalog
        self.frame = kwargs.pop("frame")
        self.catalog = kwargs.pop("catalog")

        # Masks
        self.special_mask = kwargs.pop("special_mask", None)
        self.ignore_mask = kwargs.pop("ignore_mask", None)
        self.bad_mask = kwargs.pop("bad_mask", None)

        # Create an empty frame for the segments
        self.segments = Frame.zeros_like(self.frame)

        # Initialize the table
        self.table = ExtendedSourceTable()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the extended source finder ...")

        # Create a new list
        self.sources = []

        # Clear the frame
        self.frame = None

    # -----------------------------------------------------------------

    @property
    def principal(self):

        """
        This function ...
        :return:
        """

        for source in self.sources:
            if source is None: continue
            if source.principal: return source
        return None

    # -----------------------------------------------------------------

    @property
    def companions(self):

        """
        This function ...
        :return:
        """

        companions = []
        for source in self.sources:
            if source is None: continue
            if source.companion: companions.append(source)
        return companions

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the extended sources from the catalog ...")

        # Loop over the entries in the catalog
        for index in range(len(self.catalog)):

            # Get position
            position = self.catalog.get_position(index)

            # If the source falls outside of the frame, skip it
            if not self.frame.contains(position): source = None
            else:

                # Calculate the pixel position of the galaxy in the frame
                pixel_position = position.to_pixel(self.frame.wcs)

                # Create a source
                source = self.catalog.create_source(index)

                # Enable track record if requested
                #if self.config.track_record: galaxy.enable_track_record()

                # Set attributes based on masks (special and ignore)
                if self.special_mask is not None: source.special = self.special_mask.masks(pixel_position)
                if self.ignore_mask is not None: source.ignore = self.ignore_mask.masks(pixel_position)

                # If the input mask masks this galaxy's position, set to None
                if self.bad_mask is not None and self.bad_mask.masks(pixel_position) and not source.principal: source = None

            # Add the new source to the list
            self.sources.append(source)

        # Debug messages
        if self.principal is not None:
            log.debug(self.principal.name + " is the principal galaxy in the frame")
            log.debug("The following galaxies are its companions: " + str(self.principal.companions))
        else: log.warning("No principal galaxy found in the frame")

    # -----------------------------------------------------------------

    def detect_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Detecting the sources in the frame ...")

        # Loop over all sources
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this source should be ignored, skip it
            if source.ignore: continue

            # If the galaxy is the principal galaxy and a region file is given
            if source.principal and self.config.principal_region is not None:

                # Load the principal galaxy region file
                region = SkyRegionList.from_file(self.config.principal_region)
                shape = region[0].to_pixel(self.frame.wcs)

                # Create a detection for the galaxy from the shape in the region file
                outer_factor = self.config.detection.background_outer_factor
                source.detection_from_shape(self.frame, shape, outer_factor)

            else:

                # If requested, use the galaxy extents obtained from the catalog to create the source
                if self.config.detection.use_d25 and source.has_extent:

                    outer_factor = self.config.detection.background_outer_factor
                    expansion_factor = self.config.detection.d25_expansion_factor
                    source.detection_from_parameters(self.frame, outer_factor, expansion_factor)

                else:

                    # Detect the source
                    try: source.detect(self.frame, self.config.detection)
                    except Exception as e:
                        #import traceback
                        log.error("Error during detection")
                        print(type(e))
                        print(e)
                        traceback.print_exc()
                        #if self.config.plot_track_record_if_exception:
                        #    if source.has_track_record: source.track_record.plot()
                        #    else: log.warning("Track record is not enabled")

            # If a source was not found for the principal or companion galaxies, force it
            outer_factor = self.config.detection.background_outer_factor
            if source.principal and not source.has_detection: source.detection_from_parameters(self.frame, outer_factor)
            elif source.companion and not source.has_detection and source.has_extent: source.detection_from_parameters(self.frame, outer_factor)

        # Inform the user
        log.info("Found a detection for {0} out of {1} objects ({2:.2f}%)".format(self.have_detection, len(self.sources), self.have_source/len(self.sources)*100.0))

    # -----------------------------------------------------------------

    def set_detections(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the detections ...")

        # Loop over all sources
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this source should be ignored, skip it
            if source.ignore: continue

            # If the galaxy is the principal galaxy and a region file is given
            if source.principal and self.config.principal_region is not None:

                # Load the principal galaxy region file
                region = SkyRegionList.from_file(self.config.principal_region)
                shape = region[0].to_pixel(self.frame.wcs)

                # Create a detection for the galaxy from the shape in the region file
                outer_factor = self.config.detection.background_outer_factor
                source.detection_from_shape(self.frame, shape, outer_factor)

            else:

                # If requested, use the galaxy extents obtained from the catalog to create the source
                if source.has_extent:

                    outer_factor = self.config.detection.background_outer_factor
                    expansion_factor = self.config.detection.d25_expansion_factor
                    source.detection_from_parameters(self.frame, outer_factor, expansion_factor)

                # The galaxy has no extent, use a standard radius of 20 pixels
                else:

                    default_radius = self.config.region.default_radius

                    outer_factor = self.config.detection.background_outer_factor
                    shape = source.ellipse(self.frame.wcs, default_radius)
                    source.detection_from_shape(self.frame, shape, outer_factor)

    # -----------------------------------------------------------------

    def find_contours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing elliptical contours to encompass the detected galaxies ...")

        # Loop over all galaxies
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this galaxy should be ignored, skip it
            if source.ignore: continue

            # If the galaxy does not have a source, continue
            if source.has_source: source.find_contour(self.frame, self.config.apertures)

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

        x_radius = 0.5 * self.principal.major.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value
        y_radius = 0.5 * self.principal.minor.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value
        radius = PixelStretch(x_radius, y_radius)

        # Create and return an ellipse
        return PixelEllipseRegion(center, radius, angle)

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
        sky_ellipse = SkyEllipseRegion.from_pixel(ellipse, self.frame.wcs)

        # Return the sky ellipse
        return sky_ellipse

    # -----------------------------------------------------------------

    @property
    def principal_mask(self):

        """
        This function ...
        :return:
        """

        if self.principal is None: return Mask.empty_like(self.frame)

        #return self.galaxies.get_principal_mask(self.frame)

        # Create a new mask with the dimensions of the frame
        mask = Mask.empty_like(self.frame)

        # Add the principal galaxy's mask to the total mask
        mask[self.principal.detection.cutout.y_slice, self.principal.detection.cutout.x_slice] = self.principal.detection.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @property
    def companion_mask(self):

        """
        This function ...
        :return:
        """

        #return self.galaxies.get_companion_mask(self.frame)

        # Create a new mask with the dimension of the frame
        mask = Mask.empty_like(self.frame)

        # Loop over all companion galaxies
        for source in self.companions:

            try:
                # Check if the galaxy has a source and add its mask to the total mask
                if source.has_detection: mask[source.detection.cutout.y_slice, source.detection.cutout.x_slice] = source.detection.mask
            except IndexError: pass

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def create_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating regions ...")

        # Initialize the region
        self.regions = PixelRegionList()

        # Loop over all sources
        for source in self.sources:

            if source is None: continue

            # Get the center in pixel coordinates
            center = source.pixel_position(self.frame.wcs)

            # Set the angle
            angle = source.pa_for_wcs(self.frame.wcs).to("deg") if source.pa is not None else Angle(0.0, "deg")

            if source.major is None:

                color = "red"

                x_radius = self.config.region.default_radius
                y_radius = self.config.region.default_radius

            elif source.minor is None or source.pa is None:

                color = "green"

                x_radius = 0.5 * source.major.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value
                y_radius = x_radius

            else:

                color = "green"

                x_radius = 0.5 * source.major.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value
                y_radius = 0.5 * source.minor.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value

            radius = PixelStretch(x_radius, y_radius)

            # Create a coordinate for the center and add it to the region
            meta = {"point": "x"}
            self.regions.append(PixelPointRegion(center.x, center.y, meta=meta))

            text = source.name
            if source.principal: text += " (principal)"

            # If hand-drawn principal region
            if source.principal and self.config.principal_region is not None: shape = source.shape

            # Create an ellipse for the galaxy
            else: shape = PixelEllipseRegion(center, radius, angle, meta=meta)

            # Set meta information
            meta = {"text": text, "color": color, "index": source.index}
            shape.meta = meta

            # Add the shape to the region list
            self.regions.append(shape)

    # -----------------------------------------------------------------

    def create_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the segmentation map ...")

        # Loop over all source
        for source in self.sources:

            # No source: continue
            if source is None: continue

            # Skip galaxies without source
            if not source.has_detection: continue

            # Determine the label for the galaxy
            if source.principal: label = 1
            elif source.companion: label = 2
            else: label = 3

            # Add the galaxy mask to the segmentation map
            self.segments[source.detection.y_slice, source.detection.x_slice][source.detection.mask] = label

    # -----------------------------------------------------------------

    def create_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the table of extended sources ...")

        # Loop over the sources
        for source in self.sources:

            # No source?
            #if source is None: continue

            # Add source
            self.table.add_source(source)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write table
        self.write_table()

        # Write regions
        self.write_regions()

        # Write the segmentation maps
        self.write_segments()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table ...")

        # Determine the path
        path = self.output_path_file("extended_sources.dat")

        # Write
        self.table.saveto(path)

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_cutouts(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the cutouts directory
        #directory_path = self.full_output_path(self.config.writing.cutouts_path)

        directory_path = fs.join(self.config.output_path, self.config.writing.cutouts_path)

        # Inform the user
        log.info("Writing cutout boxes to " + directory_path + " ...")

        # Keep track of the number of stars encountered
        principals = 0
        companions = 0
        with_source = 0

        # Loop over all sources
        for source in self.sources:

            if source is None: continue

            # Check if this is the principal galaxy source
            if source.principal:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_principal_" + str(principals) + ".fits")
                source.detection.saveto(path, origin=self.name)

                # Increment the counter of the number of principal galaxies (there should only be one, really...)
                principals += 1

            # Check if this is a companion galaxy
            elif source.companion:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_companion_" + str(companions) + ".fits")
                source.detection.saveto(path, origin=self.name)

                # Increment the counter of the number of companion galaxies
                companions += 1

            # Check if this galaxy has a source
            elif source.has_detection:

                # Save the cutout as a FITS file
                path = fs.join(directory_path, "galaxy_source_" + str(principals) + ".fits")
                source.detection.saveto(path, origin=self.name)

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
        for source in self.sources:
            if source is None: continue
            count += 1
        return count

    # -----------------------------------------------------------------

    @property
    def have_detection(self):

        """
        This function ...
        :return:
        """

        count = 0
        for source in self.sources:
            if source is None: continue
            if source.ignore: continue
            count += source.has_detection
        return count

# -----------------------------------------------------------------
