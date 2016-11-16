#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.star Contains the abstract Star class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical units
from astropy import units as u
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .skyobject import SkyObject
from ..core.box import CutoutMask, Box
from ..core.source import Source
from ..tools import statistics, fitting, masks, plotting
from ..analysis import sources
from ..region.ellipse import PixelEllipseRegion
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Star(SkyObject):

    """
    This class ...
    """

    def __init__(self, index, catalog=None, id=None, position=None, ra_error=None, dec_error=None, magnitudes=None, magnitude_errors=None, on_galaxy=False):

        """
        The constructor ...
        :return:
        """

        # Set the attributes
        self.index = index
        self.catalog = catalog
        self.id = id
        self.ra_error = ra_error
        self.dec_error = dec_error
        self.magnitudes = magnitudes
        self.magnitude_errors = magnitude_errors
        self.on_galaxy = on_galaxy

        self.confidence_level = 1

        # Set the model attribute to None initially
        self.model = None

        # The saturation source
        self.saturation = None

        # Call the constructor of the base class
        super(Star, self).__init__(position)

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.model is not None

    # -----------------------------------------------------------------

    @property
    def has_saturation(self):

        """
        This function ...
        :return:
        """

        return self.saturation is not None

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # Return the fwhm value of the model
        return fitting.fwhm(self.model)

    # -----------------------------------------------------------------

    @property
    def flux(self):

        """
        This function ...
        :return:
        """

        # Return the flux of the source
        return self.source.flux

    # -----------------------------------------------------------------

    def get_flux(self, without_background=False):

        """
        This function ...
        :param without_background:
        :return:
        """

        return self.source.get_flux(without_background)

    # -----------------------------------------------------------------

    def ellipse(self, wcs, default_radius):

        """
        This function ...
        :param wcs:
        :param default_radius:
        :return:
        """

        center, radius, angle = self.ellipse_parameters(wcs, default_radius)
        return PixelEllipseRegion(center, radius, angle)

    # -----------------------------------------------------------------

    def ellipse_parameters(self, wcs, default_radius):

        """
        This function ...
        :param wcs:
        :param default_radius:
        :return:
        """

        # Return the parameters
        return self.pixel_position(wcs), default_radius, Angle(0.0, u.Unit("deg"))

    # -----------------------------------------------------------------

    def fit_model(self, config, source=None, debug=False):

        """
        This function ...
        :param config:
        :param source:
        :param debug:
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("fitting")

        # Fit model to the source, in a loop over different analytical forms for the model
        for level in range(len(config.model_names)):

            if source is None:

                # Do the fitting
                source, model = sources.fit_model_to_source(self.source, config, self.track_record, level=level, special=debug)

            else:

                source, model = sources.fit_model_to_source(source, config, self.track_record, level=level)

            # If a model was found, set the attributes of the star object and exit the loop
            if model is not None:

                self.source = source
                self.model = model
                break

    # -----------------------------------------------------------------

    def source_from_shape(self, frame, shape, outer_factor):

        """
        This function ...
        :param frame:
        :param shape:
        :param outer_factor:
        :return:
        """

        # Create the source
        self.source = Source.from_shape(frame, shape, outer_factor)

    # -----------------------------------------------------------------

    def source_at_sigma_level(self, frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=False, shape=None):

        """
        This function ...
        :param frame:
        :param default_fwhm:
        :param sigma_level:
        :param outer_factor:
        :param use_default_fwhm:
        :param shape:
        :return:
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius of the contour in which the star will be removed
        if self.model is None or use_default_fwhm: radius = default_sigma * sigma_level
        else: radius = fitting.sigma(self.model) * sigma_level

        # Determine the center position of the source (center of model if present, otherwise position of the star)
        if self.source is not None:

            # If the star has been modeled succesfully, use the center position of the model
            # Otherwise, use the source's peak
            if self.model is not None: center = fitting.center(self.model)
            elif self.source.has_peak: center = self.source.peak
            else:

                log.warning("Star source does not have peak")
                center = self.pixel_position(frame.wcs)

        else:

            # Calculate the pixel coordinate of the star's position
            center = self.pixel_position(frame.wcs)

        # Create the new source
        ellipse = PixelEllipseRegion(center, radius)
        source = Source.from_ellipse(frame, ellipse, outer_factor, shape=shape)

        # Set peak to that of the previous source
        source.peak = self.source.peak if self.source is not None else None

        # Set the model to that of the previous source
        if self.model is not None:

            x_min = self.source.x_min
            y_min = self.source.y_min
            x_shift = x_min - source.x_min
            y_shift = y_min - source.y_min
            shifted_model = fitting.shifted_model(self.model, x_shift, y_shift)

            # Set the new model
            source.model = shifted_model

        # Return the new source
        return source

    # -----------------------------------------------------------------

    def remove(self, frame, mask, config, default_fwhm, force=False):

        """
        This function removes the star from a given frame
        :param frame:
        :param mask:
        :param config:
        :param default_fwhm:
        :param force:
        :return:
        """

        # Check which removal method to use, depending on the case
        # (star has model, star has no model but source, star has neither)
        if self.has_model: removal_method = config.method[0]
        elif self.has_source: removal_method = config.method[1]
        else: removal_method = config.method[2]

        # Star that is 'forced' to be removed
        if removal_method is None and force: removal_method = "interpolation"

        # Stars from the DustPedia catalog should always be removed (because we trust this catalog)
        # New: only enable this for optical and NIR (some stars are not present in UV maps and MIR maps)
        if frame.wavelength is None or (frame.wavelength > 0.39 * u.Unit("micron") and frame.wavelength < 10.0 * u.Unit("micron")):
            if self.catalog == "DustPedia" and removal_method is None: removal_method = "interpolation"

        # Remove the star by subtracting the model if a model was found and the method is set to 'model'
        if removal_method == "model":

            # Check whether this star has a model
            if not self.has_model: raise ValueError("Cannot use 'model' mode for stars without a model")

            # Add a new stage to the track record
            if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            self.source = self.source_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Evaluate the model in the cutout of the star's source
            evaluated = self.source.cutout.evaluate_model(self.model)

            # Determine the value at the peak for both the source and the model
            rel_peak = self.source.cutout.rel_position(self.source.peak)

            # Create a box where the model has been subtracted
            subtracted = self.source.cutout - evaluated

            # To plot the difference between the source and the fitted model
            if self.special: plotting.plot_star(self.source.cutout, rel_peak, self.model, "Star about to be removed by subtracting model")

            # Add the evaluated and subtracted boxes to the track record
            if self.has_track_record: self.track_record.append(evaluated)
            if self.has_track_record: self.track_record.append(subtracted)

            # Replace the frame with the subtracted box
            subtracted.replace(frame, where=self.source.mask)

            # Set the subtracted cutout as the background of the source
            self.source.background = subtracted

            # Update the mask
            mask[self.source.cutout.y_slice, self.source.cutout.x_slice] += self.source.mask

        # If a segment was found that can be identified with a source
        elif removal_method == "interpolation":

            # Add a new stage to the track record
            if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            self.source = self.source_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Determine whether we want the background to be sigma-clipped when interpolating over the source
            if self.on_galaxy and config.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = config.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            # NEW: only enable this for optical and NIR (galaxy has smooth emission there but not in UV and MIR)
            # We take 0.39 micron and 20 micron as the limits for 'smoothness'
            if frame.wavelength is None or (frame.wavelength > 0.39 * u.Unit("micron") and frame.wavelength < 10.0 * u.Unit("micron")):
                if self.on_galaxy and config.polynomial_on_galaxy: method = "polynomial"
                else: method = config.interpolation_method
            else: method = config.interpolation_method

            # Estimate the background
            self.source.estimate_background(method, sigma_clip)

            # FOR PLOTTING THE REMOVAL
            if self.special:
                cutout_interpolated = self.source.cutout.copy()
                cutout_interpolated[self.source.mask] = self.source.background[self.source.mask]
                # Do the plotting
                plotting.plot_removal(self.source.cutout, self.source.mask, self.source.background, cutout_interpolated)

            # Add the source to the track record
            if self.has_track_record: self.track_record.append(self.source)

            # Replace the frame with the estimated background
            self.source.background.replace(frame, where=self.source.mask)

            # Update the mask
            mask[self.source.cutout.y_slice, self.source.cutout.x_slice] += self.source.mask

        # None is a valid removal method
        elif removal_method is None: return
        else: raise ValueError("The valid options for removal methods are 'model', 'interpolation' or None")

    # -----------------------------------------------------------------

    def find_saturation(self, frame, config, default_fwhm, star_mask=None):

        """
        This function ...
        :param frame:
        :param config:
        :param default_fwhm:
        :param star_mask:
        :return:
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius for the saturation detection
        model = self.model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Make sure the radius is never smaller than 4 pixels
        radius = max(radius, 4.)

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        ellipse = PixelEllipseRegion(self.pixel_position(frame.wcs), radius, Angle(0.0, "deg"))

        #frame_star_erased = frame.copy()
        #frame_star_erased[self.source.y_slice, self.source.x_slice][self.source.mask] = 0.0

        #saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special)
        #saturation_source = sources.find_source_segmentation(frame_star_erased, ellipse, config, track_record=self.track_record, special=self.special)

        mask_cutout = CutoutMask(self.source.mask, self.source.x_min, self.source.x_max, self.source.y_min, self.source.y_max)
        saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special)

        # Check if the found source segment is larger than the PSF source
        if saturation_source is not None:

            mask_saturation = CutoutMask(saturation_source.mask, saturation_source.x_min, saturation_source.x_max, saturation_source.y_min, saturation_source.y_max)
            mask_saturation_as_cutout = mask_saturation.as_cutout(mask_cutout)
            if self.source.mask.covers(mask_saturation_as_cutout): saturation_source = None

        # If a 'saturation' source was found
        if saturation_source is not None:

            if self.special: log.debug("Initial saturation source found")

            # Calculate the elliptical contour
            # contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)
            x_min = saturation_source.x_min
            x_max = saturation_source.x_max
            y_min = saturation_source.y_min
            y_max = saturation_source.y_max
            contour = sources.find_contour(Box(saturation_source.mask.astype(int), x_min, x_max, y_min, y_max), saturation_source.mask, config.apertures.sigma_level) # determine the segment properties of the actual mask segment

            # Check whether the source centroid matches the star position
            if config.check_centroid:

                if self.special: log.debug("Checking contour parameters ...")

                # Calculate the offset
                difference = contour.center - self.pixel_position(frame.wcs)

                star_mask_cutout = star_mask[saturation_source.cutout.y_slice, saturation_source.cutout.x_slice]

                # Remove the mask of this star from the star_mask_cutout
                x_min_cutout = saturation_source.cutout.x_min
                x_max_cutout = saturation_source.cutout.x_max
                y_min_cutout = saturation_source.cutout.y_min
                y_max_cutout = saturation_source.cutout.y_max

                x_min_source = self.source.cutout.x_min
                x_max_source = self.source.cutout.x_max
                y_min_source = self.source.cutout.y_min
                y_max_source = self.source.cutout.y_max

                try:
                    #plotting.plot_box(star_mask_cutout, title="before removing central source")
                    star_mask_cutout[y_min_source-y_min_cutout:y_max_source-y_min_cutout, x_min_source-x_min_cutout:x_max_source-x_min_cutout][self.source.mask] = False
                    #plotting.plot_box(star_mask_cutout, title="after removing central source")
                except IndexError: pass
                    #plotting.plot_box(frame[saturation_source.y_slice, saturation_source.x_slice])
                    #plotting.plot_box(saturation_source.mask)
                    #plotting.plot_box(star_mask_cutout)
                    #print(star_mask_cutout.shape)
                    #plotting.plot_box(saturation_source.cutout)
                    #print(saturation_source.cutout.shape)
                    #plotting.plot_box(self.source.mask)
                    #print(self.source.mask.shape)
                    #print(y_min_source, y_min_cutout, y_max_source, y_min_cutout, x_min_source, x_min_cutout, x_max_source, x_min_cutout)
                    #print(y_min_source-y_min_cutout) # becomes negative!
                    #print(y_max_source-y_min_cutout)
                    #print(x_min_source-x_min_cutout) # becomes negative !
                    #print(x_max_source-x_min_cutout)
                    #print(star_mask_cutout[y_min_source-y_min_cutout:y_max_source-y_min_cutout, x_min_source-x_min_cutout:x_max_source-x_min_cutout].shape)

                    #source_mask_smaller = self.source.mask[y_min_cutout-y_min_source:,x_min_cutout-x_min_source]

                    #star_mask_cutout[0:y_max_source-y_min_cutout][0:x_max_source-x_min_cutout][source_mask_smaller] = False

                    # TODO: fix this problem ! (how can it be that the source box is not inside the saturation box??)
                    # saturation sources are created by expanding the initial source box ??

                # Discard this saturation source if the centroid offset or the ellipticity is too large
                if not masks.overlap(saturation_source.mask, star_mask_cutout):
                    if self.special: log.debug("Checking offset and ellipticity")
                    if difference.norm > config.max_centroid_offset or contour.ellipticity > config.max_centroid_ellipticity:
                        if self.special: log.debug("Found to large offset or ellipticity: not a saturation source")
                        return
                else:
                    if self.special: log.debug("Saturation mask overlaps other stars, so contour parameters will not be checked")

            if config.second_segmentation:

                # Find all of the saturation light in a second segmentation step
                saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special, sigma_level=config.second_sigma_level)
                #contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)
                contour = sources.find_contour(saturation_source.mask.astype(int), saturation_source.mask, config.apertures.sigma_level) # determine the segment properties of the actual mask segment

                # Check whether the source centroid matches the star position
                if config.check_centroid:

                    # Calculate the offset
                    difference = contour.center - self.pixel_position(frame.wcs)

                    star_mask_cutout = star_mask[saturation_source.cutout.y_slice, saturation_source.cutout.x_slice]

                    # Remove the mask of this star from the star_mask_cutout
                    x_min_cutout = saturation_source.cutout.x_min
                    x_max_cutout = saturation_source.cutout.x_max
                    y_min_cutout = saturation_source.cutout.y_min
                    y_max_cutout = saturation_source.cutout.y_max

                    x_min_source = self.source.cutout.x_min
                    x_max_source = self.source.cutout.x_max
                    y_min_source = self.source.cutout.y_min
                    y_max_source = self.source.cutout.y_max

                    #plotting.plot_box(star_mask_cutout, title="before removing central source")
                    star_mask_cutout[y_min_source-y_min_cutout:y_max_source-y_min_cutout, x_min_source-x_min_cutout:x_max_source-x_min_cutout][self.source.mask] = False
                    #plotting.plot_box(star_mask_cutout, title="after removing central source")

                    # Discard this saturation source if the centroid offset or the ellipticity is too large
                    if not masks.overlap(saturation_source.mask, star_mask_cutout):
                        if difference.norm > config.max_centroid_offset or contour.ellipticity > config.max_centroid_ellipticity: return

            # Replace the pixels of the cutout box by the pixels of the original frame (because the star itself is already removed)
            #saturation_source.cutout = frame.box_like(saturation_source.cutout)

            # TODO: check with classifier to verify this is actually a saturation source!

            if self.special: saturation_source.plot(title="Final saturation source")

            # Replace the source by a source that covers the saturation
            self.saturation = saturation_source
            self.contour = contour

    # -----------------------------------------------------------------

    def remove_saturation(self, frame, mask, config):

        """
        This function ...
        :param frame:
        :param mask:
        :param config:
        """

        # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
        if self.on_galaxy and config.no_sigma_clip_on_galaxy: sigma_clip = False
        else: sigma_clip = config.sigma_clip

        # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
        # NEW: only enable this for optical and IR (galaxy has smooth emission there but not in UV)
        if frame.wavelength is None or frame.wavelength > 0.39 * u.Unit("micron"):
            if self.on_galaxy and config.polynomial_on_galaxy: interpolation_method = "polynomial"
            else: interpolation_method = config.interpolation_method
        else: interpolation_method = config.interpolation_method

        # Estimate the background
        self.saturation.estimate_background(interpolation_method, sigma_clip)

        # FOR PLOTTING THE REMOVAL
        if self.special:
            cutout_interpolated = self.saturation.cutout.copy()
            cutout_interpolated[self.saturation.mask] = self.saturation.background[self.saturation.mask]
            # Do the plotting
            plotting.plot_removal(self.saturation.cutout, self.saturation.mask, self.saturation.background, cutout_interpolated)

        # Replace the frame with the estimated background
        self.saturation.background.replace(frame, where=self.saturation.mask)

        # Update the mask
        mask[self.saturation.cutout.y_slice, self.saturation.cutout.x_slice] += self.saturation.mask

# -----------------------------------------------------------------
