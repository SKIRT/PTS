#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.star Contains the abstract Star class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical units
from astropy import units as u
from astropy.coordinates import Angle

# Import the relevant AstroMagic classes and modules
from .skyobject import SkyObject
from ..core import Source
from ..tools import statistics, fitting, masks, plotting
from ..analysis import sources
from ..basics import Ellipse

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

    def ellipse(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param wcs:
        :param pixelscale:
        :param default_radius:
        :return:
        """

        center, radius, angle = self.ellipse_parameters(wcs, pixelscale, default_radius)
        return Ellipse(center, radius, angle)

    # -----------------------------------------------------------------

    def ellipse_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param wcs:
        :param pixelscale:
        :param initial_radius:
        :return:
        """

        # Return the parameters
        return self.pixel_position(wcs), default_radius, Angle(0.0, u.deg)

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

    def source_at_sigma_level(self, frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=False, shape=None):

        """
        This function ...
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

                print("WARNING: star source does not have peak")
                center = self.pixel_position(frame.wcs)

        else:

            # Calculate the pixel coordinate of the star's position
            center = self.pixel_position(frame.wcs)

        # Create the new source
        ellipse = Ellipse(center, radius, Angle(0.0, u.deg))
        source = Source.from_ellipse(frame, ellipse, outer_factor, shape=shape)

        # Set peak to that of the previous source
        source.peak = self.source.peak if self.source is not None else None

        # Return the new source
        return source

    # -----------------------------------------------------------------

    def remove(self, frame, mask, config, default_fwhm, force=False):

        """
        This function removes the star from a given frame
        :param frame:
        :param mask:
        :return:
        """

        # Check which removal method to use, depending on the case
        # (star has model, star has no model but source, star has neither)
        if self.has_model: removal_method = config.method[0]
        elif self.has_source: removal_method = config.method[1]
        else: removal_method = config.method[2]

        # Star that is 'forced' to be removed
        if removal_method is None and force: removal_method = "interpolation"

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
            if self.on_galaxy and config.polynomial_on_galaxy: method = "polynomial"
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

    def find_saturation(self, frame, original_frame, config, default_fwhm, galaxy_mask=None, star_mask=None):

        """
        This function ...
        :param frame:
        :param original_frame:
        :param config:
        :param default_fwhm:
        :param galaxy_mask:
        :return:
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius for the saturation detection
        model = self.model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        ellipse = Ellipse(self.pixel_position(frame.wcs), radius, Angle(0.0, u.Unit("deg")))
        saturation_source = sources.find_source_segmentation(original_frame, ellipse, config, track_record=self.track_record, special=self.special)
        #saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special)

        # If a 'saturation' source was found
        if saturation_source is not None:

            # Calculate the elliptical contour
            contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)

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

            if config.second_segmentation:

                # Find all of the saturation light in a second segmentation step
                saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special, sigma_level=config.second_sigma_level)
                contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)

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
            saturation_source.cutout = original_frame.box_like(saturation_source.cutout)

            # TODO: check with classifier to verify this is actually a saturation source!

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
        if self.on_galaxy and config.polynomial_on_galaxy: interpolation_method = "polynomial"
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
