#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.pointsource Contains the PointSource class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .source import Source
from .detection import Detection
from ..core.cutout import CutoutMask, Cutout
from ..tools import statistics, fitting, masks, plotting
from ..analysis import sources
from ..region.ellipse import PixelEllipseRegion
from ...core.basics.log import log
from ..basics.stretch import PixelStretch
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

class PointSource(Source):

    """
    This class...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSource, self).__init__(**kwargs)

        # Set other properties
        self.catalog = kwargs.pop("catalog", None)
        self.id = kwargs.pop("id", None)
        self.ra_error = kwargs.pop("ra_error", None)
        self.dec_error = kwargs.pop("dec_error", None)
        self.confidence = kwargs.pop("confidence", None)
        self.magnitudes = kwargs.pop("magnitudes", dict())
        self.magnitude_errors = kwargs.pop("magnitude_errors", dict())
        self.on_galaxy = kwargs.pop("on_galaxy", False)

        # PSF model
        self.psf_model = None

        # Saturation detection
        self.saturation = None

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.psf_model is not None

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        return fitting.fwhm(self.psf_model) if self.has_model else None

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
    def flux(self):

        """
        This function ...
        :return:
        """

        # Return the flux of the source
        return self.detection.flux

    # -----------------------------------------------------------------

    def get_flux(self, without_background=False):

        """
        This function ...
        :param without_background:
        :return:
        """

        return self.detection.get_flux(without_background)

    # -----------------------------------------------------------------

    def find_contour(self, frame, config, saturation=False):

        """
        This function ...
        :param frame:
        :param config:
        :param saturation:
        :return:
        """

        # Determine which box and mask
        if saturation:

            box = self.saturation.cutout
            mask = self.saturation.mask

            # Implementation
            self._find_contour_impl(frame.wcs, box, mask, config)

        # Call the base class implementation
        else: super(PointSource, self).find_contour(frame, config)

    # -----------------------------------------------------------------

    def detection_from_shape(self, frame, shape, outer_factor):

        """
        This function ...
        :param frame:
        :param shape:
        :param outer_factor:
        :return:
        """

        # Create the detection
        self.detection = Detection.from_shape(frame, shape, outer_factor)

    # -----------------------------------------------------------------

    def detection_at_sigma_level(self, frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=False, shape=None):

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
        if self.psf_model is None or use_default_fwhm: radius = default_sigma * sigma_level
        else: radius = fitting.sigma(self.psf_model) * sigma_level

        # Determine the center position of the detection (center of model if present, otherwise position of the star)
        if self.detection is not None:

            # If the star has been modeled succesfully, use the center position of the model
            # Otherwise, use the source's peak
            if self.psf_model is not None: center = fitting.center(self.psf_model)
            elif self.detection.has_peak: center = self.detection.peak
            else:

                log.warning("Star source does not have peak")
                center = self.pixel_position(frame.wcs)

        # Calculate the pixel coordinate of the star's position
        else: center = self.pixel_position(frame.wcs)

        # Create the new source
        radius = PixelStretch(radius, radius)
        ellipse = PixelEllipseRegion(center, radius)
        detection = Detection.from_ellipse(frame, ellipse, outer_factor, shape=shape)

        # Set peak to that of the previous source
        detection.peak = self.detection.peak if self.detection is not None else None

        # Set the model to that of the previous source
        if self.psf_model is not None:

            x_min = self.detection.x_min
            y_min = self.detection.y_min
            x_shift = x_min - detection.x_min
            y_shift = y_min - detection.y_min
            shifted_model = fitting.shifted_model(self.psf_model, x_shift, y_shift)

            # Set the new model
            detection.model = shifted_model

        # Return the new detection
        return detection

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
        return self.pixel_position(wcs), default_radius, Angle(0.0, "deg")

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
        elif self.has_detection: removal_method = config.method[1]
        else: removal_method = config.method[2]

        # Star that is 'forced' to be removed
        if removal_method is None and force: removal_method = "interpolation"

        # Stars from the DustPedia catalog should always be removed (because we trust this catalog)
        # New: only enable this for optical and NIR (some stars are not present in UV maps and MIR maps)
        if frame.wavelength is None or (frame.wavelength > 0.39 * u("micron") and frame.wavelength < 10.0 * u("micron")):
            if self.catalog == "DustPedia" and removal_method is None: removal_method = "interpolation"

        # Remove the star by subtracting the model if a model was found and the method is set to 'model'
        if removal_method == "model":

            # Check whether this star has a model
            if not self.has_model: raise ValueError("Cannot use 'model' mode for stars without a model")

            # Add a new stage to the track record
            #if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            self.detection = self.detection_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Evaluate the model in the cutout of the star's source
            evaluated = self.detection.cutout.evaluate_model(self.psf_model)

            # Determine the value at the peak for both the source and the model
            rel_peak = self.detection.cutout.rel_position(self.detection.peak)

            # Create a box where the model has been subtracted
            subtracted = self.detection.cutout - evaluated

            # To plot the difference between the source and the fitted model
            if self.special: plotting.plot_star(self.detection.cutout, rel_peak, self.psf_model, "Star about to be removed by subtracting model")

            # Add the evaluated and subtracted boxes to the track record
            #if self.has_track_record: self.track_record.append(evaluated)
            #if self.has_track_record: self.track_record.append(subtracted)

            # Replace the frame with the subtracted box
            subtracted.replace(frame, where=self.detection.mask)

            # Set the subtracted cutout as the background of the source
            self.detection.background = subtracted

            # Update the mask
            mask[self.detection.cutout.y_slice, self.detection.cutout.x_slice] += self.detection.mask

        # If a segment was found that can be identified with a source
        elif removal_method == "interpolation":

            # Add a new stage to the track record
            #if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            self.detection = self.detection_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Determine whether we want the background to be sigma-clipped when interpolating over the source
            if self.on_galaxy and config.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = config.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            # NEW: only enable this for optical and NIR (galaxy has smooth emission there but not in UV and MIR)
            # We take 0.39 micron and 20 micron as the limits for 'smoothness'
            if frame.wavelength is None or (frame.wavelength > 0.39 * u("micron") and frame.wavelength < 10.0 * u("micron")):
                if self.on_galaxy and config.polynomial_on_galaxy: method = "polynomial"
                else: method = config.interpolation_method
            else: method = config.interpolation_method

            # Estimate the background
            self.detection.estimate_background(method, sigma_clip)

            # FOR PLOTTING THE REMOVAL
            if self.special:

                cutout_interpolated = self.detection.cutout.copy()
                cutout_interpolated[self.detection.mask] = self.detection.background[self.detection.mask]

                # Do the plotting
                plotting.plot_removal(self.detection.cutout, self.detection.mask, self.detection.background, cutout_interpolated)

            # Add the source to the track record
            #if self.has_track_record: self.track_record.append(self.source)

            # Replace the frame with the estimated background
            self.detection.background.replace(frame, where=self.detection.mask)

            # Update the mask
            mask[self.detection.cutout.y_slice, self.detection.cutout.x_slice] += self.detection.mask

        # None is a valid removal method
        elif removal_method is None: return
        else: raise ValueError("The valid options for removal methods are 'model', 'interpolation' or None")

    # -----------------------------------------------------------------

    def fit_model(self, config, detection=None):

        """
        This function ...
        :param config:
        :param detection:
        :param debug:
        :return:
        """

        # Add a new stage to the track record
        #if self.has_track_record: self.track_record.set_stage("fitting")
        track_record = None

        # Fit model to the source, in a loop over different analytical forms for the model
        for level in range(len(config.model_names)):

            # Do the fitting
            if detection is None: detection, model = sources.fit_model_to_source(self.detection, config, track_record, level=level, special=self.special)
            else: detection, model = sources.fit_model_to_source(detection, config, track_record, level=level)

            # If a model was found, set the attributes of the star object and exit the loop
            if model is not None:

                self.detection = detection
                self.psf_model = model
                break

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
        model = self.psf_model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Make sure the radius is never smaller than 4 pixels
        radius = max(radius, 4.)

        # Add a new stage to the track record
        #if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        radius_ellipse = PixelStretch(radius, radius)
        ellipse = PixelEllipseRegion(self.pixel_position(frame.wcs), radius_ellipse)

        # frame_star_erased = frame.copy()
        # frame_star_erased[self.source.y_slice, self.source.x_slice][self.source.mask] = 0.0

        # saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=self.track_record, special=self.special)
        # saturation_source = sources.find_source_segmentation(frame_star_erased, ellipse, config, track_record=self.track_record, special=self.special)

        mask_cutout = CutoutMask(self.detection.mask, self.detection.x_min, self.detection.x_max, self.detection.y_min, self.detection.y_max)
        track_record = None
        saturation_source = sources.find_source_segmentation(frame, ellipse, config, track_record=track_record, special=self.special)

        # Check if the found source segment is larger than the PSF source
        if saturation_source is not None:

            mask_saturation = CutoutMask(saturation_source.mask, saturation_source.x_min,
                                         saturation_source.x_max, saturation_source.y_min,
                                         saturation_source.y_max)
            mask_saturation_as_cutout = mask_saturation.as_cutout(mask_cutout)
            if self.detection.mask.covers(mask_saturation_as_cutout): saturation_source = None

        # If a 'saturation' source was found
        if saturation_source is not None:

            if self.special: log.debug("Initial saturation source found")

            x_min = saturation_source.x_min
            x_max = saturation_source.x_max
            y_min = saturation_source.y_min
            y_max = saturation_source.y_max

            # DEBLEND FIRST
            if config.deblend:

                import numpy as np
                from photutils.segmentation import deblend_sources

                # from astropy.convolution import Kernel2D
                # Kernel2D._model = self.psf_model
                # if self.psf_model is not None:
                #     kernelsize = 2 * int(round(fitting.sigma(self.psf_model) * 3.))
                #     print("kernelsize", kernelsize)
                #     kernel = Kernel2D(x_size=kernelsize)
                # else: kernel = None
                kernel = None

                segments = deblend_sources(saturation_source.cutout, saturation_source.mask.astype(int),
                                           npixels=config.deblending.min_npixels, contrast=config.deblending.contrast,
                                           mode=config.deblending.mode, nlevels=config.deblending.nlevels,
                                           filter_kernel=kernel)

                plotting.plot_box(segments)

                smallest_distance = None
                smallest_distance_mask = None
                for index in np.unique(segments)[1:]:

                    where = segments == index
                    fake_box = Cutout(where.astype(int), x_min, x_max, y_min, y_max)
                    contour = sources.find_contour(fake_box, where, sigma_level=1)

                    difference = contour.center - self.pixel_position(frame.wcs)
                    distance = difference.norm

                    if smallest_distance is None or distance < smallest_distance:
                        smallest_distance = distance
                        smallest_distance_mask = where

                        # print(index, difference.norm)

                # SET NEW MASK
                saturation_source.mask = smallest_distance_mask

            # AFTER DEBLENDING, CALCULATE CONTOUR
            # Calculate the elliptical contour
            # contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)
            contour = sources.find_contour(Cutout(saturation_source.mask.astype(int), x_min, x_max, y_min, y_max),
                                           saturation_source.mask,
                                           config.apertures.sigma_level)  # determine the segment properties of the actual mask segment

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

                x_min_source = self.detection.cutout.x_min
                x_max_source = self.detection.cutout.x_max
                y_min_source = self.detection.cutout.y_min
                y_max_source = self.detection.cutout.y_max

                try:
                    # plotting.plot_box(star_mask_cutout, title="before removing central source")
                    star_mask_cutout[y_min_source - y_min_cutout:y_max_source - y_min_cutout,
                    x_min_source - x_min_cutout:x_max_source - x_min_cutout][self.detection.mask] = False
                    # plotting.plot_box(star_mask_cutout, title="after removing central source")
                except IndexError:
                    pass
                # plotting.plot_box(frame[saturation_source.y_slice, saturation_source.x_slice])
                # plotting.plot_box(saturation_source.mask)
                # plotting.plot_box(star_mask_cutout)
                # print(star_mask_cutout.shape)
                # plotting.plot_box(saturation_source.cutout)
                # print(saturation_source.cutout.shape)
                # plotting.plot_box(self.source.mask)
                # print(self.source.mask.shape)
                # print(y_min_source, y_min_cutout, y_max_source, y_min_cutout, x_min_source, x_min_cutout, x_max_source, x_min_cutout)
                # print(y_min_source-y_min_cutout) # becomes negative!
                # print(y_max_source-y_min_cutout)
                # print(x_min_source-x_min_cutout) # becomes negative !
                # print(x_max_source-x_min_cutout)
                # print(star_mask_cutout[y_min_source-y_min_cutout:y_max_source-y_min_cutout, x_min_source-x_min_cutout:x_max_source-x_min_cutout].shape)

                # source_mask_smaller = self.source.mask[y_min_cutout-y_min_source:,x_min_cutout-x_min_source]

                # star_mask_cutout[0:y_max_source-y_min_cutout][0:x_max_source-x_min_cutout][source_mask_smaller] = False

                # TODO: fix this problem ! (how can it be that the source box is not inside the saturation box??)
                # saturation sources are created by expanding the initial source box ??

                # Discard this saturation source if the centroid offset or the ellipticity is too large
                if not masks.overlap(saturation_source.mask, star_mask_cutout):
                    if self.special: log.debug("Checking offset and ellipticity")
                    if difference.norm > config.max_centroid_offset or contour.ellipticity > config.max_centroid_ellipticity:
                        if self.special: log.debug(
                            "Found to large offset or ellipticity: not a saturation source")
                        return
                else:
                    if self.special: log.debug("Saturation mask overlaps other stars, so contour parameters will not be checked")

            if config.second_segmentation:

                # Find all of the saturation light in a second segmentation step
                track_record = None
                saturation_source = sources.find_source_segmentation(frame, ellipse, config,
                                                                     track_record=track_record,
                                                                     special=self.special,
                                                                     sigma_level=config.second_sigma_level)
                # contour = sources.find_contour(saturation_source.cutout, saturation_source.mask, config.apertures.sigma_level)
                contour = sources.find_contour(saturation_source.mask.astype(int), saturation_source.mask, config.apertures.sigma_level)  # determine the segment properties of the actual mask segment

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

                    x_min_source = self.detection.cutout.x_min
                    x_max_source = self.detection.cutout.x_max
                    y_min_source = self.detection.cutout.y_min
                    y_max_source = self.detection.cutout.y_max

                    # plotting.plot_box(star_mask_cutout, title="before removing central source")
                    star_mask_cutout[y_min_source - y_min_cutout:y_max_source - y_min_cutout,
                    x_min_source - x_min_cutout:x_max_source - x_min_cutout][self.detection.mask] = False
                    # plotting.plot_box(star_mask_cutout, title="after removing central source")

                    # Discard this saturation source if the centroid offset or the ellipticity is too large
                    if not masks.overlap(saturation_source.mask, star_mask_cutout):
                        if difference.norm > config.max_centroid_offset or contour.ellipticity > config.max_centroid_ellipticity: return

            # Replace the pixels of the cutout box by the pixels of the original frame (because the star itself is already removed)
            # saturation_source.cutout = frame.box_like(saturation_source.cutout)

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
        if frame.wavelength is None or frame.wavelength > 0.39 * u("micron"):
            if self.on_galaxy and config.polynomial_on_galaxy:
                interpolation_method = "polynomial"
            else: interpolation_method = config.interpolation_method
        else: interpolation_method = config.interpolation_method

        # Estimate the background
        self.saturation.estimate_background(interpolation_method, sigma_clip)

        # FOR PLOTTING THE REMOVAL
        if self.special:

            cutout_interpolated = self.saturation.cutout.copy()
            cutout_interpolated[self.saturation.mask] = self.saturation.background[self.saturation.mask]

            # Do the plotting
            plotting.plot_removal(self.saturation.cutout, self.saturation.mask, self.saturation.background,
                                  cutout_interpolated)

        # Replace the frame with the estimated background
        self.saturation.background.replace(frame, where=self.saturation.mask)

        # Update the mask
        mask[self.saturation.cutout.y_slice, self.saturation.cutout.x_slice] += self.saturation.mask

# -----------------------------------------------------------------
