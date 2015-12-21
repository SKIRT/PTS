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
from ..tools import statistics
from ..tools import fitting
from ..analysis import sources

# -----------------------------------------------------------------

class Star(SkyObject):

    """
    This class ...
    """

    def __init__(self, catalog=None, id=None, position=None, ra_error=None, dec_error=None, magnitudes=None, magnitude_errors=None, on_galaxy=False):

        """
        The constructor ...
        :return:
        """

        # Set the attributes
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

        # Set the has_saturation flag to False initially
        self.has_saturation = False

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

    def fit_model(self, config, source=None):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("fitting")

        # Fit model to the source, in a loop over different analytical forms for the model
        for level in range(len(config.model_names)):

            if source is None:

                # Do the fitting
                source, model = sources.fit_model_to_source(self.source, config, self.track_record, level=level)

            else:

                source, model = sources.fit_model_to_source(source, config, self.track_record, level=level)

            # If a model was found, set the attributes of the star object and exit the loop
            if model is not None:

                self.source = source
                self.model = model
                break

    # -----------------------------------------------------------------

    def source_at_sigma_level(self, frame, default_fwhm, sigma_level, outer_factor):

        """
        This function ...
        :return:
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius of the contour in which the star will be removed
        radius = fitting.sigma(self.model) * sigma_level if self.model is not None else default_sigma * sigma_level

        # Determine the center position of the source (center of model if present, otherwise position of the star)
        if self.source is not None:

            # If the star has been modeled succesfully, use the center position of the model
            # Otherwise, use the source's peak
            if self.model is not None: center = fitting.center(self.model)
            else: center = self.source.peak

        else:

            # Calculate the pixel coordinate of the star's position
            center = self.pixel_position(frame.wcs)

        # Create a source and return it
        return Source(frame, center, radius, Angle(0.0, u.deg), outer_factor)

    # -----------------------------------------------------------------

    def remove(self, frame, mask, config, default_fwhm):

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

        # Remove the star by subtracting the model if a model was found and the method is set to 'model'
        if removal_method == "model":

            # Check whether this star has a model
            if not self.has_model: raise ValueError("Cannot use 'model' mode for stars without a model")

            # Add a new stage to the track record
            if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            source = self.source_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Evaluate the model in the cutout of the star's source
            evaluated = source.cutout.evaluate_model(self.model)

            # Determine the value at the peak for both the source and the model
            rel_peak = source.cutout.rel_position(self.source.peak)

            # Create a box where the model has been subtracted
            subtracted = source.cutout - evaluated

            # To plot the difference between the source and the fitted model
            #from ..tools import plotting
            #plotting.plot_star(source.cutout, rel_peak, self.model, "Star about to be removed by subtracting model")

            # Add the evaluated and subtracted boxes to the track record
            if self.has_track_record: self.track_record.append(evaluated)
            if self.has_track_record: self.track_record.append(subtracted)

            # Replace the frame with the subtracted box
            subtracted.replace(frame, where=source.mask)

            # Update the mask
            mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

        # If a segment was found that can be identified with a source
        elif removal_method == "interpolation":

            # Add a new stage to the track record
            if self.has_track_record: self.track_record.set_stage("removal")

            # Create a source for the desired sigma level and outer factor
            source = self.source_at_sigma_level(frame, default_fwhm, config.sigma_level, config.outer_factor)

            # Determine whether we want the background to be sigma-clipped when interpolating over the source
            if self.on_galaxy and config.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = config.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            if self.on_galaxy and config.polynomial_on_galaxy: method = "polynomial"
            else: method = config.interpolation_method

            # Estimate the background
            source.estimate_background(method, sigma_clip)

            # FOR PLOTTING THE REMOVAL
            #import copy
            #cutout_interpolated = copy.deepcopy(source.cutout)
            #cutout_interpolated[source.mask] = source.background[source.mask]
            #from ..tools import plotting
            # Do the plotting
            #plotting.plot_removal(source.cutout, source.mask, source.background, cutout_interpolated)

            # Add the source to the track record
            if self.has_track_record: self.track_record.append(source)

            # Replace the frame with the estimated background
            source.background.replace(frame, where=source.mask)

            # Update the mask
            mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

        # None is a valid removal method
        elif removal_method is None: return
        else: raise ValueError("The valid options for removal methods are 'model', 'interpolation' or None")

    # -----------------------------------------------------------------

    def remove_saturation(self, frame, mask, config, default_fwhm):

        """
        This function ...
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius for the saturation detection
        model = self.model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        source = sources.find_source_segmentation(frame, self.pixel_position(frame.wcs), radius, Angle(0.0, u.deg), config, track_record=self.track_record, special=self.special)

        # If a 'saturation' source was found
        if source is not None:

            # Check whether the source centroid matches the star position
            if config.check_centroid:

                from photutils import segment_properties, properties_table
                from ..basics import Position

                # Get the segment properties
                # Since there is only one segment in the source.mask (the center segment), the props
                # list contains only one entry (one galaxy)
                props = segment_properties(source.cutout, source.mask)
                properties = props[0]

                x_shift = source.cutout.x_min
                y_shift = source.cutout.y_min

                # Obtain the position, orientation and extent
                position = Position(properties.xcentroid.value + x_shift, properties.ycentroid.value + y_shift)
                a = properties.semimajor_axis_sigma.value
                b = properties.semiminor_axis_sigma.value
                theta = properties.orientation.value

                # Create the aperture
                #self.aperture = EllipticalAperture(position, a, b, theta=theta)

                difference = position - self.pixel_position(frame.wcs)

                with open(config.centroid_table_path, 'a') as centroid_file:
                    centroid_file.write(str(difference.norm) + "  " + str(a/b) + "\n")

            # Replace the source by a source that covers the saturation
            self.source = source

            # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
            if self.on_galaxy and config.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = config.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            if self.on_galaxy and config.polynomial_on_galaxy: interpolation_method = "polynomial"
            else: interpolation_method = config.interpolation_method

            # Estimate the background
            self.source.estimate_background(interpolation_method, sigma_clip)

            # FOR PLOTTING THE REMOVAL
            #import copy
            #cutout_interpolated = copy.deepcopy(source.cutout)
            #cutout_interpolated[source.mask] = source.background[source.mask]
            #from ..tools import plotting
            # Do the plotting
            #plotting.plot_removal(source.cutout, source.mask, source.background, cutout_interpolated)

            # Replace the frame with the estimated background
            self.source.background.replace(frame, where=self.source.mask)

            # Update the mask
            mask[self.source.cutout.y_slice, self.source.cutout.x_slice] += self.source.mask

            # Indicate that a saturation source was found
            return True

        # Otherwise, return False
        else: return False

# -----------------------------------------------------------------
