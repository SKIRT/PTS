#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.source Contains the Source class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks
from photutils import detect_sources
from photutils import detect_threshold
from photutils import EllipticalAperture
from astropy.coordinates import Angle
from astropy import units as u
from astropy.convolution import convolve, convolve_fft

# Import the relevant AstroMagic classes and modules
from . import Image, Frame
from .box import Box
from ..basics import Position, Extent, Mask, Ellipse
from ..tools import plotting, statistics

# -----------------------------------------------------------------

class Source(object):

    """
    This class...
    """

    def __init__(self, center, radius, angle, factor, cutout, mask, background=None, removed=None, peak=None):

        """
        The constructor ...
        """

        # Set attributes
        self.center = center
        self.radius = radius
        self.angle = angle
        self.factor = factor
        self.cutout = cutout
        self.mask = mask
        self.background = background
        self.removed = removed
        self.peak = peak

        # The aperture
        a = self.radius.x if isinstance(self.radius, Extent) else self.radius
        b = self.radius.y if isinstance(self.radius, Extent) else self.radius
        self.aperture = EllipticalAperture((self.center.x, self.center.y), a, b, theta=self.angle.radian)

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, frame, center, radius, angle, factor, shape=None):

        """
        This function ...
        :param frame:
        :param center:
        :param radius:
        :param angle:
        :param factor:
        :return:
        """

        # Create the cutout box
        cutout = Box.from_ellipse(frame, center, radius * factor, angle, shape=shape)

        # Calculate the relative coordinate of the center for the cutout box
        rel_center = cutout.rel_position(center)

        # Create masks
        #mask = Mask.from_ellipse(cutout.xsize, cutout.ysize, rel_center, radius, angle) # Old way

        ellipse = Ellipse(rel_center, radius, angle)
        mask = Mask.from_ellipse(cutout.xsize, cutout.ysize, ellipse)
        #plotting.plot_difference(mask, mask2)

        # Set (estimated) background and removed to None
        background = None
        removed = None

        # Set peak position to None initially
        peak = None

        # Create and return a new Source instance
        return cls(center, radius, angle, factor, cutout, mask, background, removed, peak)

    # -----------------------------------------------------------------

    @classmethod
    def from_aperture(cls, frame, aperture, factor):

        """
        This function ...
        :param frame:
        :param aperture:
        :return:
        """

        # TODO: make this work with apertures other than EllipticalAperture

        # Get the parameters of the elliptical aperture
        x_center, y_center = aperture.positions[0]
        center = Position(x=x_center, y=y_center)

        major = aperture.a
        minor = aperture.b

        radius = Extent(x=major, y=minor)

        # theta is in radians
        angle = Angle(aperture.theta, u.rad)

        # Return a new Mask object
        return cls.from_ellipse(frame, center, radius, angle, factor)

    # -----------------------------------------------------------------

    @classmethod
    def from_image(cls, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Get the interesting properties from the header
        x_min = int(image.metadata["x_min"])
        x_max = int(image.metadata["x_max"])
        y_min = int(image.metadata["y_min"])
        y_max = int(image.metadata["y_max"])
        center = Position(float(image.metadata["center_x"]), float(image.metadata["center_y"]))
        if "radius" in image.metadata: radius = float(image.metadata["radius"])
        elif "radius_x" in image.metadata: radius = Extent(float(image.metadata["radius_x"]), float(image.metadata["radius_y"]))
        else: RuntimeError("Radius information is missing in metadata")
        angle = Angle(float(image.metadata["angle"]), u.deg)
        factor = float(image.metadata["factor"])
        if "peak_x" in image.metadata: peak = Position(float(image.metadata["peak_x"]), float(image.metadata["peak_y"]))
        else: peak = None

        # Get the cutout and mask
        cutout = Box(image.frames.cutout, x_min, x_max, y_min, y_max)
        mask = Mask(image.frames.mask)

        # Check whether a background frame is present
        if "background" in image.frames: background = Box(image.frames.background, x_min, x_max, y_min, y_max)
        else: background = None

        # Check whether a frame is present that represents the removed source
        if "removed"in image.frames: removed = Box(image.frames.removed, x_min, x_max, y_min, y_max)
        else: removed = None

        # Create and return a new Source instance
        return cls(center, radius, angle, factor, cutout, mask, background, removed, peak)

    # -----------------------------------------------------------------

    @property
    def background_mask(self):

        """
        This function ...
        :return:
        """

        return self.mask.inverse()

    # -----------------------------------------------------------------

    @property
    def subtracted(self):

        """
        This function ...
        :return:
        """

        # Return the box with the background subtracted
        return self.cutout - self.background

    # -----------------------------------------------------------------

    @property
    def has_background(self):

        """
        This function ...
        :return:
        """

        return self.background is not None

    # -----------------------------------------------------------------

    @property
    def has_removed(self):

        """
        This function ...
        :return:
        """

        return self.removed is not None

    # -----------------------------------------------------------------

    @property
    def has_peak(self):

        """
        This function ...
        :return:
        """

        return self.peak is not None

    # -----------------------------------------------------------------

    @property
    def flux(self):

        """
        This function ...
        :return:
        """

        # Calculate the relative coordinate of the center in the background box
        #rel_center = self.background.rel_position(self.center)

        # Method with photutils:
        #aperture = CircularAperture((x_center, y_center), r=x_radius)
        #phot_table = aperture_photometry(self.frames[frame_name], aperture, mask=np.isnan(self.frames[frame_name]))
        #flux = phot_table[0]["aperture_sum"]

        # Sum the pixel values

        #plotting.plot_box(np.ma.masked_array(np.asarray(self.subtracted), mask=self.cutout_mask.inverse()))

        value = np.ma.sum(np.ma.masked_array(np.asarray(self.subtracted), mask=self.background_mask))

        if np.isnan(value): return 0.0
        else: return value

    # -----------------------------------------------------------------

    def estimate_background(self, method, sigma_clip=True, sigma_level=3.0):

        """
        This function ...
        :return:
        """

        # Perform sigma-clipping on the background if requested
        if sigma_clip: mask = statistics.sigma_clip_mask(self.cutout, sigma_level=sigma_level, mask=self.mask)
        else: mask = self.mask

        # Perform the interpolation
        self.background = self.cutout.interpolated(mask, method)

    # -----------------------------------------------------------------

    def find_center_segment(self, sigma_level, kernel=None, min_pixels=5):

        """
        This function ...
        :return:
        """

        # If a subtracted box is present, use it to find the center segment
        box = self.subtracted if self.has_background else self.cutout

        if not np.all(self.mask):

            #print(box)
            #print(self.mask)

            # Calculate threshold for segmentation
            try:
                mean, median, stddev = statistics.sigma_clipped_statistics(box, mask=self.mask)
                threshold = mean + stddev * sigma_level
            except TypeError:

                print(box)
                print(self.mask)

                print("not_nan=", np.sum(np.logical_not(np.isnan(box))))

                exit()

        else:

            threshold = detect_threshold(box, snr=2.0) #snr=2.0

        # Perform the segmentation
        segments = detect_sources(box, threshold, npixels=min_pixels, filter_kernel=kernel).data

        # To plot the multiple segments that are detected
        #if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Get the label of the center segment
        rel_center = self.cutout.rel_position(self.center)
        label = segments[rel_center.y, rel_center.x]

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: return Mask(np.zeros_like(self.cutout, dtype=bool))

        # Create a mask of the center segment
        else: return Mask((segments == label))

    # -----------------------------------------------------------------

    def locate_peaks(self, sigma_level, kernel=None):

        """
        This function ...
        :return:
        """

        # If a subtracted box is present, use it to locate the peaks
        box = self.subtracted if self.has_background else self.cutout

        # Calculate the sigma-clipped statistics of the box
        mean, median, stddev = statistics.sigma_clipped_statistics(box, sigma=3.0, mask=self.background_mask) # Sigma 3.0 for clipping is what photutils uses in detect_threshold
        threshold = median + (sigma_level * stddev)

        # Convolve the box with the given kernel, if any
        if kernel is not None: box = convolve_fft(box, kernel, normalize_kernel=True)

        # Find peaks
        #threshold = detect_threshold(box, snr=2.0) # other method
        peaks = find_peaks(box, threshold, box_size=5, mask=self.background_mask)

        # For some reason, once in a while, an ordinary list comes out of the find_peaks routine instead of an
        # Astropy Table instance. We assume we need an empty table in this case
        if type(peaks) is list: peaks = Table([[], []], names=('x_peak', 'y_peak'))

        # Initialize a list to contain the peak positions
        positions = []

        # Loop over the peaks
        for peak in peaks:

            # Calculate the absolute x and y coordinate of the peak
            x = peak['x_peak'] + self.cutout.x_min
            y = peak['y_peak'] + self.cutout.y_min

            # Add the coordinates to the positions list
            positions.append(Position(x, y))

        # If exactly one peak was found, set the self.peak attribute accordingly
        if len(positions) == 1: self.peak = positions[0]

        # Return the list of peak positions
        return positions

    # -----------------------------------------------------------------

    def zoom(self, factor):

        """
        This function ...
        :return:
        """

        # Create a copy of this object
        source = copy.deepcopy(self)

        # Zoom in on the cutout
        source.cutout = self.cutout.zoom(self.center, factor)

        # Calculate the relative coordinate of the center for the cutout and background boxes
        rel_center = source.cutout.rel_position(self.center)

        # Decrease the radius
        source.radius = self.radius / factor

        # Create smaller mask
        ellipse = Ellipse(rel_center, source.radius, source.angle)
        source.mask = Mask.from_ellipse(source.cutout.xsize, source.cutout.ysize, ellipse)

        # Set other properties to None
        source.background = None
        source.removed = None

        # Return the new source
        return source

    # -----------------------------------------------------------------

    def plot(self, title=None, peaks=None, show=True, scale="sqrt", frame=None):

        """
        This function ...
        :return:
        """

        if peaks is not None:

            x_positions = []
            y_positions = []

            # Loop over all peaks
            for peak in peaks:

                rel_peak = self.cutout.rel_position(peak)
                x_positions.append(rel_peak.x)
                y_positions.append(rel_peak.y)

            peak_coordinates = [x_positions, y_positions]

        elif self.peak is not None:

            rel_peak = self.cutout.rel_position(self.peak)
            peak_coordinates = [[rel_peak.x], [rel_peak.y]]

        else: peak_coordinates = None

        # If the background has been estimated for this source
        if self.has_background:

            # Do the plotting
            plotting.plot_source(self.cutout, self.mask, self.background, peaks=peak_coordinates, title=title, show=show, scale=scale, frame=frame)

        # Else, we just have a background and cutout box
        else:

            # Do the plotting
            plotting.plot_background_center(self.cutout, self.mask, peaks=peak_coordinates, title=title, show=show, scale=scale)

    # -----------------------------------------------------------------

    def save(self, path, frame=None):

        """
        This function ...
        :param path:
        :param frame: replacement frame
        :return:
        """

        cutout = frame.box_like(self.cutout) if frame is not None else self.cutout

        # Create an image to contain the cutout box and source mask
        image = Image()
        image.add_frame(Frame(cutout), "cutout")

        # Add the background-subtracted cutout, if present
        if self.has_background: image.add_frame(Frame(self.background), "background")

        # Add the source-subtracted cutout, if present
        if self.has_removed: image.add_frame(Frame(self.removed), "removed")

        # Add the source mask to the image
        image.add_frame(Frame(self.mask.astype(int)), "mask")

        # Add meta information specifically for sources
        image.metadata["x_min"] = self.cutout.x_min
        image.metadata["x_max"] = self.cutout.x_max
        image.metadata["y_min"] = self.cutout.y_min
        image.metadata["y_max"] = self.cutout.y_max
        image.metadata["center_x"] = float(self.center.x) # is array, why?
        image.metadata["center_y"] = float(self.center.y) # is array, why?
        try: # x and y are sometimes Quantity ??
            image.metadata["radius_x"] = self.radius.x.value
            image.metadata["radius_y"] = self.radius.y.value
        except AttributeError: # x and y are sometimes float ??
            try:
                image.metadata["radius_x"] = self.radius.x
                image.metadata["radius_y"] = self.radius.y
            except AttributeError:  # radius is sometimes float (for circles)
                image.metadata["radius"] = self.radius
        image.metadata["angle"] = self.angle.degree
        image.metadata["factor"] = self.factor

        if self.has_peak:

            image.metadata["peak_x"] = self.peak.x
            image.metadata["peak_y"] = self.peak.y

        # Save the image
        image.frames.select_all()
        image.save(path)

# -----------------------------------------------------------------
