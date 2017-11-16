#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.detection Contains the Detection class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy
from skimage import morphology

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks
from photutils import detect_sources
from photutils import detect_threshold
from astropy.coordinates import Angle
from astropy.convolution import convolve, convolve_fft

# Import the relevant PTS classes and modules
from .image import Image
from .frame import Frame
from .cutout import Cutout
from ..basics.mask import Mask
from ..region.ellipse import PixelEllipseRegion
from ..region.rectangle import PixelRectangleRegion
from ..region.circle import PixelCircleRegion
from ..tools import plotting, statistics
from ..basics.stretch import PixelStretch
from ..basics.coordinate import PixelCoordinate

# -----------------------------------------------------------------

class Detection(object):

    """
    This class...
    """

    def __init__(self, center, radius, angle, factor, cutout, mask, background=None, removed=None, peak=None, shape=None):

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
        self.shape = shape

        self.special = False
        self.contamination = None

        # The elliptical contour
        self.contour = PixelEllipseRegion(self.center, self.radius, self.angle)

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def around_coordinate(cls, frame, coordinate, radius, factor):

        """
        This function ...
        :param frame:
        :param coordinate:
        :param radius:
        :param factor:
        :return:
        """

        circle = PixelCircleRegion(coordinate, radius)
        return cls.from_shape(frame, circle, factor)

    # -----------------------------------------------------------------

    @classmethod
    def from_shape(cls, frame, shape, factor):

        """
        This function ...
        :param frame:
        :param shape:
        :param factor:
        :return:
        """

        # Create a source from the bounding box of the shape
        source = cls.from_rectangle(frame, shape.bounding_box, factor)

        #print(shape.center, shape.radius)
        #print(source.shift)
        #shape_min = shape - source.shift
        #print(shape_min.center, shape_min.radius)

        # Calculate the relative coordinate of the center for the cutout box
        #rel_center = cutout.rel_position(center)

        # Create the shifted shape
        relative_shape = shape - source.shift

        # Create a mask based on the shape, shifted into the source cutout
        mask = relative_shape.to_mask(source.cutout.xsize, source.cutout.ysize)

        #plotting.plot_mask(mask, title="mask")

        # Set the source mask
        source.mask = mask

        # Set the source shape
        source.shape = relative_shape

        # Return the source
        return source

    # -----------------------------------------------------------------

    @classmethod
    def from_rectangle(cls, frame, rectangle, factor):

        """
        This function ...
        :param frame:
        :param rectangle:
        :param factor:
        :return:
        """

        # Get original rectangle properties
        center = rectangle.center
        radius = rectangle.radius
        angle = rectangle.angle

        # Create cutout box
        rectangle = PixelRectangleRegion(center, radius * factor, angle) # new expanded rectangle
        cutout = Cutout.from_rectangle(frame, rectangle)

        # Calculate the relative coordinate of the center for the cutout box
        rel_center = cutout.rel_position(center)

        # Create source mask
        rectangle = PixelRectangleRegion(rel_center, radius, angle)
        mask = Mask.from_shape(rectangle, cutout.xsize, cutout.ysize)

        # Set (estimated) background and removed to None
        background = None
        removed = None

        # Set peak position to None initially
        peak = None

        # Create and return a new Source instance
        return cls(center, radius, angle, factor, cutout, mask, background, removed, peak, rectangle)

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, frame, ellipse, factor, shape=None):

        """
        This function ...
        :param frame:
        :param ellipse:
        :param factor:
        :param shape:
        :return:
        """

        # Get original ellipse properties
        center = ellipse.center
        radius = ellipse.radius
        angle = ellipse.angle

        #print(center, radius.axis1, radius.axis2, angle)

        # Create cutout box
        ellipse = PixelEllipseRegion(center, radius * factor, angle) # new, expanded ellipse
        cutout = Cutout.from_ellipse(frame, ellipse, shape)

        # Calculate the relative coordinate of the center for the cutout box
        rel_center = cutout.rel_position(center)

        # Create source mask
        ellipse = PixelEllipseRegion(rel_center, radius, angle)
        mask = Mask.from_shape(ellipse, cutout.xsize, cutout.ysize)

        # Set (estimated) background and removed to None
        background = None
        removed = None

        # Set peak position to None initially
        peak = None

        # Create and return a new Source instance
        return cls(center, radius, angle, factor, cutout, mask, background, removed, peak, ellipse)

    # -----------------------------------------------------------------

    @classmethod
    def from_mask(cls, frame, mask, factor):

        """
        This function ...
        :param frame:
        :param mask:
        :param factor:
        :return:
        """

        from ..analysis.sources import find_contour_mask

        # Create contour for mask
        ellipse = find_contour_mask(mask)
        #print(ellipse)

        #bounding_box = ellipse.bounding_box
        #x_min = int(round(bounding_box.x_min - 20))
        #x_max = int(round(bounding_box.x_max + 20))
        #y_min = int(round(bounding_box.y_min - 20))
        #y_max = int(round(bounding_box.y_max + 20))
        #ellipse.center.x -= x_min
        #ellipse.center.y -= y_min
        #plotting.plot_box(mask[y_min:y_max, x_min:x_max], region=ellipse)

        # Create from ellipse
        detection = cls.from_ellipse(frame, ellipse, factor=factor)

        #print(type(mask))
        #print(mask)
        #plotting.plot_mask(mask)
        #plotting.plot_mask(mask.data[detection.y_slice, detection.x_slice])

        # Set the mask
        mask_cutout = Mask(mask.data[detection.y_slice, detection.x_slice])
        #plotting.plot_mask(mask_cutout, title="mask cutout")
        detection.mask = mask_cutout

        # Return the detection
        return detection

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
        center = PixelCoordinate(float(image.metadata["center_x"]), float(image.metadata["center_y"]))

        if "radius" in image.metadata: radius = float(image.metadata["radius"])
        elif "radius_x" in image.metadata: radius = PixelStretch(float(image.metadata["radius_x"]), float(image.metadata["radius_y"]))
        else: raise RuntimeError("Radius information is missing in metadata")

        angle = Angle(float(image.metadata["angle"]), "deg")
        factor = float(image.metadata["factor"])
        if "peak_x" in image.metadata: peak = PixelCoordinate(float(image.metadata["peak_x"]), float(image.metadata["peak_y"]))
        else: peak = None

        # Get the cutout and mask
        cutout = Cutout(image.frames["cutout"], x_min, x_max, y_min, y_max)
        mask = Mask(image.frames["mask"])

        # Check whether a background frame is present
        if "background" in image.frames: background = Cutout(image.frames["background"], x_min, x_max, y_min, y_max)
        else: background = None

        # Check whether a frame is present that represents the removed source
        if "removed" in image.frames: removed = Cutout(image.frames["removed"], x_min, x_max, y_min, y_max)
        else: removed = None

        # Create and return a new Source instance
        return cls(center, radius, angle, factor, cutout, mask, background, removed, peak)

    # -----------------------------------------------------------------

    @property
    def x_slice(self):

        """
        This property ...
        :return:
        """

        return self.cutout.x_slice

    # -----------------------------------------------------------------

    @property
    def y_slice(self):

        """
        This property ...
        :return:
        """

        return self.cutout.y_slice

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.cutout.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.cutout.ysize

    # -----------------------------------------------------------------

    @property
    def shift(self):

        """
        This function ...
        :return:
        """

        return PixelStretch(self.x_min, self.y_min)

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This property ...
        :return:
        """

        return self.cutout.x_min

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return self.cutout.x_max

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return self.cutout.y_min

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return self.cutout.y_max

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

    def get_flux(self, without_background=False):

        """
        This function ...
        :param without_background:
        :return:
        """

        if without_background: data = self.cutout
        else: data = self.subtracted

        value = np.ma.sum(np.ma.masked_array(np.asarray(data), mask=self.background_mask))

        if np.isnan(value): return 0.0
        else: return value

    # -----------------------------------------------------------------

    def estimate_background(self, method, sigma_clip=True, sigma_level=3.0, sigma=None, fwhm=None, max_iterations=10,
                            not_converge="keep", plot_interpolation=False):

        """
        This function ...
        :param method:
        :param sigma_clip:
        :param sigma_level:
        :param sigma:
        :param fwhm:
        :param max_iterations:
        :param not_converge:
        :param plot_interpolation:
        :return:
        """

        # Make a distinction between the "PTS" way of estimating the background and all other methods.
        # For the PTS way, in case sigma clipping is disabled, this means it is disabled only for the 'polynomial fitting step'
        # of the background estimation, so provide two distinct masks to the interpolated() method: the clipped mask for
        # the 'background noise' estimation and the non-clipped mask for the 'polynomial fitting step'.
        if method == "pts":

            if sigma_clip:
                try:
                    if self.mask is not None and np.sum(~self.mask) == 0: raise TypeError('Invalid mask')
                    mask = statistics.sigma_clip_mask(self.cutout, sigma_level=sigma_level, mask=self.mask)
                except TypeError:
                    #plotting.plot_box(self.cutout)
                    #plotting.plot_mask(self.mask)
                    #print("xsize", self.cutout.xsize, self.cutout.ysize)
                    radius = int(round(0.25 * self.cutout.xsize))
                    #print("radius", 0.25*self.cutout.xsize, radius)
                    disk = morphology.disk(radius, dtype=bool)
                    mask = Mask.empty_like(self.cutout)
                    x_min = int(round(0.5 * (self.cutout.xsize - disk.shape[1])))
                    y_min = int(round(0.5 * (self.cutout.ysize - disk.shape[0])))
                    #plotting.plot_mask(mask)
                    mask[y_min:y_min+disk.shape[0], x_min:x_min+disk.shape[1]] = disk
                    #plotting.plot_mask(mask)

                no_clip_mask = None
            else:
                if self.mask is not None and np.sum(~self.mask) == 0:

                    radius = int(round(0.50 * self.cutout.xsize))
                    disk = morphology.disk(radius, dtype=bool)
                    mask = Mask.empty_like(self.cutout)
                    
                    x_min = int(round(0.5 * (self.cutout.xsize - disk.shape[1])))
                    y_min = int(round(0.5 * (self.cutout.ysize - disk.shape[0])))
                    mask[y_min:y_min+disk.shape[0], x_min:x_min+disk.shape[1]] = disk
                    self.mask = mask

                else: mask = statistics.sigma_clip_mask(self.cutout, sigma_level=sigma_level, mask=self.mask)

                no_clip_mask = self.mask

        # Other methods
        else:

            # Perform sigma-clipping on the background if requested
            if sigma_clip:
                try: mask = statistics.sigma_clip_mask(self.cutout, sigma_level=sigma_level, mask=self.mask)
                except: mask = self.mask
            else: mask = self.mask
            no_clip_mask = None

        # Add contamination
        if self.contamination is not None:
            mask = mask + self.contamination
            if no_clip_mask is not None: no_clip_mask = no_clip_mask + self.contamination

        # Plot mask
        if self.special: plotting.plot_mask(mask, title="interpolation mask")

        # Perform the interpolation
        if self.special: plot_interpolation = True
        self.background = self.cutout.interpolated(mask, method, no_clip_mask=no_clip_mask, plot=plot_interpolation,
                                                   sigma=sigma, fwhm=fwhm, max_iterations=max_iterations,
                                                   not_converge=not_converge)

        # Special: plot
        if self.special: self.plot(title="background estimated")

        # Return the interpolation mask
        return mask

    # -----------------------------------------------------------------

    def find_center_segment(self, sigma_level, kernel=None, min_pixels=5):

        """
        This function ...
        :param sigma_level:
        :param kernel:
        :param min_pixels:
        :return:
        """

        # If a subtracted box is present, use it to find the center segment
        box = self.subtracted if self.has_background else self.cutout

        if not np.all(self.mask):

            mean, median, stddev = statistics.sigma_clipped_statistics(box, mask=self.mask.data)
            threshold = mean + stddev * sigma_level

        else: threshold = detect_threshold(box, snr=2.0) #snr=2.0

        # Perform the segmentation
        segments = detect_sources(box, threshold, npixels=min_pixels, filter_kernel=kernel).data

        # To plot the multiple segments that are detected
        #if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Get the label of the center segment
        rel_center = self.cutout.rel_position(self.center)
        label = segments[int(round(rel_center.y)), int(round(rel_center.x))]

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: return Mask(np.zeros_like(self.cutout, dtype=bool))

        # Create a mask of the center segment
        else: return Mask((segments == label))

    # -----------------------------------------------------------------

    def locate_peaks(self, sigma_level, kernel=None):

        """
        This function ...
        :param sigma_level:
        :param kernel:
        :return:
        """

        # If a subtracted box is present, use it to locate the peaks
        box = self.subtracted if self.has_background else self.cutout

        # Calculate the sigma-clipped statistics of the box
        mean, median, stddev = statistics.sigma_clipped_statistics(box, sigma=3.0, mask=self.background_mask.data) # Sigma 3.0 for clipping is what photutils uses in detect_threshold
        #sigma_level = 1.5   # I once tried to investigate why some clear peaks were not detected, did not have time ..
        threshold = median + (sigma_level * stddev)

        # Convolve the box with the given kernel, if any
        if kernel is not None: box = convolve_fft(box, kernel, normalize_kernel=True)

        # Find peaks
        #threshold = detect_threshold(box, snr=2.0) # other method (snr corresponds to sigma_level as I use it above)
        peaks = find_peaks(box, threshold, box_size=5, mask=self.background_mask)

        # For some reason, once in a while, an ordinary list comes out of the find_peaks routine instead of an
        # Astropy Table instance. We assume we need an empty table in this case
        if type(peaks) is list: peaks = Table([[], []], names=('x_peak', 'y_peak'))

        # Initialize a list to contain the peak positions
        positions = []

        # Loop over the peaks
        for peak in peaks:

            # Calculate the absolute x and y coordinate of the peak
            x_rel = peak['x_peak']
            y_rel = peak['y_peak']
            x = x_rel + self.cutout.x_min
            y = y_rel + self.cutout.y_min

            # Check whether the peak position falls in the box and then add it to the list
            peak_position = PixelCoordinate(x, y)
            if self.cutout.contains(peak_position):
                # Add the coordinates to the positions list
                positions.append(peak_position)
            else: print("DEBUG: peak position", peak_position, "falls outside of box with shape", box.shape)

        # If exactly one peak was found, set the self.peak attribute accordingly
        if len(positions) == 1: self.peak = positions[0]

        # Return the list of peak positions
        return positions

    # -----------------------------------------------------------------

    def zoom(self, factor, min_xpixels=None, min_ypixels=None, min_pixels=None):

        """
        This function ...
        :param factor:
        :param min_xpixels:
        :param min_ypixels:
        :param min_pixels:
        :return:
        """

        # Create a copy of this object
        source = copy.deepcopy(self)

        # Zoom in on the cutout
        source.cutout = self.cutout.zoom(self.center, factor, min_xpixels=min_xpixels, min_ypixels=min_ypixels, min_pixels=min_pixels)

        # Calculate the relative coordinate of the center for the cutout and background boxes
        rel_center = source.cutout.rel_position(self.center)

        # Decrease the radius
        try: ## SOMETIMES, THIS FAILS (WHY??)
            source.radius = self.radius / factor
        except TypeError:
            #print(type(self.radius), type(factor))
            #print(type(self.radius.x), type(self.radius.y))
            source.radius = self.radius * (1.0 / factor)

        # Create smaller mask
        ellipse = PixelEllipseRegion(rel_center, source.radius, source.angle)
        source.mask = Mask.from_shape(ellipse, source.cutout.xsize, source.cutout.ysize)

        # Set other properties to None
        source.background = None
        source.removed = None

        # Return the new source
        return source

    # -----------------------------------------------------------------

    def zoom_out(self, factor, frame, keep_original_mask=False):

        """
        This function ...
        :param factor:
        :param frame:
        :param keep_original_mask:
        :return:
        """

        new_radius = self.radius * factor

        # Create cutout box
        ellipse = PixelEllipseRegion(self.center, new_radius * self.factor, self.angle) # new, expanded ellipse
        cutout = Cutout.from_ellipse(frame, ellipse)

        # Calculate the relative coordinate of the center for the cutout box
        rel_center = cutout.rel_position(self.center)

        # Create source mask
        ellipse = PixelEllipseRegion(rel_center, new_radius, self.angle)
        mask = Mask.from_shape(ellipse, cutout.xsize, cutout.ysize)

        # Create the new source
        new_source = Detection(self.center, new_radius, self.angle, self.factor, cutout, mask, peak=self.peak)

        if keep_original_mask:

            original_x_min = self.cutout.x_min
            original_y_min = self.cutout.y_min
            original_x_max = self.cutout.x_max
            original_y_max = self.cutout.y_max

            new_source.mask = Mask(np.zeros(new_source.cutout.shape))

            rel_x_min = original_x_min - new_source.cutout.x_min
            rel_y_min = original_y_min - new_source.cutout.y_min
            rel_x_max = original_x_max - new_source.cutout.x_min
            rel_y_max = original_y_max - new_source.cutout.y_min

            # Replace source's mask by found center segment mask
            new_source.mask[rel_y_min:rel_y_max, rel_x_min:rel_x_max] = self.mask

        # Return the zoomed-out source
        return new_source

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # If the position does not lie within the cutout box of the galaxy's source, return False
        if not self.cutout.contains(position): return False

        # If it does, check whether the pixel position is masked by the mask of the galaxy's source
        return self.mask.masks(self.cutout.rel_position(position))

    # -----------------------------------------------------------------

    def plot(self, title=None, peaks=None, show=True, scale="sqrt", frame=None):

        """
        This function ...
        :param title:
        :param peaks:
        :param show:
        :param scale:
        :param frame:
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

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path exists
        if self.path is None: raise RuntimeError("Path is not defined for this source")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path, frame=None, origin=None):

        """
        This function ...
        :param path:
        :param frame: replacement frame
        :param origin:
        :return:
        """

        cutout = frame.box_like(self.cutout) if frame is not None else self.cutout

        # Create an image to contain the cutout box and source mask
        image = Image("source")
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
        image.saveto(path, origin=origin)

        # Update the path
        self.path = path

# -----------------------------------------------------------------
