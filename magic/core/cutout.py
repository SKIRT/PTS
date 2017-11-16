#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.cutout Contains the Cutout and CutoutMask classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from scipy import ndimage

# Import the relevant PTS classes and modules
from ..basics.vector import Position, Extent
from ..region.rectangle import PixelRectangleRegion
from ..tools import cropping, fitting, interpolation, plotting
from ...core.basics.log import log
from ..basics.coordinate import PixelCoordinate
from .mask import Mask
from .alpha import AlphaMask
from ..basics.mask import Mask as oldMask
from ..basics.mask import MaskBase

# -----------------------------------------------------------------

interpolation_methods = ["polynomial", "local_mean", "idw", "mean", "median", "biharmonic", "pts", "kernel"]

# -----------------------------------------------------------------

class CutoutMask(np.ndarray):

    """
    This class ...
    """

    def __new__(cls, data, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param cls:
        :param data:
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        if isinstance(data, MaskBase): data = data.data

        obj = np.asarray(data, dtype=bool).view(cls)
        obj.x_min = x_min
        obj.x_max = x_max
        obj.y_min = y_min
        obj.y_max = y_max

        return obj

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.shape[0]

    # -----------------------------------------------------------------

    def as_cutout(self, cutout, padding_value=0.0):

        """
        This function ...
        :param cutout:
        :param padding_value:
        :return:
        """

        rel_x_min = cutout.x_min - self.x_min
        rel_x_max = cutout.x_max - self.x_min

        rel_y_min = cutout.y_min - self.y_min
        rel_y_max = cutout.y_max - self.y_min

        if rel_x_min < 0 or rel_y_min < 0 or rel_x_max > self.xsize or rel_y_max > self.ysize:

            # Create box
            new = np.zeros((cutout.ysize, cutout.xsize))
            if padding_value != 0: new[:,:] = padding_value

            # Normal:
            a = 0
            b = cutout.ysize
            c = 0
            d = cutout.xsize
            aa = rel_y_min
            bb = rel_y_max
            cc = rel_x_min
            dd = rel_x_max

            # Cross-over on the lower x border
            if rel_x_min < 0:
                c = 0 - rel_x_min
                cc = 0

            # Cross-over on the lower y border
            if rel_y_min < 0:
                a = 0 - rel_y_min
                aa = 0

            # Cross-over on the upper x border
            if rel_x_max > self.xsize:

                d = self.xsize - rel_x_min
                dd = self.xsize

            # Cross-over on the upper y border
            if rel_y_max > self.ysize:

                b = self.ysize - rel_y_min
                bb = self.ysize

            new[a:b, c:d] = self[aa::bb, cc::dd]

            # Return the new CutoutMask
            return CutoutMask(new, cutout.x_min, cutout.x_max, cutout.y_min, cutout.y_max)

        # Return a new CutoutMask
        else: return CutoutMask(self[rel_y_min:rel_y_max, rel_x_min:rel_x_max], cutout.x_min, cutout.x_max, cutout.y_min, cutout.y_max)

# -----------------------------------------------------------------

class Cutout(np.ndarray):

    """
    This class ...
    """

    def __new__(cls, data, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param cls:
        :param data:
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # Create an object of the Box class
        obj = np.asarray(data).view(cls)

        # Set attributes of the object
        obj.x_min = x_min
        obj.x_max = x_max
        obj.y_min = y_min
        obj.y_max = y_max

        # Return the object
        return obj

    # -----------------------------------------------------------------

    @classmethod
    def cutout(cls, frame, center, x_radius, y_radius=None):

        """
        This class method ...
        :param frame:
        :param center:
        :param x_radius:
        :param y_radius:
        :return:
        """

        if y_radius is None: y_radius = x_radius

        # Crop the frame
        cropped, x_min, x_max, y_min, y_max = cropping.crop(frame, center.x, center.y, x_radius, y_radius)

        # Check that the center position lies within the box
        assert (x_min <= center.x < x_max and y_min <= center.y < y_max)

        # Return a new box
        return cls(cropped, x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    @classmethod
    def cutout_limits(cls, frame, x_min, x_max, y_min, y_max, absolute=False):

        """
        This function ...
        :param frame:
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :param absolute:
        :return:
        """

        # Crop the frame
        if absolute: cropped = cropping.crop_absolute(frame, x_min, x_max, y_min, y_max)
        else: cropped, x_min, x_max, y_min, y_max = cropping.crop_direct(frame, x_min, x_max, y_min, y_max)

        # Return a new box
        return cls(cropped, x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    @classmethod
    def from_rectangle(cls, frame, rectangle, absolute=False):

        """
        This function ...
        :param frame:
        :param rectangle:
        :param absolute:
        :return:
        """

        # TODO: fix this for rotated rectangles

        # Convert into integers
        x_min = int(round(rectangle.x_min))
        x_max = int(round(rectangle.x_max))
        y_min = int(round(rectangle.y_min))
        y_max = int(round(rectangle.y_max))

        if x_min == x_max:
            x_min = int(math.floor(rectangle.x_min))
            x_max = int(math.ceil(rectangle.x_max))
        if y_min == y_max:
            y_min = int(math.floor(rectangle.y_min))
            y_max = int(math.ceil(rectangle.y_max))

        if x_max - x_min < 2: x_max = x_min + 2
        if y_max - y_min < 2: y_max = y_min + 2

        # Create cutout
        return cls.cutout_limits(frame, x_min, x_max, y_min, y_max, absolute)

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, frame, ellipse, shape=None):

        """
        This function ...
        :param frame:
        :param ellipse:
        :param shape:
        :return:
        """

        # Get bouding box
        rectangle = ellipse.bounding_box

        # ...
        if shape is not None: rectangle = PixelRectangleRegion(rectangle.center, Extent(0.5 * shape.x, 0.5 * shape.y))

        return cls.from_rectangle(frame, rectangle, absolute=(shape is not None))

    # -----------------------------------------------------------------

    @property
    def origin(self):

        """
        This function ...
        :return:
        """

        return PixelCoordinate(self.x_min, self.y_min)

    # -----------------------------------------------------------------

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        rel_y_min = box.y_min - self.y_min
        rel_y_max = box.y_max - self.y_min

        rel_x_min = box.x_min - self.x_min
        rel_x_max = box.x_max - self.x_min

        data = self[rel_y_min:rel_y_max, rel_x_min:rel_x_max]

        # Create the new box and return it
        return Cutout(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # -----------------------------------------------------------------

    def zoom(self, center, factor, min_xpixels=None, min_ypixels=None, min_pixels=None):

        """
        This function ...
        :param center:
        :param factor:
        :param min_xpixels:
        :param min_ypixels:
        :param min_pixels:
        :return:
        """

        # Calculate the size of the smaller box
        new_xsize = int(round(self.xsize / factor))
        new_ysize = int(round(self.ysize / factor))

        # Check min pixels conditions
        if min_xpixels is not None:
            if min_xpixels > self.xsize: raise ValueError("Cannot have a minimum number of x pixels that is larger than the original number of x pixels of the cutout")
            new_xsize = max(new_xsize, min_xpixels)
        if min_ypixels is not None:
            if min_ypixels > self.ysize: raise ValueError("Cannot have a minimum number of y pixels that is larger than the original number of y pixels of the cutout")
            new_ysize = max(new_ysize, min_ypixels)
        if min_pixels is not None:
            original_npixels = self.xsize * self.ysize
            if min_pixels > original_npixels: raise ValueError("Cannot have a minimum number of pixels that is larger than the original number of pixels of the cutout")
            factor_2d = min_pixels / original_npixels
            factor_1d = np.sqrt(factor_2d)
            new_xsize = int(round(self.xsize / factor_1d))
            new_ysize = int(round(self.ysize / factor_1d))
            if new_xsize > self.xsize: raise ValueError("Cannot zoom with a minimum pixels of '" + str(min_pixels))
            if new_ysize > self.ysize: raise ValueError("Cannot zoom with a minimum pixels of '" + str(min_pixels))

        # Calculate the relative coordinate of the center
        rel_center = self.rel_position(center)

        # Calculate radius
        new_xradius = 0.5 * new_xsize
        new_yradius = 0.5 * new_ysize

        # Create a smaller box
        data, rel_x_min, rel_x_max, rel_y_min, rel_y_max = cropping.crop(self, rel_center.x, rel_center.y, new_xradius, new_yradius)

        # Create the new box
        return Cutout(data, rel_x_min+self.x_min, rel_x_max+self.x_min, rel_y_min+self.y_min, rel_y_max+self.y_min)

    # -----------------------------------------------------------------

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        # see InfoArray.__array_finalize__ for comments
        if obj is None: return

        self.x_min = getattr(obj, 'x_min', None)
        self.x_max = getattr(obj, 'x_max', None)
        self.y_min = getattr(obj, 'y_min', None)
        self.y_max = getattr(obj, 'y_max', None)

    # -----------------------------------------------------------------

    def plot(self, frame=None):

        """
        This function ...
        :param frame:
        :return:
        """

        plotting.plot_box(self)

        # If frame is not None, plot 'zoom' plot

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This property ...
        :return:
        """

        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This property ...
        :return:
        """

        return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def x_slice(self):

        """
        This property ...
        :return:
        """

        return slice(self.x_min, self.x_max)

    # -----------------------------------------------------------------

    @property
    def y_slice(self):

        """
        This property ...
        :return:
        """

        return slice(self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def rel_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the relative position
        return position.__class__(position.x - self.x_min, position.y - self.y_min)

    # -----------------------------------------------------------------

    def abs_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the absolute position
        return position.__class__(position.x + self.x_min, position.y + self.y_min)

    # -----------------------------------------------------------------

    def fit_polynomial(self, order, mask=None):

        """
        This function ...
        :param order:
        :param mask:
        :return:
        """

        from .mask import Mask
        if isinstance(mask, Mask): mask = mask.data

        #print(self)
        #print(mask)

        # Do the fitting
        polynomial = fitting.fit_polynomial(self, order, mask=mask)

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

        # Return a new box
        return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def gaussian_filter(self, sigma):

        """
        This function ...
        :param sigma:
        :return:
        """

        # Apply the filter
        data = ndimage.filters.gaussian_filter(self, sigma=sigma)

        # Return a new box
        return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def fit_model(self, center, model_name, sigma=None, max_center_offset=None, amplitude=None, max_sigma_offset=None):

        """
        This function ...
        :param center:
        :param model_name:
        :param sigma:
        :param max_center_offset:
        :param amplitude:
        :param max_sigma_offset:
        :return:
        """

        # Calculate the relative coordinate of the center
        rel_center = self.rel_position(center)

        # Fit a 2D Gaussian to the data
        if model_name == "Gaussian":

            # Do the fitting
            model = fitting.fit_2D_Gaussian(self, rel_center, sigma=sigma, max_center_offset=max_center_offset, amplitude=amplitude, max_sigma_offset=max_sigma_offset)

        # Fit an Airy Disk model to the data
        elif model_name == "Airy":

            # See https://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile and
            # http://astropy.readthedocs.org/en/latest/api/astropy.modeling.functional_models.AiryDisk2D.html#astropy.modeling.functional_models.AiryDisk2D
            #sigma = 0.42 * airy.radius * 0.81989397882
            radius = sigma / (0.42 * 0.81989397882) if sigma is not None else None

            # Do the fitting
            model = fitting.fit_2D_Airy(self, rel_center, radius=radius, max_center_offset=max_center_offset, amplitude=amplitude)

        # Fit a 2D (vertically) shifted Gaussian to the data
        elif model_name == "ShiftedGaussian":

            # Do the fitting
            model = fitting.fit_2D_ShiftedGaussian(self, rel_center, sigma=sigma, max_center_offset=max_center_offset, amplitude=amplitude)

        # Unknown model name
        else: raise ValueError("Model name should be 'Gaussian' or 'Airy'")

        # Shift the position of the model so that it represents its absolute position in the frame
        fitting.shift_model(model, self.x_min, self.y_min)

        # Return the model
        return model

    # -----------------------------------------------------------------

    def evaluate_model(self, model):

        """
        This function ...
        :param model:
        :return:
        """

        # Make a local copy of the model so that we can adapt its position to be relative to this box
        rel_model = fitting.shifted_model(model, -self.x_min, -self.y_min)

        # Create x and y meshgrid for evaluating
        y_values, x_values = np.mgrid[:self.ysize, :self.xsize]

        # Evaluate the model
        data = rel_model(x_values, y_values)

        # Return a new box
        return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def interpolated(self, mask, method, no_clip_mask=None, plot=False, sigma=None, fwhm=None, max_iterations=10,
                     not_converge="keep", dilation_radius_steps=2):

        """
        This function ...
        :param mask:
        :param method:
        :param no_clip_mask:
        :param plot:
        :param sigma:
        :param fwhm:
        :param max_iterations:
        :param not_converge:
        :param dilation_radius_steps:
        :return:
        """

        # Fit a polynomial to the data
        if method == "polynomial":
            try:
                return self.fit_polynomial(3, mask=mask)
            except TypeError:
                #from ..tools import plotting
                log.debug("Error while fitting polynomial to the box ...")
                log.debug("cutout = " + str(type(self)) + " " + str(mask.shape))
                log.debug("mask = " + str(type(mask)) + " " + str(mask.shape))
                #print("HERE!")
                #plotting.plot_box(self)
                #plotting.plot_box(mask)
                #plotting.plot_box(np.ma.masked_array(self, mask=mask))
                #print("AND HEREEE!!!")
                mask = mask.eroded_rc(connectivity=2, iterations=1)
                #plotting.plot_box(np.ma.masked_array(self, mask=mask))
                try: return self.fit_polynomial(3, mask=mask)
                except TypeError: return self.fit_polynomial(3)

        # Interpolate using the local mean method
        elif method == "local_mean":

            try:
                # Calculate the interpolated data
                data = interpolation.in_paint(self, mask)
            except IndexError:
                log.debug("Error while inpainting using the local_mean method ...")
                data = np.zeros((self.ysize, self.xsize))

            # Create and return a new box
            return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Interpolate using inverse distance weighing
        elif method == "idw":

            # Calculate the interpolated data
            data = interpolation.in_paint(self, mask, method="idw")

            # Create and return a new box
            return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Calculate the mean value of the data
        elif method == "mean":

            mean = np.ma.mean(np.ma.masked_array(self, mask=mask))
            data = np.full(self.shape, mean)
            return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Calculate the median value of the data
        elif method == "median":

            median = float(np.ma.median(np.ma.masked_array(self, mask=mask)))
            data = np.full(self.shape, median)
            return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Use the biharmonic method of the Scikit-image package
        elif method == "biharmonic":

            data = interpolation.inpaint_biharmonic(self, mask)

            # Check whether the data is NaN-free, because this method will sometimes fail to interpolate and
            # just put NaNs for the pixels that are covered by the mask. (I still have to check under which
            # conditions but it has something to do with being close to the boundary)
            if np.any(np.isnan(data)):

                #from ..tools import plotting
                #plotting.plot_box(normalized)
                #plotting.plot_box(mask)

                # Interpolate by local_mean, this does not leave nans
                data = interpolation.in_paint(data, mask)

            return Cutout(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Use the kernel  method
        elif method == "kernel":

            from ..tools import statistics
            from astropy.convolution import Gaussian2DKernel

            # Get the sigma
            if fwhm is not None:
                sigma_from_fwhm = fwhm * statistics.fwhm_to_sigma
                if sigma is not None and not np.isclose(sigma, sigma_from_fwhm): raise ValueError("Sigma and FWHM do not correspond")
                sigma = sigma_from_fwhm
            elif sigma is None: raise ValueError("Either FWHM or sigma has to be specified")

            # Debugging
            log.debug("The sigma for kernel interpolation is " + str(sigma) + " pixels")

            # We smooth with a Gaussian kernel with stddev passed by the user
            # Create the kernel
            kernel = Gaussian2DKernel(stddev=sigma)

            # Create data with NaNs
            data = self.copy()

            # Convert to new mask object
            if isinstance(mask, oldMask): mask = Mask(mask)

            # Interpolate
            result = interpolate_kernel(data, mask, kernel, max_iterations=max_iterations, not_converge=not_converge, plot=plot)

            # Interpolate around the mask
            previous_mask = mask
            while True:

                new_data = data.copy()
                new_mask = previous_mask.disk_dilated(radius=dilation_radius_steps)
                if new_mask.all_masked: break

                # Interpolate
                new_result = interpolate_kernel(new_data, new_mask, kernel, max_iterations=max_iterations, not_converge="error", ignore_nans_in=previous_mask, plot=plot)

                # Where new mask and not previous mask
                where = new_mask * np.logical_not(previous_mask)

                # Set data
                result[where] = new_result[where]

                # Set previous mask
                previous_mask = new_mask

            # Return the result
            return Cutout(result, self.x_min, self.x_max, self.y_min, self.y_max)

        # Use the 'PTS' method
        elif method == "pts":

            #test = no_clip_mask is not None

            # POLYNOMIAL WITHOUT SIGMA-CLIPPING

            """
            polynomial_fit_mask = mask if no_clip_mask is None else no_clip_mask  # determine the mask to use for fitting the polynomial
            order = 3
            try:
                polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
            except TypeError:
                try:
                    polynomial_fit_mask = polynomial_fit_mask.eroded(connectivity=2, iterations=1)
                    polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
                except TypeError:
                    log.warning("Cannot interpolate over this box ...")
                    return self.copy()

            # Evaluate the polynomial
            poly_data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

            subtracted = self - poly_data

            if test: plotting.plot_difference(self, poly_data, title="without sigma-clipping")
            """

            #no_clip_mask = None
            #polynomial_fit_mask = mask if no_clip_mask is None else no_clip_mask # determine the mask to use for fitting the polynomial

            polynomial_fit_mask = mask

            order = 3
            try: polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
            except TypeError:
                try:
                    polynomial_fit_mask = polynomial_fit_mask.eroded(connectivity=2, iterations=1)
                    polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
                except TypeError:
                    log.warning("Cannot interpolate over this box ...")
                    return self.copy()

            # Evaluate the polynomial
            poly_data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

            subtracted = self - poly_data

            if plot: plotting.plot_difference(self, poly_data, title="with sigma-clipping")

            #if no_clip_mask is not None:
            #    plotting.plot_difference(self, poly_data)
            #    plotting.plot_box(subtracted)
            #    plotting.plot_box(np.ma.array(subtracted, mask=mask))
            #    from ..tools import statistics
            #    new_mask = statistics.sigma_clip_mask(subtracted, sigma_level=2.0, mask=mask)
            #    plotting.plot_box(np.ma.array(subtracted, mask=new_mask))
            #    background_pixels = subtracted[new_mask.inverse()]
            #else:

            #if test:
            #    new_mask = statistics.sigma_clip_mask(subtracted, sigma_level=2.0, mask=mask)
            #    background_pixels = subtracted[new_mask.inverse()]
            #else: background_pixels = subtracted[mask.inverse()]

            background_pixels = subtracted[mask.inverse()]

            #distribution = Distribution.from_values(background_pixels)
            #gaussian = distribution.fit_gaussian()
            #center = gaussian.mean
            #stddev = gaussian.stddev
            #distribution.plot(model=gaussian)

            center = np.mean(background_pixels)
            stddev = np.std(background_pixels)

            #print("amplitude", gaussian.amplitude)
            #print("center", center)
            #print("stddev", stddev)

            # Generate random pixel values
            random = np.random.normal(center, stddev, self.shape)

            if plot: plotting.plot_box(random, title="random (center value=" + str(center) + " stddev=" + str(stddev) + ")")

            # Fit polynomial again but with no-sigma-clipped-mask

            if no_clip_mask is not None:

                polynomial_fit_mask = no_clip_mask
                order = 3
                try:
                    polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
                except TypeError:
                    try:
                        polynomial_fit_mask = polynomial_fit_mask.eroded(connectivity=2, iterations=1)
                        polynomial = fitting.fit_polynomial(self, order, mask=polynomial_fit_mask)
                    except TypeError:
                        log.warning("Cannot interpolate over this box ...")
                        return self.copy()

                # Evaluate the polynomial
                poly_data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

                if plot: plotting.plot_difference(self, poly_data, title="polynomial with sigma-clipping")

            return Cutout(random + poly_data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Invalid option
        else: raise ValueError("Unknown interpolation method")

    # -----------------------------------------------------------------

    def value(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Get the x and y coordinate of the corresponding pixel
        x = int(round(position.x - self.x_min))
        y = int(round(position.y - self.y_min))

        # Return the pixel value
        return self[y, x]

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Convert to relative position
        rel_position = self.rel_position(position)

        # Calculate the pixel coordinates of the position
        x_pixel = int(round(rel_position.x))
        y_pixel = int(round(rel_position.y))

        # Check whether this box contains the position
        if x_pixel < 0 or y_pixel < 0 or x_pixel >= self.xsize or y_pixel >= self.ysize: return False
        else: return True

    # -----------------------------------------------------------------

    def replace(self, frame, where=None):

        """
        This function ...
        :param frame:
        :param where:
        :return:
        """

        # Replace the pixel values in the frame
        if where is None: frame[self.y_min:self.y_max, self.x_min:self.x_max] = self
        elif isinstance(where, AlphaMask):

            fractional = where.as_real()
            frame[self.y_min:self.y_max, self.x_min:self.x_max] = self * fractional + frame[self.y_min:self.y_max, self.x_min:self.x_max] * (1. - fractional)

        else: frame[self.y_min:self.y_max, self.x_min:self.x_max][where] = self[where]

        #from ..tools import plotting
        #plotting.plot_box(frame[self.y_min:self.y_max, self.x_min:self.x_max])
        #plotting.plot_box(np.ma.masked_array(frame[self.y_min:self.y_max, self.x_min:self.x_max], mask=where))

# -----------------------------------------------------------------

def interpolate_kernel(data, mask, kernel, max_iterations=10, not_converge="keep", ignore_nans_in=None, plot=False):

    """
    Thisf unction ...
    :param data:
    :param mask:
    :param kernel:
    :param max_iterations:
    :param not_converge:
    :param ignore_nans_in:
    :param plot:
    :return:
    """

    # Set Nans
    data[mask] = float("NaN")

    # Interpolate over nans
    return interpolate_nans_kernel(data, kernel, max_iterations=max_iterations, not_converge=not_converge, ignore_nans_in=ignore_nans_in, plot=plot)

# -----------------------------------------------------------------

def interpolate_nans_kernel(data, kernel, max_iterations=10, not_converge="keep", ignore_nans_in=None, plot=False):

    """
    This function ...
    :param data:
    :param kernel:
    :param max_iterations:
    :param not_converge:
    :param ignore_nans_in:
    :param plot:
    :return:
    """

    from astropy.convolution import interpolate_replace_nans

    # Debugging
    nans = np.isnan(data)
    nnans = np.sum(nans)
    log.debug("The number of NaNs at the start is " + str(nnans))

    # Plot?
    if plot: plotting.plot_mask(nans, title="NaNs")

    # Interpolate
    result = interpolate_replace_nans(data, kernel)
    niterations = 1

    # Plot?
    if plot: plotting.plot_mask(result, title="Result after iteration 1")

    # Plot?
    if plot and ignore_nans_in is not None: plotting.plot_mask(ignore_nans_in, title="ignoring NaNs")

    # Get the current number of nans
    # previous_nnans = None  # don't get it for performance
    if ignore_nans_in is not None: nnans = np.sum(np.isnan(result[np.logical_not(ignore_nans_in)]))
    else: nnans = np.sum(np.isnan(result))

    # Plot?
    if plot: plotting.plot_mask(np.isnan(result), title="NaNs after iteration 1")

    # Debugging
    log.debug("The number of NaN values after iteration 1 is " + str(nnans))

    # Are there still NaNs?
    while nnans > 0:

        # Check number of iterations
        if max_iterations is not None and niterations == max_iterations: raise RuntimeError("The maximum number of iterations has been reached without success")

        # Debugging
        log.debug("Interpolation iteration " + str(niterations + 1) + " ...")

        # Perform next interpolation
        result = interpolate_replace_nans(result, kernel)

        # Increment the niterations counter
        niterations += 1

        # Get the current number of nans
        previous_nnans = nnans
        if ignore_nans_in is not None: nnans = np.sum(np.isnan(result[np.logical_not(ignore_nans_in)]))
        else: nnans = np.sum(np.isnan(result))

        # Plot?
        if plot: plotting.plot_mask(np.isnan(result), title="NaNs after iteration " + str(niterations))

        condition = nnans > previous_nnans if ignore_nans_in else nnans >= previous_nnans

        # Not converging
        if condition:
            if not_converge == "keep":
                log.warning("The number of NaNs could not converge to zero: " + str(nnans) + " NaN values will remain")
                break  # break the loop
            elif not_converge == "error": raise RuntimeError("The number of NaNs is not converging to zero (nnans = " + str(nnans) + ", previous nnans = " + str(previous_nnans) + ")")
            else: raise ValueError("Invalid option for 'not_converge'")

        # Debugging
        log.debug("The number of NaN values after iteration " + str(niterations) + " is " + str(nnans))

    # Return the resulting data
    return result

# -----------------------------------------------------------------
