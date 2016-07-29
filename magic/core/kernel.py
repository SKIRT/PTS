#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.kernel Contains the ConvolutionKernel class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from scipy import ndimage
from scipy.ndimage.interpolation import shift, zoom

# Import astronomical modules
from astropy.modeling import models, fitting
from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg

# Import the relevant PTS classes and modules
from .frame import Frame
from ...core.tools.logging import log
from ..tools import statistics

# -----------------------------------------------------------------

class ConvolutionKernel(Frame):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, data, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(ConvolutionKernel, self).__init__(data, *args, **kwargs)

        # Check the FWHM
        if self.fwhm is None: raise ValueError("FWHM must be specified if not present in header")

        # Set the WCS to None, but keep the pixelscale
        if self._wcs is not None:

            if self._pixelscale is None: self._pixelscale = self.wcs.pixelscale
            elif not np.isclose(self._pixelscale, self.wcs.pixelscale): raise ValueError("Pixelscale in the header does not correspond to the specified pixelscale")
            self._wcs = None

        elif self._pixelscale is None: raise ValueError("Pixelscale must be specified if not present in header")

        # Prepared
        self._prepared = False

        if "prepared" in self.meta: self._prepared = self.meta["prepared"]

        # Make sure that the data is in 64 bit floating-point precision
        self._data = self._data.astype("float64")

    # -----------------------------------------------------------------

    @property
    def prepared(self):

        """
        This function ...
        :return:
        """

        return self._prepared

    # -----------------------------------------------------------------

    @property
    def normalized(self):

        """
        This function ...
        :return:
        """

        #return np.abs(self.sum() - 1.) < 1e-7  # criterion as in astropy.convolution module = 1e-8, BUT I HAD A PROBLEM OF THE SAME FITS FILE HAVING A SLIGHTLY DIFFERNENT SOME ON DIFFERENT SYSTEMS !!! (laptop and nancy)

        # NOW I KNOW WHERE THIS PROBLEM WAS COMING FROM (SEE ASTROPY ISSUE #5176)
        # the magical number of 32-bit precision floats is 1.1920928955078125e-07

        # BUT NOW I MAKE SURE THAT A KERNEL IS ALWAYS IN 64-BIT REPRESENTATION AND SO 1e-8 CAN BE USED

        return np.abs(self.sum() - 1.) < 1e-8

    # -----------------------------------------------------------------

    def prepare_for(self, image, sigma_level=10.0):

        """
        This function ...
        :param image:
        :param sigma_level:
        :return:
        """

        # Prepare
        self.prepare(image.pixelscale, sigma_level)

    # -----------------------------------------------------------------

    def prepare_chris(self, pixelscale):

        """
        This function ...
        :return:
        """

        av_pixelscale_arcsec = pixelscale.average.to("arcsec/pix").value
        kernel_av_pixelscale_arcsec = self.average_pixelscale.to("arcsec/pix").value

        # If PSF pixel size is different to map pixel size, rescale PSF accordingly
        if (av_pixelscale_arcsec / kernel_av_pixelscale_arcsec) > 1.001 or (av_pixelscale_arcsec / kernel_av_pixelscale_arcsec) < 0.999:

            psf_in = self._data

            # 'Trim' the edges of the input PSF until the rescaled PSF has odd dimensions
            psf_even = True
            while psf_even:

                zoom_factor = float(kernel_av_pixelscale_arcsec) / float(av_pixelscale_arcsec)
                psf = zoom(psf_in, (zoom_factor, zoom_factor), mode='nearest')
                if (psf.shape[0] % 2 != 0) and (psf.shape[1] % 2 != 0):
                    psf_even = False
                else:
                    psf_in = psf_in[1:, 1:]
                    psf_in = psf_in[:-1, :-1]

            self._data = psf_in

        # Else, if pixel sizes are already the same, leave as-is

        # Normalise PSF
        self._data /= np.nansum(self._data)

    # -----------------------------------------------------------------

    def prepare(self, pixelscale, sigma_level=10.0):

        """
        This function ...
        :param pixelscale:
        :param sigma_level:
        :return:
        """

        # Inform the user
        log.info("Preparing the kernel ...")

        # Truncate
        self.truncate(sigma_level)

        # Adjust pixelscale
        self.adjust_pixelscale(pixelscale)

        # Recenter
        self.recenter()

        # Normalize
        self.normalize()

        # Set prepared flag to True
        self._prepared = True

    # -----------------------------------------------------------------

    def truncate(self, sigma_level=10.0):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Truncating the kernel to a sigma level of " + str(sigma_level) + " ...")

        # Determine the radius in number of pixels
        sigma_pix = statistics.fwhm_to_sigma * self.fwhm_pix
        radius = sigma_level * sigma_pix

        center_x = 0.5 * (self.xsize - 1.)
        center_y = 0.5 * (self.ysize - 1.)

        min_x = int(round(center_x - radius))
        max_x = int(round(center_x + radius))

        min_y = int(round(center_y - radius))
        max_y = int(round(center_y + radius))

        # Crop
        self.crop(min_x, max_x, min_y, max_y)

    # -----------------------------------------------------------------

    def adjust_pixelscale(self, pixelscale):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adjusting the pixelscale of the kernel to match that of the image ...")

        average_pixelscale = 0.5 * (pixelscale.x + pixelscale.y)

        # Calculate the zooming factor
        factor = (average_pixelscale / self.average_pixelscale).to("").value

        # Rebin to the pixelscale
        new_data = ndimage.interpolation.zoom(self._data, zoom=1.0 / factor)

        # Set the new data and pixelscale
        self._data = new_data
        self._pixelscale = pixelscale

    # -----------------------------------------------------------------

    def recenter(self, centroid_method="2dg"):

        """
        This function ...
        :return:
        """

        center_x = 0.5 * (self.xsize - 1)
        center_y = 0.5 * (self.ysize - 1)

        if centroid_method == "com": x_centroid, y_centroid = self.centroid_com()
        elif centroid_method == "fit": x_centroid, y_centroid = self.centroid_fit()
        elif centroid_method == "2dg": x_centroid, y_centroid = self.centroid_2dg()
        elif centroid_method == "aniano": x_centroid, y_centroid = self.get_maximum_aniano()
        else: raise ValueError("Invalid centroid method")

        # Debugging
        log.debug("The centroid coordinate of the kernel was found to be " + str(x_centroid) + ", " + str(y_centroid))
        log.debug("The center of the kernel image is " + str(center_x) + ", " + str(center_y))

        # Calculate shift
        shift_x = center_x - x_centroid
        shift_y = center_y - y_centroid

        # Debugging
        log.debug("The required shift is (" + str(shift_x) + ", " + str(shift_y) + ")")

        # If the shift (in ABSOLUTE VALUE) is less than 0.2 pixel, don't shift
        if abs(shift_x) < 0.2 and abs(shift_y) < 0.2:

            log.debug("Kernel is already perfectly aligned with the center: skipping recentering ...")
            return

        # Debugging
        log.debug("Shifting the kernel center by (" + str(shift_x) + ", " + str(shift_y) + ") pixels ...")

        # Shift
        self._data = shift(self._data, [shift_x, shift_y])

        # CHECK AGAIN

        if centroid_method == "com": x_centroid, y_centroid = self.centroid_com()
        elif centroid_method == "fit": x_centroid, y_centroid = self.centroid_fit()
        elif centroid_method == "2dg": x_centroid, y_centroid = self.centroid_2dg()
        elif centroid_method == "aniano": x_centroid, y_centroid = self.get_maximum_aniano()
        else: raise ValueError("Invalid centroid method")

        new_shift_x = center_x - x_centroid
        new_shift_y = center_y - y_centroid

        new_shift_x_relative = abs(new_shift_x) / abs(shift_x)
        new_shift_y_relative = abs(new_shift_y) / abs(shift_y)

        #print("new shift x relative " + str(new_shift_x_relative))
        #print("new shift y relative " + str(new_shift_y_relative))

        if new_shift_x_relative >= 0.1: raise RuntimeError("The recentering of the kernel failed: new x shift = " + str(new_shift_x) + ", previous x shift = " + str(shift_x))
        if new_shift_y_relative >= 0.1: raise RuntimeError("The recentering of the kernel failed: new y shift = " + str(new_shift_y) + ", previous y shift = " + str(shift_y))

    # -----------------------------------------------------------------

    def centroid_com(self):

        """
        This function ...
        :return:
        """

        return centroid_com(self._data)

    # -----------------------------------------------------------------

    def centroid_2dg(self):

        """
        This function ...
        :return:
        """

        return centroid_2dg(self._data)

    # -----------------------------------------------------------------

    def center_fit(self):

        """
        This function ...
        :return:
        """

        from .box import Box
        from ..basics.vector import Position

        # Box
        box = Box(self._data, 0, self.xsize, 0, self.ysize)

        # Fit model
        model = box.fit_model(Position(0.5*(self.xsize-1), 0.5*(self.ysize-1)), "Gaussian")

        # Get x and y mean
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        return x_mean, y_mean

    # -----------------------------------------------------------------

    def center_aniano(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Centering the kernel ...")

        # FROM CONVOLVE_IMAGE.PRO (G. Aniano)

        center_x = int(0.5 * (self.xsize - 1))
        center_y = int(0.5 * (self.ysize - 1))

        # get_maximun,image,x_max,y_max

        x_max, y_max = self.get_maximum()

        # ; determine the needed shifts
        shift_x = center_x - x_max
        shift_y = center_y - y_max

        # ; make the shift if nonzero
        if (shift_x != 0) or (shift_y != 0):

            # Debugging
            log.debug("Shifting the kernel center by (" + str(shift_x) + ", " + str(shift_y) + ") pixels ...")

            self._data = shift(self._data, [shift_x,shift_y])

            # Y
            self._data[:abs(shift_y),:] = 0.0
            self._data[self.ysize-1-abs(shift_y):self.ysize,:] = 0.0

            # X
            self._data[:,:abs(shift_x)] = 0.0
            self._data[:,self.xsize-1-abs(shift_x):] = 0.0

        # CHECK

        # Calculate shift again
        x_max, y_max = self.get_maximum()
        new_shift_x = center_x - x_max
        new_shift_y = center_y - y_max

        # Raise exception if there is still a shift
        if (new_shift_x != 0) or (new_shift_y != 0): raise RuntimeError("Something went wrong during the kernel centering: "
                                                                "new shift x = " + str(new_shift_x) + ", new shift y = "
                                                                + str(new_shift_y) + " (previous shift x = " + str(shift_x)
                                                                        + ", previous shift y = " + str(shift_y))

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        self.__idiv__(self.sum())

    # -----------------------------------------------------------------

    def get_maximum_aniano(self):

        """
        This function ...
        :return:
        """

        rad_to_mean = 5

        data_copy = self._data.copy()

        #
        mean_im = data_copy * 0.0

        i_range = range(-int(rad_to_mean), int(rad_to_mean)+1)
        #print("irange", i_range)

        for i in i_range:
           j_range = range(-int(math.sqrt(rad_to_mean ** 2 - i ** 2)), int(math.sqrt(rad_to_mean ** 2 - i ** 2))+1)
           #print("jrange", j_range)
           for j in j_range:
              mean_im += shift(data_copy, [i, j])

        mean_im_sum = np.sum(mean_im)

        #mx    = max(mean_im, location)
        #index = ARRAY_INDICES(mean_im, location)
        #x_max = index[0]
        #y_max = index[1]

        # Get x and y max
        max_index = np.argmax(mean_im)
        c = (max_index // len(mean_im[0]), max_index % len(mean_im[0]))
        x_max = c[1]
        y_max = c[0]

        max_value = mean_im[y_max, x_max]
        where = np.abs(mean_im - max_value) < (5e-4 * mean_im_sum)
        count = np.sum(where)

        if count > 1:

            log.debug("WARNING: The PSF has " + str(count) + "pixels with values similar to its maximum... we will take their centroid...")

            xsize = data_copy.shape[1]
            ysize = data_copy.shape[0]

            xv, yv = np.meshgrid(np.arange(xsize), np.arange(ysize))

            # Average x max
            x_max = np.sum(xv[where]) / float(count)

            # Average y max
            y_max = np.sum(xv[where]) / float(count)

        # Return xmax and ymax position
        return x_max, y_max

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the save function of the base class
        super(ConvolutionKernel, self).save(path, extra_header_info={"PREPARED": self.prepared})

# -----------------------------------------------------------------
