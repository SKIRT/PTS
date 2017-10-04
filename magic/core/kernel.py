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
#from astropy.modeling import models, fitting
from photutils.centroids import centroid_com, centroid_1dg, centroid_2dg
from astropy.modeling.models import Gaussian2D, AiryDisk2D
from astropy.convolution.kernels import Gaussian2DKernel, AiryDisk2DKernel

# Import the relevant PTS classes and modules
from .frame import Frame
from ...core.basics.log import log
from ..tools import statistics
from ...core.filter.filter import parse_filter
from ..tools import fitting

# -----------------------------------------------------------------

class ConvolutionKernel(Frame):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, data, *args, **kwargs):

        """
        This function ...
        :param data:
        :param args:
        :param kwargs:
        """

        # Set kwargs into meta so that they can be read more below in this constructor
        if "meta" not in kwargs: kwargs["meta"] = dict()
        if "from_filter" in kwargs: kwargs["meta"]["frmfltr"] = str(kwargs.pop("from_filter"))
        if "to_filter" in kwargs: kwargs["meta"]["tofltr"] = str(kwargs.pop("to_filter"))
        if "prepared" in kwargs: kwargs["meta"]["prepared"] = str(kwargs.pop("prepared"))

        # Call the constructor of the base class
        super(ConvolutionKernel, self).__init__(data, *args, **kwargs)

        # Check the FWHM
        if self.fwhm is None: raise ValueError("FWHM must be specified if not present in header")

        # Set the WCS to None, but keep the pixelscale
        if self._wcs is not None:

            if self._pixelscale is None: self._pixelscale = self.wcs.pixelscale
            elif not np.isclose(self._pixelscale, self.wcs.pixelscale): raise ValueError("Pixelscale in the header does not correspond to the specified pixelscale")
            self._wcs = None

        # If the WCS is not defined, the pixelscale must
        elif self._pixelscale is None: raise ValueError("Pixelscale must be specified if not present in header")

        # Prepared
        self._prepared = False
        if "prepared" in self.metadata: self._prepared = self.metadata["prepared"]

        # Get from filter
        if "frmfltr" in self.metadata:
            # Prevent sending None string to parse_filter
            if self.metadata["frmfltr"] != 'None':
                self.from_filter = parse_filter(self.metadata["frmfltr"])
            else:
                self.from_filter = None
        else: self.from_filter = None

        # Get target filter
        if "tofltr" in self.metadata:
            if self.metadata["tofltr"] is None or self.metadata["tofltr"] == "None":
                log.warning("The target filter for this convolution kernel is not defined")
                self.to_filter = None
                self._psf_filter = None
            else:
                self.to_filter = parse_filter(self.metadata["tofltr"])
                self._psf_filter = self.to_filter
        else:
            log.warning("The target filter for this convolution kernel is not defined")
            self.to_filter = None
            self._psf_filter = None

        # Make sure that the data is in 64 bit floating-point precision (the reason is the precision of Astropy's normalization criterion)
        self._data = self._data.astype("float64")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, fwhm=None, from_filter=None, to_filter=None):

        """
        This function ...
        :param path: 
        :param fwhm: 
        :param from_filter: 
        :param to_filter: 
        :return: 
        """

        # Set extra meta
        extra_meta = dict()
        if from_filter is not None: extra_meta["frmfltr"] = str(from_filter)
        if to_filter is not None: extra_meta["tofltr"] = str(to_filter)

        # Call the from_file function of the base class
        return super(ConvolutionKernel, cls).from_file(path, fwhm=fwhm, add_meta=True, extra_meta=extra_meta, no_filter=True, no_wcs=True)

    # -----------------------------------------------------------------

    @classmethod
    def from_model(cls, model, from_filter=None, to_filter=None):

        """
        This function ...
        :param model: 
        :param from_filter:
        :param to_filter:
        :return: 
        """

        # Get properties
        #center =
        #sigma = fitting.sigma_symmetric(model)

        # Create a kernel
        #kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)
        #kernel.normalize()  # to suppress warning

        if isinstance(model, Gaussian2D):
            x_stddev = model.x_stddev
            y_stddev = model.y_stddev
            if not np.isclose(x_stddev, y_stddev): raise ValueError("Model is not symmetric")
            kernel = Gaussian2DKernel(x_stddev)
        elif isinstance(model, AiryDisk2D):
            radius = model.radius
            kernel = AiryDisk2DKernel(radius)
        else: raise ValueError("Model not supported")
        kernel.normalize()

        # Get the FWHM
        fwhm = fitting.fwhm_symmetric(model)

        # Set metadata
        extra_meta = dict()
        if from_filter is not None: extra_meta["frmfltr"] = str(from_filter)
        if to_filter is not None: extra_meta["tofltr"] = str(to_filter)

        # Create instance of this class
        return cls(kernel.array, fwhm=fwhm, extra_meta=extra_meta, prepared=True)

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

        #return np.abs(self.sum() - 1.) < 1e-7  # criterion as in astropy.convolution module = 1e-8,
        # BUT I HAD A PROBLEM OF THE SAME FITS FILE HAVING A SLIGHTLY DIFFERNENT SOME ON DIFFERENT SYSTEMS !!! (laptop and nancy)

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

    @property
    def odd_xsize(self):

        """
        This function ...
        :return:
        """

        return self.xsize % 2 != 0

    # -----------------------------------------------------------------

    @property
    def odd_ysize(self):

        """
        This function ...
        :return:
        """

        return self.ysize % 2 != 0

    # -----------------------------------------------------------------

    @property
    def x_fwhms(self):

        """
        Thisf nction ...
        :return:
        """

        return 0.5 * float(self.xsize) / self.fwhm_pix

    # -----------------------------------------------------------------

    @property
    def y_fwhms(self):

        """
        This function ...
        :return:
        """

        return 0.5 * float(self.ysize) / self.fwhm_pix

    # -----------------------------------------------------------------

    @property
    def fwhms(self):

        """
        This function ...
        :return:
        """

        xfwhms = self.x_fwhms
        yfwhms = self.y_fwhms

        if np.isclose(xfwhms, yfwhms, rtol=0.05): return np.mean([xfwhms, yfwhms])
        else: raise ValueError("The number of FWHMs along the x and y axis differs more than 5%")

    # -----------------------------------------------------------------

    @property
    def sigma(self):

        """
        This function ...
        :return:
        """

        return self.fwhm * statistics.fwhm_to_sigma

    # -----------------------------------------------------------------

    @property
    def sigma_pix(self):

        """
        This function ...
        :return:
        """

        return (self.sigma / self.average_pixelscale).to("").value if self.sigma is not None else None

    # -----------------------------------------------------------------

    @property
    def x_sigmas(self):

        """
        Thisf unction ...
        :return:
        """

        return 0.5 * float(self.xsize) / self.sigma_pix

    # -----------------------------------------------------------------

    @property
    def y_sigmas(self):

        """
        This function ...
        :return:
        """

        return 0.5 * float(self.ysize) / self.sigma_pix

    # -----------------------------------------------------------------

    @property
    def sigmas(self):

        """
        This function ...
        :return:
        """

        xsigmas = self.x_sigmas
        ysigmas = self.y_sigmas

        if np.isclose(xsigmas, ysigmas, rtol=0.05): return np.mean([xsigmas, ysigmas])
        else: raise ValueError("The number of sigmas along the x and y axis differs more than 5%")

    # -----------------------------------------------------------------

    @property
    def sigma_level(self):

        """
        This function ...
        :return:
        """

        return self.sigmas

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

        # Check the dimensions
        self.check_dimensions()

        # Truncate
        if sigma_level is not None: self.truncate(sigma_level)

        # Adjust pixelscale
        self.adjust_pixelscale(pixelscale)

        # Recenter
        self.recenter()

        # Normalize
        self.normalize()

        # Set prepared flag to True
        self._prepared = True

    # -----------------------------------------------------------------

    def check_dimensions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the dimensions of the kernel image ...")

        if self.odd_xsize and self.odd_ysize:

            log.info("Dimensions are OK")
            return # everything normal

        elif (not self.odd_xsize) and self.odd_ysize: # only odd in y direction

            # Trim in x direction
            self._data = self._data[:, 1:]

        elif not self.odd_xsize and self.odd_ysize: # only odd in x direction

            # Trim in y direction
            self._data = self._data[1:, :]

        # trim in x and y
        else: self._data = self._data[1:, 1:]

        # Recentering
        self.recenter()

    # -----------------------------------------------------------------

    def truncate(self, sigma_level=10.0):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Truncating the kernel to a sigma level of " + str(sigma_level) + " ...")

        # Determine the current sigma level
        current_sigma_level = self.sigmas
        if current_sigma_level < sigma_level:
            log.warning("The current sigma level of the kernel (" + str(current_sigma_level) + ") is smaller than the requested sigma level")
            return

        # Determine the radius in number of pixels
        sigma_pix = statistics.fwhm_to_sigma * self.fwhm_pix
        radius = sigma_level * sigma_pix
        radius_npixels = int(round(radius))

        center_x = 0.5 * (self.xsize - 1.)
        center_y = 0.5 * (self.ysize - 1.)

        assert int(center_x) == center_x
        assert int(center_x) == center_x

        center_x = int(center_x)
        center_y = int(center_y)

        # Determine limits
        min_x = center_x - radius_npixels
        max_x = center_x + radius_npixels + 1

        min_y = center_y - radius_npixels
        max_y = center_y + radius_npixels + 1

        # Crop
        self.crop(min_x, max_x, min_y, max_y, out_of_bounds="adjust")

    # -----------------------------------------------------------------

    def adjust_pixelscale(self, pixelscale):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adjusting the pixelscale of the kernel to match that of the image ...")

        # Calculate the zooming factor
        x_zoom_factor = (self.pixelscale.x / pixelscale.x).to("").value
        y_zoom_factor = (self.pixelscale.y / pixelscale.y).to("").value

        # If PSF pixel size is different to map pixel size, rescale PSF accordingly
        if 0.999 < x_zoom_factor < 1.001 and 0.999 < y_zoom_factor < 1.001:

            # Inform the user and return
            log.info("The pixelscale of the kernel is already (almost) identical to that of the image: no adjustment done")
            return

        # Loop until xsize and ysize are odd (so that kernel center is at exactly the center pixel)
        while True:

            # Rebin to the pixelscale
            new_data = ndimage.interpolation.zoom(self._data, zoom=(y_zoom_factor, x_zoom_factor), mode="nearest") # previously, mode was not specified

            if (new_data.shape[0] % 2 != 0) and (new_data.shape[1] % 2 != 0): break
            else:

                # 'Trim' the edges of the input PSF until the rescaled PSF has odd dimensions
                self._data = self._data[1:, 1:]
                self._data = self._data[:-1, :-1]

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
        #self._data = shift(self._data, [shift_x, shift_y])
        self._data = shift(self._data, [shift_y, shift_x])

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

        #from ..tools import plotting
        #plotting.plot_box(self._data)
        return centroid_2dg(self._data)

    # -----------------------------------------------------------------

    def center_fit(self):

        """
        This function ...
        :return:
        """

        from .cutout import Cutout
        from ..basics.vector import Position

        # Box
        box = Cutout(self._data, 0, self.xsize, 0, self.ysize)

        # Fit model
        model = box.fit_model(Position(0.5*(self.xsize-1), 0.5*(self.ysize-1)), "Gaussian")

        # Get x and y mean
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        return x_mean, y_mean

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

    def save(self, path=None):

        """
        This function ...
        :param path:
        :param path:
        :return:
        """

        if path is None: path = self.path
        self.saveto(path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        # Set extra header info
        extra_header_info = dict()
        extra_header_info["PREPARED"] = self.prepared
        if self.from_filter is not None: extra_header_info["FRMFLTR"] = str(self.from_filter)
        if self.to_filter is not None: extra_header_info["TOFLTR"] = str(self.to_filter)

        # Call the save function of the base class
        super(ConvolutionKernel, self).saveto(path, extra_header_info=extra_header_info)

# -----------------------------------------------------------------

def get_fwhm(kernel_path):

    """
    Thisf untion ...
    :param path:
    :return:
    """

    from ..tools import headers
    from .fits import get_header

    # Get the header
    header = get_header(kernel_path)

    # Get and return the FWHM
    fwhm = headers.get_fwhm(header)
    return fwhm

# -----------------------------------------------------------------
