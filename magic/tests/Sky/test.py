#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from itertools import product

# Import astronomical modules
from photutils.datasets import make_100gaussians_image
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance, mad_std
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources, detect_threshold
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
#from photutils import Background2D, SigmaClip, MedianBackground # for 0.3.1
from photutils.background import Background
from astropy.utils import lazyproperty
from astropy.modeling.models import Sersic2D

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame
from pts.magic.core.mask import Mask

# -----------------------------------------------------------------

# https://photutils.readthedocs.io/en/stable/photutils/background.html

# -----------------------------------------------------------------

description = "Test the sky subtraction"

# -----------------------------------------------------------------

class SkyTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(SkyTest, self).__init__(config, interactive)

        # The sources frame
        self.original_frame = None

        # The frame
        self.frame = None

        # The sources mask
        self.sources_mask = None

        # The rotation mask
        self.rotate_mask = None

        # The real sky frame
        self.real_sky = None

        # Sky reference estimation
        self.reference_sky = None

        # Sky estimated by PTS
        self.estimated_sky = None

        # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Generate the data
        self.generate_data()

        # Rotate the data
        self.rotate()

        # Statistics
        self.statistics()

        # Sigma clip
        self.sigma_clip()

        # Mask sources
        self.mask_sources()

        # Add the background
        self.add_background()

        # Reference
        self.reference()

        # Add the galaxy
        #self.add_galaxy()

        # Subtract
        self.subtract()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SkyTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def generate_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the data ...")

        # we load a synthetic image comprised of 100 sources with a Gaussian-distributed
        # background whose mean is 5 and standard deviation is 2:

        self.real_sky = 5.

        # Create the image
        data = make_100gaussians_image()
        self.frame = Frame(data)
        self.original_frame = self.frame.copy()

        # Show
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.frame, norm=norm, origin='lower', cmap='Greys_r')
        plt.show()

    # -----------------------------------------------------------------

    def rotate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rotating the data ...")

        # Rotate
        new_data = rotate(self.frame.data, -45.)
        self.frame = Frame(new_data)

        # Set rotate mask
        self.rotate_mask = (self.frame == 0)

    # -----------------------------------------------------------------

    def statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating statistics ...")

        flattened = np.ma.array(self.frame.data, mask=self.rotate_mask.data).compressed()

        print(flattened)

        median = np.median(flattened)
        biweight_loc = biweight_location(flattened)

        biweight_midvar = biweight_midvariance(flattened)
        median_absolute_deviation = mad_std(flattened)

        print("median", median)
        print("biweigth_loc", biweight_loc)
        print("biweight_midvar", biweight_midvar)
        print("median_absolute_deviation", median_absolute_deviation)

        median = np.median(self.original_frame)
        biweight_loc = biweight_location(self.original_frame)
        biweight_midvar = biweight_midvariance(self.original_frame)
        median_absolute_deviation = mad_std(self.original_frame)

        print("median", median)
        print("biweigth_loc", biweight_loc)
        print("biweight_midvar", biweight_midvar)
        print("median_absolute_deviation", median_absolute_deviation)

    # -----------------------------------------------------------------

    def sigma_clip(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Sigma-clipping ...")

        # Sigma clip
        mean, median, std = sigma_clipped_stats(self.frame, sigma=3.0, iters=5, mask=self.rotate_mask)

        print("sigma-clip mean:", mean)
        print("sigma-clip median:", median)
        print("sigma-clip std:", std)

    # -----------------------------------------------------------------

    def mask_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Masking sources ...")

        # Create sources mask
        mask = make_source_mask(self.frame, snr=2, npixels=5, dilate_size=11, mask=self.rotate_mask)
        self.sources_mask = Mask(mask)

        # Statistics
        mean, median, std = sigma_clipped_stats(self.frame.data, sigma=3.0, mask=self.total_mask.data, iters=5)

        print("sigma-clip mean after source masking:", mean)
        print("sigma-clip median after source masking:", median)
        print("sigma_clip std after source masking:", std)

    # -----------------------------------------------------------------

    def add_background(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding background ...")

        ny, nx = self.frame.shape
        y, x = np.mgrid[:ny, :nx]

        # Add the sky
        self.real_sky += Frame(x * y / 5000.)

        # Construct the frame by adding background
        self.frame = self.frame + self.real_sky

        # Set padded pixels to zero again
        self.frame[self.rotate_mask] = 0.0

        # Show
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.frame, norm=norm, origin='lower', cmap='Greys_r')
        plt.show()

    # -----------------------------------------------------------------

    @lazyproperty
    def total_mask(self):

        """
        This function ...
        :return:
        """

        return self.sources_mask + self.rotate_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def real_subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.real_sky

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.reference_sky

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.estimated_sky

    # -----------------------------------------------------------------

    def add_galaxy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding smooth galaxy source ...")

        effective_radius = 20.

        axial_ratio = 2.
        angle_deg = 36.

        # Produce guess values
        initial_sersic_amplitude = 100.
        initial_sersic_r_eff = effective_radius / 10.0
        initial_sersic_n = 1.5
        initial_sersic_x_0 = 40.
        initial_sersic_y_0 = 123.
        initial_sersic_ellip = (axial_ratio - 1.0) / axial_ratio
        initial_sersic_theta = np.deg2rad(angle_deg)

        # Produce sersic model from guess parameters, for time trials
        sersic_x, sersic_y = np.meshgrid(np.arange(self.frame.xsize), np.arange(self.frame.ysize))
        sersic_model = Sersic2D(amplitude=initial_sersic_amplitude, r_eff=initial_sersic_r_eff,
                                                        n=initial_sersic_n, x_0=initial_sersic_x_0,
                                                        y_0=initial_sersic_y_0, ellip=initial_sersic_ellip,
                                                        theta=initial_sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

        self.galaxy = Frame(sersic_map)

        # Plot
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.galaxy, norm=norm, origin='lower', cmap='Greys_r')
        plt.show()

    # -----------------------------------------------------------------

    def reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating background with photutils ...")

        #sigma_clip = SigmaClip(sigma=3., iters=10)
        #bkg_estimator = MedianBackground()
        #bkg = Background2D(data2, (50, 50), filter_size=(3, 3), sigma_clip = sigma_clip, bkg_estimator = bkg_estimator)

        box_shape = (50,50)
        filter_size = (3,3)

        # Estimate the background
        bkg = Background(self.frame, box_shape, filter_shape=filter_size,
                         filter_threshold=None, mask=self.total_mask, method='sextractor',
                         backfunc=None, interp_order=3, sigclip_sigma=3.,
                         sigclip_iters=10)

        print("median background", bkg.background_median)

        print("rms background", bkg.background_rms_median)

        plt.imshow(bkg.background, origin='lower', cmap='Greys_r')

        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.frame - bkg.background, norm=norm, origin='lower', cmap = 'Greys_r')
        plt.show()

        # Set the sky
        self.reference_sky = bkg.background * ~self.total_mask.data

        # Show sky
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.reference_sky, norm=norm, origin='lower', cmap='Greys_r')
        plt.show()

        # Show subtracted frame
        plt.imshow(self.reference_subtracted, norm=norm, origin='lower', cmap='Greys_r')
        plt.show()


        # Will not work
        #mesh_yidx, mesh_xidx = np.unravel_index(self.mesh_idx, self._mesh_shape)
        #x, y, yx, data_coords = calc_coordinates(bkg, self.frame.shape)

        # Plot meshes
        #plt.imshow(self.frame, origin='lower', cmap='Greys_r', norm=norm)
        #plot_meshes(bkg, outlines=True, color='#1f77b4', x=x, y=y)
        #plt.show()

        # WE NEED NEWER PHOTUTILS TO MAKE THIS WORK

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky ...")

        # Settings
        settings = dict()

        # Input
        input_dict = dict()

        # Set input
        input_dict["frame"] = self.frame
        #input_dict["principal_shape"] =
        input_dict["sources_mask"] = self.sources_mask
        input_dict["extra_mask"] = self.rotate_mask

        # Create command
        command = Command("subtract_sky", "subtract the sky from an artificially created image", settings, input_dict)

        self.subtractor = self.run_command(command)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, norm=norm, origin='lower', cmap='Greys_r')

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------

def make_source_mask(data, snr, npixels, mask=None, mask_value=None,
                     filter_fwhm=None, filter_size=3, filter_kernel=None,
                     sigclip_sigma=3.0, sigclip_iters=5, dilate_size=11):
    """
    Make a source mask using source segmentation and binary dilation.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    snr : float
        The signal-to-noise ratio per pixel above the ``background`` for
        which to consider a pixel as possibly being part of a source.

    npixels : int
        The number of connected pixels, each greater than ``threshold``,
        that an object must have to be detected.  ``npixels`` must be a
        positive integer.

    mask : array_like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are ignored when computing the image background
        statistics.

    mask_value : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image background statistics.  ``mask_value`` will
        be ignored if ``mask`` is input.

    filter_fwhm : float, optional
        The full-width at half-maximum (FWHM) of the Gaussian kernel to
        filter the image before thresholding.  ``filter_fwhm`` and
        ``filter_size`` are ignored if ``filter_kernel`` is defined.

    filter_size : float, optional
        The size of the square Gaussian kernel image.  Used only if
        ``filter_fwhm`` is defined.  ``filter_fwhm`` and ``filter_size``
        are ignored if ``filter_kernel`` is defined.

    filter_kernel : array-like (2D) or `~astropy.convolution.Kernel2D`, optional
        The 2D array of the kernel used to filter the image before
        thresholding.  Filtering the image will smooth the noise and
        maximize detectability of objects with a shape similar to the
        kernel.  ``filter_kernel`` overrides ``filter_fwhm`` and
        ``filter_size``.

    sigclip_sigma : float, optional
        The number of standard deviations to use as the clipping limit
        when calculating the image background statistics.

    sigclip_iters : int, optional
       The number of iterations to perform sigma clipping, or `None` to
       clip until convergence is achieved (i.e., continue until the last
       iteration clips nothing) when calculating the image background
       statistics.

    dilate_size : int, optional
        The size of the square array used to dilate the segmentation
        image.

    Returns
    -------
    mask : 2D `~numpy.ndarray`, bool
        A 2D boolean image containing the source mask.
    """

    from scipy import ndimage

    threshold = detect_threshold(data, snr, background=None, error=None,
                                 mask=mask, mask_value=None,
                                 sigclip_sigma=sigclip_sigma,
                                 sigclip_iters=sigclip_iters)

    kernel = None
    if filter_kernel is not None:
        kernel = filter_kernel
    if filter_fwhm is not None:
        sigma = filter_fwhm * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma, x_size=filter_size,
                                  y_size=filter_size)
    if kernel is not None:
        kernel.normalize()

    segm = detect_sources(data, threshold, npixels, filter_kernel=kernel)

    selem = np.ones((dilate_size, dilate_size))
    return ndimage.binary_dilation(segm.data.astype(np.bool), selem)

# -----------------------------------------------------------------

def calc_coordinates(background, shape):

    """
    This function ...
    :param background:
    :param shape:
    :return:
    """

    """
    Calculate the coordinates to use when calling an interpolator.
    These are needed for `Background2D` and `BackgroundIDW2D`.
    Regular-grid interpolators require a 2D array of values.  Some
    require a 2D meshgrid of x and y.  Other require a strictly
    increasing 1D array of the x and y ranges.
    """

    # the position coordinates used to initialize an interpolation
    y = (background.mesh_yidx * background.box_size[0] +
              (background.box_size[0] - 1) / 2.)
    x = (background.mesh_xidx * background.box_size[1] +
              (background.box_size[1] - 1) / 2.)
    yx = np.column_stack([y, x])

    # the position coordinates used when calling an interpolator
    nx, ny = shape
    data_coords = np.array(list(product(range(ny), range(nx))))

    return x, y, yx, data_coords

# -----------------------------------------------------------------

def plot_meshes(bkg, ax=None, marker='+', color='blue', outlines=False, **kwargs):

    """
    Plot the low-resolution mesh boxes on a matplotlib Axes
    instance.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes` instance, optional
        If `None`, then the current ``Axes`` instance is used.

    marker : str, optional
        The marker to use to mark the center of the boxes.  Default
        is '+'.

    color : str, optional
        The color for the markers and the box outlines.  Default is
        'blue'.

    outlines : bool, optional
        Whether or not to plot the box outlines in addition to the
        box centers.

    kwargs
        Any keyword arguments accepted by
        `matplotlib.patches.Patch`.  Used only if ``outlines`` is
        True.
    """

    x = kwargs.pop("x")
    y = kwargs.pop("y")

    kwargs['color'] = color
    if ax is None:
        ax = plt.gca()
    ax.scatter(x, y, marker=marker, color=color)
    if outlines:
        from photutils.aperture_core import RectangularAperture
        xy = np.column_stack([x, y])
        apers = RectangularAperture(xy, bkg.box_shape[1], bkg.box_shape[0], 0.)
        apers.plot(ax=ax, **kwargs)
    return

# -----------------------------------------------------------------
