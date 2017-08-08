#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.get_fwhm Calculate the FWHM.
#  From: https://github.com/fred3m/astropyp/blob/master/astropyp/phot/psf.py
#  Has to be merged with the PTS find_sources and extract_sources codebase

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import traceback, functools
import numpy as np
import logging
from collections import OrderedDict

# Import astronomical modules
from astropy.io import fits
from astropy.extern import six
from astropy.nddata.utils import extract_array, subpixel_indices
from astropy.coordinates import SkyCoord
import astropy.units as apu
from astropy import table
from astropy.modeling import Fittable2DModel
from astropy.modeling.parameters import Parameter
from astropy.modeling.fitting import LevMarLSQFitter

# Import the relevant PTS classes and modules
from pts.core.basics.log import log

# -----------------------------------------------------------------

def trace_unhandled_exceptions(func):

    """
    Wrapper for multiprocessing pool functions. By default when a
    multiprocessing pool fails the process exists but with no
    traceback that describes the error that crashed the script.
    By using this function wrapper on multiprocess pool functions
    users get a traceback for each process that crashes.

    Use::

        @trace_unhandled_exceptions
        def my_multiprocessing_func(x,y):
            return x/y
    """

    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            print('Exception in ' + func.__name__)
            traceback.print_exc()
            raise Exception("Error in multiprocessing")

    return wrapped_func

# -----------------------------------------------------------------

psf_flags = OrderedDict([
    (128, 'Bad Pixel'),  # Bad pixels
    (64, 'Edge'),  # Edge is closer than aperture radius
    (32, '<not used>'),
    (16, '<not used>'),
    (8, 'Low Peak'),  # Amplitude is below minimum
    (4, 'Low Flux'),  # Flux is below minimum
    (2, 'Elliptical'),  # Source is elliptical, not circular
    (1, 'Crowded')  # nearby source, requires group photometry
])

def select_psf_sources(img_data, catalog, aper_radius=None,
                       min_flux=None, min_amplitude=None, min_dist=None,
                       max_ratio=None, edge_dist=None, verbose=True, units='deg',
                       badpix_flags=[], flag_max=0):
    """
    From an image and a list of sources, generate a list of high SNR
    sources with circular profiles that are likely to be good sources
    to generate the psf.

    Parameters
    ----------
    img_data: array-like
        The image data
    catalog: `~astropyp.catalog.ImageCatalog`
        Catalog of sources. The properties
        x, y, ra, dec, aper_flux, peak, a, b
        are all required to be mapped to the appropriate columns in the
        catalog.sources table.
    aper_radius: integer
        Radius to use to calculate psf flux.
    min_flux: integer, optional
        Minimum flux to be considered an valid source. If no
        ``min_flux`` is given then no cuts will be made based on
        the total flux
    min_aplitude: integer, optional
        Minimum peak flux of a pixel in for a psf source. If no
        ``min_amplitude`` is given then no cuts will be made based
        on maximum flux in a pixel.
    min_dist: integer, optional
        Minimum distance (in pixels) to the nearest neighbor
        for a PSF source. If no ``min_dist`` is specified the
        function will set ``min_dist=3*aperture_radius``.
    max_ratio: float
        Maximum ratio of a sources two axes to be considered
        sufficiently elliptical. To calculate max_ratio,
        ``a`` and ``b`` must be specified. The default value
        is the median(a/b)+stddev(a/b) for all of the sources
    edge_dist: integer, optional
        Maximum distance (in pixels) from the center of a source to the
        edge of the image. This defaults to the ``aperture_radius``+1.
    verbose: bool, optional
        Whether or not to log the results of the various cuts.
        *Default value is True*
    units: string or `~astropy.units.unit`
        Units of the ra and dec columns. *Default value is 'deg'*
    badpix_flags: list of strings or list of array-like, optional
        A list of column names in ``sources`` or a list of arrays
        with the badpixel flags for each source.
        *Default is an empty list*
    flag_max: integer, optional
        Maximum value of a bad pixel flag to be considered good.
        *Default value is 0*

    Returns
    -------
    psf_idx: `~nump.ndarray`
        Boolean array that matches psf sources selected by the function
    flags: `~astropy.table.Table`
        Table of flags for all of the sources
    """
    from scipy import spatial

    rows = catalog.shape[0]

    # Only select high signal to noise sources
    min_flux_flag = ~np.isnan(catalog.aper_flux)
    if min_flux is not None:
        min_flux_flag = min_flux_flag & (catalog.aper_flux > min_flux)
    if min_amplitude is not None:
        min_amp_flag = catalog.peak > min_amplitude

    # Eliminate Bad Pixel Sources
    bad_src_flag = np.ones((rows,), dtype=bool)
    for flags in badpix_flags:
        if isinstance(flags, six.string_types):
            flags = catalog[flags] <= flag_max
        bad_src_flag = bad_src_flag & flags

    # Cut out elliptical sources (likely galaxies)
    if catalog.a is not None and catalog.b is not None:
        a_b = catalog.b / catalog.a
        if max_ratio is None:
            max_ratio = np.median(a_b) - np.std(a_b) * 1.5
        elliptical_flag = catalog.b / catalog.a > max_ratio

    # Cut out all sources with a nearby neighbor
    if min_dist is None:
        min_dist = 3 * aper_radius
    KDTree = spatial.cKDTree
    pts = zip(catalog.x, catalog.y)
    kdt = KDTree(pts)
    d2, idx = kdt.query(pts, 2)
    distance_flag = d2[:, 1] > min_dist

    # Cut out all source near the edge
    if edge_dist is None:
        edge_dist = aper_radius + 1
    edge_flag = ((catalog.x > edge_dist) & (catalog.y > edge_dist) &
                 (catalog.x < img_data.shape[1] - edge_dist) &
                 (catalog.y < img_data.shape[0] - edge_dist))

    # Apply all of the cuts to the table of sources
    psf_idx = (
        min_flux_flag &
        min_amp_flag &
        bad_src_flag &
        elliptical_flag &
        distance_flag &
        edge_flag)

    flags = (
        ~bad_src_flag * 128 +
        ~edge_flag * 64 +
        ~min_flux_flag * 4 +
        ~elliptical_flag * 2 +
        ~distance_flag)
    catalog.sources['pipeline_flags'] = flags

    # Combine the flags into a table that can be used later
    flags = table.Table(
        [min_flux_flag, min_amp_flag, bad_src_flag,
         elliptical_flag, distance_flag, edge_flag],
        names=('min_flux', 'min_amp', 'bad_pix', 'ellipse', 'dist', 'edge')
    )

    # If verbose, print information about the cuts
    if verbose:
        # level = log.getEffectiveLevel()
        # log.setLevel(logging.INFO)
        log.info('Total sources: {0}'.format(rows)),
        log.info('Sources with low flux: {0}'.format(
            rows - np.sum(min_flux_flag)))
        log.info('Sources with low amplitude: {0}'.format(
            rows - np.sum(min_amp_flag)))
        log.info('Sources with bad pixels: {0}'.format(
            rows - np.sum(bad_src_flag)))
        log.info('Elliptical sources: {0}'.format(
            rows - np.sum(elliptical_flag)))
        log.info('Source with close neighbors: {0}'.format(
            rows - np.sum(distance_flag)))
        log.info('Sources near an edge: {0}'.format(
            rows - np.sum(edge_flag)))
        log.info('Sources after cuts: {0}'.format(np.sum(psf_idx))),
        # log.setLevel(level)
    return psf_idx, flags


def build_psf(img_data, aper_radius, sources=None, x='x', y='y',
              subsampling=5, combine_mode='median', offset_buffer=3):
    """
    Build a subsampled psf by stacking a list of psf sources and
    a its source image.

    Parameters
    ----------
    img_data: array-like
        Image containing the sources
    aper_radius: integer
        Radius of the aperture to use for the psf source. This must be an
        odd number.
    sources: `~astropy.table.Table`, optional
        Table with x,y positions of psf sources. If sources is passed to the
        function, ``x`` and ``y`` should be the names of the x and y
        coordinates in the table.
    x: string or array-like, optional
        ``x`` can either be a string, the name of the x column in ``sources``,
        or an array of x coordinates of psf sources. *Defaults to 'x'*
    y: string or array-like, optional
        ``y`` can either be a string, the name of the y column in ``sources``,
        or an array of y coordinates of psf sources. *Defaults to 'y'*
    subsampling: integer, optional
        Number of subdivided pixel in the psf for each image pixel.
        *Defaults value is 5*
    combine_mode: string, optional
        Method to combine psf sources to build psf. Can be either
        "median" or "mean". *Defaults to "median"*
    offset_buffer: integer
        Error in image pixels for the center of a given source.
        This allows the code some leeway to re-center an image
        based on its subpixels. *Default is 3*
    """
    from astropyp import utils

    if isinstance(x, six.string_types):
        x = sources[x]
    if isinstance(y, six.string_types):
        y = sources[y]
    if hasattr(x, 'shape'):
        rows = x.shape[0]
    else:
        rows = len(x)

    if combine_mode == 'median':
        combine = np.ma.median
    elif combine_mode == 'mean':
        combine = np.ma.mean
    else:
        raise Exception(
            "Invalid combine_mode, you must choose from ['median','mean']")

    # Load a patch centered on each PSF source and stack them to build the
    # final PSF
    src_width = 2 * aper_radius + 1
    patches = []
    obj_shape = (src_width, src_width)
    for n in range(rows):
        patch, X, Y, new_yx = utils.misc.get_subpixel_patch(
            img_data, (y[n], x[n]), obj_shape,
            max_offset=offset_buffer, subsampling=subsampling,
            normalize=True)
        if patch is not None:
            patches.append(patch)
        else:
            log.info("{0} contains NaN values, skipping".format(
                (new_yx[1], new_yx[0])))
    psf = combine(np.dstack(patches), axis=2)

    # Add a mask to hide values outside the aperture radius
    radius = aper_radius * subsampling
    ys, xs = psf.shape
    yc, xc = (np.array([ys, xs]) - 1) / 2
    y, x = np.ogrid[-yc:ys - yc, -xc:xs - xc]
    mask = x ** 2 + y ** 2 <= radius ** 2
    circle_mask = np.ones(psf.shape)
    circle_mask[mask] = 0
    circle_mask
    psf_array = np.ma.array(psf, mask=circle_mask)
    return psf_array


def perform_psf_photometry(data, catalog, psf, separation=None,
                           verbose=False, fit_position=True, pos_range=0, indices=None,
                           kd_tree=None, exptime=None, pool_size=None):
    """
    Perform PSF photometry on all of the sources in the catalog,
    or if indices is specified, a subset of sources.

    Parameters
    ----------
    data: `~numpy.ndarray`
        Image containing the sources
    catalog: `~astropyp.catalog.Catalog`
        Catalog of sources to use. This must contain the
        x,y keys.
    psf: `SinglePSF`
        PSF to use for the fit
    separation: float, optional
        Separation (in pixels) for members to be considered
        part of the same group *Default=1.5*psf width*
    verbose: bool, optional
        Whether or not to show info about the fit progress.
        *Default=False*
    fit_position: bool, optional
        Whether or not to fit the position along with the
        amplitude of each source. *Default=True*
    pos_range: int, optional
        Maximum number of pixels (in image pixels) that
        a sources position can be changed. If ``pos_range=0``
        no bounds will be set. *Default=0*
    indices: `~numpy.ndarray` or string, optional
        Indices for sources to calculate PSF photometry.
        It is often advantageous to remove sources with
        bad pixels and sublinear flux to save processing time.
        All sources not included in indices will have their
        psf flux set to NaN. This can either be an array of
        indices for the positions in self.catalog or
        the name of a saved index in self.indices
    kd_tree: `scipy.spatial.KDTree`, optional
        KD tree containing all of the sources in the catalog
    exptime: float, optional
        Exposure time of the image
        (needed to calculate the psf magnitude).
    """
    import astropyp.catalog

    # Get the positions and estimated amplitudes of
    # the sources to fit
    if indices is not None:
        positions = zip(catalog.x[indices],
                        catalog.y[indices])
        amplitudes = catalog.peak[indices]
    else:
        positions = zip(catalog.x, catalog.y)
        amplitudes = catalog.peak

    src_count = len(positions)

    src_indices = np.arange(0, len(catalog.sources), 1)
    all_positions = np.array(zip(catalog.x, catalog.y))
    all_amplitudes = catalog.peak
    total_sources = len(all_amplitudes)
    src_psfs = []

    psf_flux = np.ones((total_sources,)) * np.nan
    psf_flux_err = np.ones((total_sources,)) * np.nan
    psf_x = np.ones((total_sources,)) * np.nan
    psf_y = np.ones((total_sources,)) * np.nan
    new_amplitudes = np.ones((total_sources,)) * np.nan

    # Find nearest neighbors to avoid contamination by nearby sources
    if separation is not None:
        if kd_tree is None:
            from scipy import spatial
            KDTree = spatial.cKDTree
            kd_tree = KDTree(all_positions)
        idx, nidx = astropyp.catalog.find_neighbors(separation,
                                                    kd_tree=kd_tree)
    # Set the number of processors to use
    if pool_size is None:
        import multiprocessing
        pool_size = multiprocessing.cpu_count()

    # Fit each source to the PSF, calcualte its flux, and its new
    # position
    if pool_size == 1:
        for n in range(len(positions)):
            if verbose:
                # level = log.getEffectiveLevel()
                log.setLevel(logging.INFO)
                log.info("Fitting {0}".format(group_id))
                # log.setLevel(level)
            src_idx = src_indices[indices][n]
            if separation is not None:
                n_indices = nidx[idx == src_idx]
                neighbor_positions = all_positions[n_indices]
                neighbor_amplitudes = all_amplitudes[n_indices]
            else:
                neighbor_positions = []
                neighbor_amplitudes = []
            src_psf = SinglePSF(psf._psf_array,
                                amplitudes[n], positions[n][0], positions[n][1],
                                subsampling=psf._subsampling,
                                neighbor_positions=neighbor_positions,
                                neighbor_amplitudes=neighbor_amplitudes)
            psf_flux[src_idx], psf_flux_err[src_idx], residual, new_pos = \
                src_psf.fit(data)
            new_amplitudes[src_idx] = src_psf.amplitude.value
            psf_x[src_idx], psf_y[src_idx] = new_pos
    else:
        import multiprocessing
        if pool_size is None:
            pool_size = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(
            processes=pool_size,
            initializer=_init_psf_phot_process,
            initargs=(data, all_positions, all_amplitudes,
                      idx, nidx, psf, amplitudes, positions))
        pool_args = []
        for n in range(len(positions)):
            src_idx = src_indices[indices][n]
            pool_args.append((n, src_idx))
        result = pool.map(_psf_phot_worker, pool_args)
        pool.close()
        pool.join()

        result = zip(*result)
        psf_indices = np.where(indices)
        psf_flux[psf_indices] = result[0]
        psf_flux_err[psf_indices] = result[1]
        new_amplitudes[psf_indices] = result[2]
        psf_x[psf_indices] = result[3]
        psf_y[psf_indices] = result[4]

    # Save the psf derived quantities in the catalog
    # Ignore divide by zero errors that occur when sources
    # have zero psf flux (i.e. bad sources)
    np_err = np.geterr()
    np.seterr(divide='ignore')
    catalog.sources['psf_flux'] = psf_flux
    catalog.sources['psf_flux_err'] = psf_flux_err
    if exptime is not None:
        psf_mag = -2.5 * np.log10(psf_flux / exptime)
        catalog.sources['psf_mag'] = psf_mag
        catalog.sources['psf_mag_err'] = psf_flux_err / psf_flux
        catalog.sources['psf_x'] = psf_x
        catalog.sources['psf_y'] = psf_y
    np.seterr(**np_err)
    return catalog, src_psfs, kd_tree

# -----------------------------------------------------------------

@trace_unhandled_exceptions
def _init_psf_phot_process(data, all_positions, all_amplitudes,
                           idx, nidx, psf, amplitudes, positions):
    """
    Global variables for psf photometry processes
    """
    import multiprocessing
    global gbl_data
    global gbl_all_positions
    global gbl_all_amplitudes
    global gbl_idx
    global gbl_nidx
    global gbl_psf
    global gbl_amplitudes
    global gbl_positions

    gbl_data = data
    gbl_all_positions = all_positions
    gbl_all_amplitudes = all_amplitudes
    gbl_idx = idx
    gbl_nidx = nidx
    gbl_psf = psf
    gbl_amplitudes = amplitudes
    gbl_positions = positions
    log.info("Initializing process {0}".format(
        multiprocessing.current_process().name))

# -----------------------------------------------------------------

@trace_unhandled_exceptions
def _psf_phot_worker(args):
    """
    Worker to perform psf phot on individual images
    """
    n, src_idx = args
    if gbl_nidx is not None:
        n_indices = gbl_nidx[gbl_idx == src_idx]
        neighbor_positions = gbl_all_positions[n_indices]
        neighbor_amplitudes = gbl_all_amplitudes[n_indices]
    else:
        neighbor_positions = []
        neighbor_amplitudes = []

    src_psf = SinglePSF(gbl_psf._psf_array,
                        gbl_amplitudes[n], gbl_positions[n][0], gbl_positions[n][1],
                        subsampling=gbl_psf._subsampling,
                        neighbor_positions=neighbor_positions,
                        neighbor_amplitudes=neighbor_amplitudes)
    psf_flux, psf_flux_err, residual, new_pos = \
        src_psf.fit(gbl_data)
    amplitudes = src_psf.amplitude.value
    psf_x, psf_y = new_pos
    return psf_flux, psf_flux_err, amplitudes, psf_x, psf_y

# -----------------------------------------------------------------

class SinglePSF(Fittable2DModel):
    """
    A discrete PSF model for a single source. Before creating a
    SinglePSF model it is necessary to run `build_psf` create a
    ``psf_array`` that is used to calibrate the source by modifying
    its amplitude and (optionally) its position.

    Parameters
    ----------
    psf_array: `numpy.ndarray`
        PSF array using ``subsampling`` subpixels for each pixel in the image
    amplitude: float
        Amplitude of the psf function
    x0: float
        X coordinate in the image
    y0: float
        Y coordinate in the image
    subsampling: int, optional
        Number of subpixels used in the ``psf_array`` for each pixel in the
        image
    recenter: bool
        Whether or not to recenter the data patch on the maximum
        value. *Default=True*
    fit_position: bool
        Whether or not to fit the positions and the amplitude or only the
        amplitude. *Default=True, fit positions and amplitude*
    pos_range: float
        +- bounds for the position (if ``fit_position=True``).
        If ``pos_range=0`` then no bounds are used and the
        x0 and y0 parameters are free. *Default=0*
    """
    amplitude = Parameter('amplitude')
    x0 = Parameter('x0')
    y0 = Parameter('y0')
    _param_names = ()

    # -----------------------------------------------------------------

    def __init__(self, psf_array, amplitude, x0, y0, subsampling=5,
                 recenter=True, dx=None, dy=None, amp_err=None,
                 neighbor_positions=[], neighbor_amplitudes=[]):
        # Set parameters defined by the psf to None so that they show
        # up as class attributes. self.set_psf_array will set all
        # of these parameters automatically
        self._subpixel_width = None
        self._width = None
        self._radius = None
        self._subsampling = None
        self._psf_array = None

        self.recenter = recenter
        self.fitter = LevMarLSQFitter()

        # New CrowdedPSF attributes
        if len(neighbor_positions) != len(neighbor_amplitudes):
            raise Exception("neighbor_positions and neighbors amplitudes "
                            "lengths {0},{1} are nor equal".format(
                len(neighbor_positions), len(neighbor_amplitudes)))
        self.src_count = len(neighbor_positions)
        self.x0_names = ['nx_{0}'.format(n)
                         for n in range(self.src_count)]
        self.y0_names = ['ny_{0}'.format(n)
                         for n in range(self.src_count)]
        self.amp_names = ['namp{0}'.format(n)
                          for n in range(self.src_count)]
        self._param_names = tuple(
            ['amplitude', 'x0', 'y0'] + self.x0_names + self.y0_names + self.amp_names)

        super(SinglePSF, self).__init__(n_models=1, x0=x0, y0=y0,
                                        amplitude=amplitude)

        # Choose whether or not the position of the PSF can be moved
        if False:  # fit_position:
            self.x0.fixed = False
            self.y0.fixed = False
            if pos_range > 0:
                self.x0.bounds = (self.x0 - pos_range, self.x0 + pos_range)
                self.y0.bounds = (self.y0 - pos_range, self.y0 + pos_range)
        else:
            pass
            # self.x0.fixed = True
            # self.y0.fixed = True
        self.set_psf_array(psf_array, subsampling=subsampling)
        amplitudes = neighbor_amplitudes
        # Set parameters for neighboring sources (if any exist)
        if self.src_count > 0:
            pos_array = np.array(neighbor_positions)
            x = pos_array[:, 0]
            y = pos_array[:, 1]
            kwargs = OrderedDict()

            for n in range(self.src_count):
                x0_name = self.x0_names[n]
                y0_name = self.y0_names[n]
                amp_name = self.amp_names[n]

                setattr(self, x0_name, x[n])
                setattr(self, y0_name, y[n])
                setattr(self, amp_name, amplitudes[n])

                if dx is not None:
                    x0 = getattr(self, x0_name)
                    x0.bounds = (x[n] - dx, x[n] + dx)
                if dy is not None:
                    y0 = getattr(self, y0_name)
                    y0.bounds = (y[n] - dy, y[n] + dy)
                if amp_err is not None:
                    amp = getattr(self, amp_name)
                    amp.bounds = (
                        amplitudes[n] * (1 - amp_err),
                        amplitudes[n] * (1 + amp_err))

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        Shape of the PSF image.
        """

        return self._psf_array.shape

    # -----------------------------------------------------------------

    def get_flux(self, amplitude=None):
        """
        PSF Flux of the source
        """
        if amplitude is None:
            amplitude = self.amplitude
        flux = np.sum(self._psf_array * amplitude) / self._subsampling ** 2
        return flux

    # -----------------------------------------------------------------

    def evaluate(self, X, Y, amplitude, x0, y0, *args):

        """
        Evaluate the SinglePSF model.
        """

        x = X[0, :]
        y = Y[:, 0]
        result = amplitude * self._psf_func(y - y0, x - x0)
        # If the source has any neighbors, add their
        # fluxes as well
        nx, ny, namp = self._parse_eval_args(args)
        for n in range(self.src_count):
            result += namp[n] * self._psf_func(y - ny[n], x - nx[n])
        result[self._psf_array.mask] = 0
        return result

    # -----------------------------------------------------------------

    def _parse_eval_args(self, args):

        """
        Separate the x0, y0, and amplitudes into lists
        """

        x0 = args[:self.src_count]
        y0 = args[self.src_count:2 * self.src_count]
        amp = args[2 * self.src_count:]
        return x0, y0, amp

    # -----------------------------------------------------------------

    @property
    def param_names(self):
        return self._param_names

    # -----------------------------------------------------------------

    def __getattr__(self, attr):

        if self._param_names and attr in self._param_names:
            return Parameter(attr, default=0.0, model=self)
        raise AttributeError(attr)

    # -----------------------------------------------------------------

    def __setattr__(self, attr, value):

        if attr[0] != '_' and self._param_names and attr in self._param_names:
            param = Parameter(attr, default=0.0, model=self)
            # This is a little hackish, but we can actually reuse the
            # Parameter.__set__ method here
            param.__set__(self, value)
        else: super(SinglePSF, self).__setattr__(attr, value)

    # -----------------------------------------------------------------

    def fit(self, img_data=None, fit_position=True, pos_range=0, indices=None,
            patch=None, X=None, Y=None):
        """
        Fit the PSF to the data.

        Parameters
        ----------
        img_data: array-like, optional
            Image data containing the source to fit. Either img_data must
            be specified or patch.
        fit_position: bool, optional
            Whether or not to fit the position. If ``fit_position=False``
            only the amplitude will be fit
        pos_range: float, optional
            Maximum distance that the position is allowed to shift
            from ``x0,y0`` initial. This is most useful when fitting
            a group of sources where the code might try to significantly
            move one of the sources for the fit. The default range
            is 0, which does not set any bounds at all.
        indices: list of arrays, optional
            ``indices`` is a list that contains the X and Y indices
            to the img_data provided. The default is None, which
            uses the image scale and size of the psf to set the indices.
        patch: subsampled image patch, optional
            Subset of data to use for the fit. If not included (default)
            then img_data must be specified. This must be the same shape
            as the psf.
        X,Y: array-like, optional
            Positions for the coordinates on the x,y axes associated with
            the patch. These must be given if patch is not None.
        """
        import astropyp.utils
        # Extract sub array with data of interest
        position = (self.y0.value, self.x0.value)
        if patch is None:
            patch, X, Y, new_pos = astropyp.utils.misc.get_subpixel_patch(
                img_data, position, (self._width, self._width),
                subsampling=self._subsampling,
                window_sampling=300, order=5,  # Need to open these options up to the class init
                normalize=False)
            # If the source was too close to an edge a patch
            # cannot be loaded, source the source cannot be fit
            if patch is None:
                return np.nan, np.nan, np.nan, (np.nan, np.nan)
            self.y0.value, self.x0.value = new_pos
        else:
            new_pos = (self.y0.value, self.x0.value)
            x_radius = self.shape[1] / self._subsampling
            y_radius = self.shape[0] / self._subsampling
            x0 = self.x0.value
            y0 = self.y0.value
            X = np.linspace(x0 - x_radius, x0 + x_radius, self.shape[1])
            Y = np.linspace(y0 - y_radius, y0 + y_radius, self.shape[0])

        # Set the values outside the aperture to zero
        # so they will not affect the fit
        patch[self._psf_array.mask] = 0
        X, Y = np.meshgrid(X, Y)

        # Fit only if PSF is completely contained in the image and no NaN
        # values are present
        if (patch.shape == self.shape and
                not np.isnan(patch).any()):
            m = self.fitter(self, X, Y, patch)

            for p in self._param_names:
                setattr(self, p, getattr(m, p))

            residual = self.get_residual(X, Y, patch)
            self.psf_error = self.calculate_psf_error(residual)
            return (self.get_flux(), self.psf_error,
                    residual, (self.x0.value, self.y0.value))
        else:
            return (0, 0, patch, (self.x0.value, self.y0.value))

    # -----------------------------------------------------------------

    def get_residual(self, X, Y, patch):

        fit = self.__call__(X, Y)
        residual = patch - fit
        return residual

    # -----------------------------------------------------------------

    def calculate_psf_error(self, residual=None, X=None, Y=None, patch=None):
        """
        Given either a residual or a patch, calculate the error in PSF flux
        """
        if residual is None:
            if patch is not None and X is not None and Y is not None:
                residual = self.get_residual(X, Y, patch)
            else:
                raise Exception(
                    "You must either supply a residual or "
                    "X,Y,patch to calculate psf error")
        psf_error = np.abs(np.sum(residual / self._subsampling ** 2))
        return psf_error

    def set_psf_array(self, psf_array, subsampling=None):
        """
        If a change is made to the psf function, it must
        be updated here so that all of the derived parameters
        (_width, _radius, _psf_array) can be updated
        """
        try:
            from scipy import interpolate
        except ImportError:
            raise Exception(
                "You must have scipy installed to use the SinglePSF class")
        if subsampling is not None:
            self._subsampling = subsampling

        # Set the width and radius of the psf
        self._subpixel_width = max(psf_array.shape[0], psf_array.shape[1])
        self._width = self._subpixel_width / self._subsampling
        self._radius = (self._width - 1) / 2.
        self._psf_array = psf_array

        # Set the function to determine the psf (up to an amplitude)
        X = np.linspace(-self._radius, self._radius,
                        self._subpixel_width)
        Y = np.linspace(-self._radius, self._radius,
                        self._subpixel_width)
        self._psf_func = interpolate.RectBivariateSpline(Y, X, self._psf_array)

# -----------------------------------------------------------------

class GroupPSF:

    """
    This represents the PSFs of a group of sources. In general
    a `GroupPSF` should only be created by the `ImagePSF` class.
    """

    def __init__(self, group_id, psf, positions=None, amplitudes=None,
                 bounds=None, mask_img=True, show_plots=True, **kwargs):
        if isinstance(psf, GroupPSF):
            self.__dict__ = psf.__dict__.copy()
        else:
            self.group_id = group_id
            self.psf = psf.copy()
            self.positions = positions
            self.amplitudes = amplitudes
            self.mask_img = mask_img
            self.show_plots = show_plots
            self.bounds = bounds
            self.combined_psf = None
        for k, v in kwargs:
            setattr(self, k, v)

    def get_patch(self, data, positions=None, patch_bounds=None,
                  show_plots=None):
        """
        Given a list of positions, get the patch of the data that
        contains all of the positions and their PSF radii and mask
        out all of the other pixels (to prevent sources outside the
        group from polluting the fit).

        Parameters
        ----------
        data: ndarray
            Image array data
        positions: list or array (optional)
            List of positions to include in the patch. If no
            positions are passed the function will use
            ``GroupPSF.positions``.
        patch_bounds: list or array (optional)
            Boundaries of the data patch of the form
            [ymin,ymax,xmin,xmax]
        """
        from scipy import interpolate

        if positions is None:
            positions = np.array(self.positions)
        if show_plots is None:
            show_plots = self.show_plots

        radius = int((self.psf._width - 1) / 2)
        x = positions[:, 0]
        y = positions[:, 1]

        # Extract the patch of data that includes all of the sources
        # and their psf radii
        if patch_bounds is None:
            if self.bounds is None:
                xc = np.round(x).astype('int')
                yc = np.round(y).astype('int')
                self.bounds = [
                    max(min(yc) - radius, 0),  # ymin
                    min(max(yc) + radius + 1, data.shape[0]),  # ymax
                    max(min(xc) - radius, 0),  # xmin
                    min(max(xc) + radius + 1, data.shape[1])  # xmax
                ]
            patch_bounds = self.bounds
        ymin, ymax, xmin, xmax = patch_bounds
        X_img = np.arange(xmin, xmax + 1, 1)
        Y_img = np.arange(ymin, ymax + 1, 1)
        Z = data[np.ix_(Y_img, X_img)]
        X = np.arange(xmin, xmax + 1, 1 / self.psf._subsampling)
        Y = np.arange(ymin, ymax + 1, 1 / self.psf._subsampling)
        data_func = interpolate.RectBivariateSpline(Y_img, X_img, Z)
        sub_data = data_func(Y, X, Z)

        # If the group is large enough, sources not contained
        # in the group might be located in the same square patch,
        # so we mask out everything outside of the radius the
        # individual sources PSFs
        if self.mask_img:
            sub_data = np.ma.array(sub_data)
            mask = np.ones(data.shape, dtype='bool')
            mask_X = np.arange(data.shape[1])
            mask_Y = np.arange(data.shape[0])
            mask_X, mask_Y = np.meshgrid(mask_X, mask_Y)
            for xc, yc in zip(x, y):
                mask = mask & (
                    (mask_X - xc) ** 2 + (mask_Y - yc) ** 2 >= (radius) ** 2)
            # Expand the mask to be the same size as the data patch
            subsampling = self.psf._subsampling
            sub_data.mask = np.kron(mask[np.ix_(Y_img, X_img)],
                                    np.ones((subsampling, subsampling), dtype=int))
            sub_data = sub_data.filled(0)

        # Optionally plot the mask and data patch
        if show_plots:
            try:
                import matplotlib
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d.axes3d import Axes3D
            except ImportError:
                raise Exception(
                    "You must have matplotlib installed"
                    " to create plots")
            # Plot mask
            if self.mask_img:
                plt.imshow(mask[ymin:ymax, xmin:xmax])

            # Plot masked patch used for fit
            # Xplt = np.arange(0, sub_data.shape[1], 1)
            # Yplt = np.arange(0, sub_data.shape[0], 1)
            Xplt, Yplt = np.meshgrid(X, Y)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_wireframe(Xplt, Yplt, sub_data,
                              rstride=5, cstride=5)
            plt.show()
        return sub_data, X, Y

    # -----------------------------------------------------------------

    def fit(self, data, positions=None, psfs=None, amplitudes=None,
            patch_bounds=None, fitter=None):
        """
        Simultaneously fit all of the sources in a PSF group. This
        functions makes a copy of the PSF for each source and creates
        an astropy `astropy.models.CompoundModel` that is fit but the PSF's
        ``fitter`` function.

        Parameters
        ----------
        data: ndarray
            Image data array
        positions : List or array (optional)
            List of positions in pixel coordinates
            where to fit the PSF/PRF. Ex: [(0.0,0.0),(1.0,2.0), (10.3,-3.2)]
        psfs : list of PSFs (optional)
            It is possible to choose a different PSF for each model,
            in which case ``psfs`` should have the same length as positions
        amplitudes: list of floats, optional
            Amplitudes (peak values) for the sources. If amplitudes is
            none then the pixel value at the position is used
        patch_bounds: list or array (optional)
            Boundaries of the data patch of the form
            [ymin,ymax,xmin,xmax]
        fitter: `~astropy.modeling.Fitter`
            Fitting class to use to fit the data. *Default is self.fitter*
        """
        if positions is None:
            positions = np.array(self.positions)

        x = positions[:, 0]
        y = positions[:, 1]

        # If not estimated ampltiudes are given for the sources,
        # use the pixel value at their positions
        if amplitudes is None:
            amplitudes = data[y.astype(int), x.astype(int)]

        if len(positions) == 1:
            self.psf.x_0, self.psf.y_0 = positions[0]
            result = [self.psf.fit(data)]
        else:
            if psfs is None:
                psfs = [self.psf.copy() for p in range(len(positions))]
            if fitter is None:
                if self.psf is not None:
                    fitter = self.psf.fitter
                else:
                    fitter = psfs[0].fitter
            sub_data, X, Y = self.get_patch(data, positions, patch_bounds)
            X, Y = np.meshgrid(X, Y)

            # Create a CompoundModel that is a combination of the
            # individual PRF's
            combined_psf = None
            for x0, y0, single_psf, amplitude in zip(x, y, psfs, amplitudes):
                single_psf.x_0 = x0
                single_psf.y_0 = y0
                single_psf.amplitude = amplitude
                if combined_psf is None:
                    combined_psf = single_psf
                else:
                    combined_psf += single_psf
            # Fit the combined PRF
            self.combined_psf = fitter(combined_psf, X, Y, sub_data)

            # Return the list of fluxes for all of the sources in the group
            # and the combined PRF
            result = [getattr(self.combined_psf, 'amplitude_' + str(n)).value
                      for n in range(len(x))]
        return result

# -----------------------------------------------------------------
