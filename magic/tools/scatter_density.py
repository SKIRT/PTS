#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.scatter_density Contains the ScatterDensityAxis class and related stuff (adapted from mpl_scatter_density package).

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections import register_projection
from math import ceil, log10
from matplotlib.image import AxesImage, _interpd_, _image
from matplotlib.transforms import IdentityTransform, TransformedBbox, BboxTransformFrom, Bbox, Affine2D
import matplotlib.colors as mcolors

# Import other modules
from fast_histogram import histogram2d

# -----------------------------------------------------------------

class ScatterDensityAxes(plt.Axes):

    """
    This class ...
    """

    name = 'scatter_density'

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        plt.Axes.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    def scatter_density(self, x, y, dpi=72, downres_factor=4, color=None, cmap=None,
                        alpha=1.0, norm=None, **kwargs):
        """
        Make a density plot of the (x, y) scatter data.

        Parameters
        ----------
        x, y : iterable
            The data to plot
        dpi : int or `None`
            The number of dots per inch to include in the density map. To use
            the native resolution of the drawing device, set this to None.
        downres_factor : int
            For interactive devices, when panning, the density map will
            automatically be made at a lower resolution and including only a
            subset of the points. The new dpi of the figure when panning will
            then be dpi / downres_factor, and the number of elements in the
            arrays will be reduced by downres_factor**2.
        cmap : `matplotlib.colors.Colormap`
            The colormap to use for the density map.
        color : str or tuple
            The color to use for the density map. This can be any valid
            Matplotlib color. If specified, this takes precedence over the
            colormap.
        alpha : float
            Transparency of the density map
        norm : `matplotlib.colors.Normalize`
            The normalization class for the density map.
        """

        self.set_xlim(np.min(x), np.max(x))
        self.set_ylim(np.min(y), np.max(y))

        scatter = ScatterDensityArtist(self, x, y, dpi=dpi, downres_factor=downres_factor,
                                       color=color, cmap=cmap,
                                       alpha=alpha, norm=norm, **kwargs)
        self.add_artist(scatter)

        return scatter

# -----------------------------------------------------------------

register_projection(ScatterDensityAxes)

# -----------------------------------------------------------------

def make_cmap(color):

    r, g, b = mcolors.colorConverter.to_rgb(color)
    cdict = {'red': [(0.0, r, r),
                     (1.0, r, r)],
             'green': [(0.0, g, g),
                       (1.0, g, g)],
             'blue': [(0.0, b, b),
                      (1.0, b, b)],
             # 'alpha': [(0.0, 0.0, 0.0),
             #           (0.1, 0.2, 0.2),
             #           (0.25, 0.4, 0.4),
             #           (0.5, 0.55, 0.55),
             #           (0.8, 0.9, 0.9),
             #           (1.0, 1.0, 1.0)]}
             # original:
             'alpha': [(0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0)]}

    # Create and return the color map
    return mcolors.LinearSegmentedColormap('custom', cdict)

# -----------------------------------------------------------------

EMPTY_IMAGE = np.array([[np.nan]])
IDENTITY = IdentityTransform()

# -----------------------------------------------------------------

SUPPORTS_RESIZE = []

try: from matplotlib.backends.backend_tkagg import FigureCanvasTk
except ImportError: pass
else: SUPPORTS_RESIZE.append(FigureCanvasTk)

try: from matplotlib.backends.backend_qt5 import FigureCanvasQT
except ImportError: pass
else: SUPPORTS_RESIZE.append(FigureCanvasQT)
SUPPORTS_RESIZE = tuple(SUPPORTS_RESIZE)

# -----------------------------------------------------------------

class GenericDensityArtist(AxesImage):

    """
    Matplotlib artist to make a density plot given a helper histogram function.

    This is a more generic form of ``ScatterDensityArtist``. Here, we can
    initialize the class with a histogram function that just takes bins and the
    range of values, and returns a density array. This is useful for cases where
    the data might be changing dynamically over time.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        The axes to plot the artist into.
    dpi : int or `None`
        The number of dots per inch to include in the density map. To use
        the native resolution of the drawing device, set this to None.
    cmap : `matplotlib.colors.Colormap`
        The colormap to use for the density map.
    color : str or tuple
        The color to use for the density map. This can be any valid
        Matplotlib color. If specified, this takes precedence over the
        colormap.
    alpha : float
        Overall transparency of the density map.
    norm : `matplotlib.colors.Normalize`
        The normalization class for the density map.
    vmin, vmax : float or func
        The lower and upper levels used for scaling the density map. These can
        optionally be functions that take the density array and returns a single
        value (e.g. a function that returns the 5% percentile, or the minimum).
        This is useful since when zooming in/out, the optimal limits change.
    histogram2d_func : callable, optional
        The function (or callable instance) to use for computing the 2D
        histogram - this should take the arguments ``bins`` and ``range`` as
        defined by :func:`~numpy.histogram2d` as well as a ``pressed`` keyword
        argument that indicates whether the user is currently panning/zooming.
    kwargs
        Any additional keyword arguments are passed to AxesImage.
    """

    def __init__(self, ax, dpi=72, color=None, vmin=None, vmax=None, norm=None, histogram2d_func=None,
                 update_while_panning=True, c_with_alpha=False, c_alpha_log_tresh=0.05, c_alpha_log_enhance=50.,
                 c_alpha_logdensity=False, **kwargs):

        super(GenericDensityArtist, self).__init__(ax, **kwargs)

        self._histogram2d_func = histogram2d_func

        self._make_image_called = False
        self._pressed = False
        self._density_vmin = np.nanmin
        self._density_vmax = np.nanmax

        self._ax = ax
        self._ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self._ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

        self._update_while_panning = update_while_panning

        self.set_dpi(dpi)

        self.on_release()
        self.set_array(EMPTY_IMAGE)

        # Set special flags and options
        self.c_with_alpha = c_with_alpha
        self.c_alpha_log_tresh = c_alpha_log_tresh
        self.c_alpha_log_enhance = c_alpha_log_enhance
        self.c_alpha_logdensity = c_alpha_logdensity

        if color is not None: self.set_color(color)
        if norm is not None: self.set_norm(norm)
        if vmin is not None or vmax is not None: self.set_clim(vmin, vmax)

        # Not all backends support timers properly, so we explicitly whitelist
        # backends for which they do. In these cases, we avoid recomputing the
        # density map during resizing.
        if isinstance(self._ax.figure.canvas, SUPPORTS_RESIZE):
            self._ax.figure.canvas.mpl_connect('resize_event', self._resize_start)
            self._timer = self._ax.figure.canvas.new_timer(interval=500)
            self._timer.single_shot = True
            self._timer.add_callback(self._resize_end)
        else: self._timer = None

        # Density field
        self._density = None

    # -----------------------------------------------------------------

    def _resize_start(self, event=None):
        if not self._make_image_called:
            # Only handle resizing once the map has been shown at least once
            # to avoid 'blinking' at the start.
            return
        self.on_press(force=True)
        self._timer.start()

    # -----------------------------------------------------------------

    def _resize_end(self, event=None):
        self.on_release()
        self.stale = True
        self._ax.figure.canvas.draw()

    # -----------------------------------------------------------------

    def set_color(self, color):
        if color is not None:
            self.set_cmap(make_cmap(color))

    # -----------------------------------------------------------------

    def set_dpi(self, dpi):
        self._dpi = dpi

    # -----------------------------------------------------------------

    def on_press(self, event=None, force=False):
        if not force:
            try:
                mode = self._ax.figure.canvas.toolbar.mode
            except AttributeError:  # pragma: nocover
                return
            if mode != 'pan/zoom':
                return
        self._pressed = True
        self.stale = True

    # -----------------------------------------------------------------

    def on_release(self, event=None):
        self._pressed = False
        self.stale = True

    # -----------------------------------------------------------------

    def get_extent(self):

        if not self._update_while_panning and self._pressed:
            return self._extent

        xmin, xmax = self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()

        self._extent = xmin, xmax, ymin, ymax

        return self._extent

    # -----------------------------------------------------------------

    def get_transform(self):

        # If we don't override this, the transform includes LogTransforms
        # and the final image gets warped to be 'correct' in data space
        # since Matplotlib 2.x:
        #
        #   https://matplotlib.org/users/prev_whats_new/whats_new_2.0.0.html#non-linear-scales-on-image-plots
        #
        # However, we want pixels to always visually be the same size, so we
        # override the transform to not include the LogTransform components.

        xmin, xmax = self._ax.get_xlim()
        ymin, ymax = self._ax.get_ylim()

        bbox = BboxTransformFrom(TransformedBbox(Bbox([[xmin, ymin], [xmax, ymax]]),
                                                 IDENTITY))

        return bbox + self._ax.transAxes

    # -----------------------------------------------------------------

    def make_image(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        if not self._update_while_panning and self._pressed:
            return super(GenericDensityArtist, self).make_image(*args, **kwargs)

        xmin, xmax = self._ax.get_xlim()
        ymin, ymax = self._ax.get_ylim()

        if self._dpi is None:
            dpi = self.axes.figure.get_dpi()
        else:
            dpi = self._dpi

        width = (self._ax.get_position().width *
                 self._ax.figure.get_figwidth())
        height = (self._ax.get_position().height *
                  self._ax.figure.get_figheight())

        nx = int(round(width * dpi))
        ny = int(round(height * dpi))

        flip_x = xmin > xmax
        flip_y = ymin > ymax

        if flip_x:
            xmin, xmax = xmax, xmin

        if flip_y:
            ymin, ymax = ymax, ymin

        bins = (ny, nx)

        array, density = self._histogram2d_func(bins=bins, range=((ymin, ymax), (xmin, xmax)), return_density=True)

        # No auxilary axis
        if self._histogram2d_func._c is None and self.c_alpha_logdensity:
            #print(np.nanmin(array), np.nanmax(array))
            array = alpha_from_densities(array, tresfrac=self.c_alpha_log_tresh, scaleval=self.c_alpha_log_enhance, logscale=True)
            #print(np.nanmin(array), np.nanmax(array))
            #print(array)
            array -= np.nanmin(array)

        if flip_x or flip_y:
            if flip_x and flip_y:
                array = array[::-1, ::-1]
                if density is not None: density = density[::-1, ::-1]
            elif flip_x:
                array = array[:, ::-1]
                if density is not None: density = density[:, ::-1]
            else:
                array = array[::-1, :]
                if density is not None: density = density[::-1, :]

        if self.origin == 'upper':
            array = np.flipud(array)
            if density is not None: density = np.flipud(density)

        if callable(self._density_vmin):
            vmin = self._density_vmin(array)
        else:
            vmin = self._density_vmin

        if callable(self._density_vmax):
            vmax = self._density_vmax(array)
        else:
            vmax = self._density_vmax

        # Set data determining the color (linked to colormap)
        self.set_data(array)

        # Set the density for when auxilary axis is used
        self._density = density

        super(GenericDensityArtist, self).set_clim(vmin, vmax)

        self._make_image_called = True

        if self.c_with_alpha: return self.make_image_special(*args, **kwargs)
        else: return super(GenericDensityArtist, self).make_image(*args, **kwargs)

    # -----------------------------------------------------------------

    def make_image_special(self, renderer, magnification=1.0, unsampled=False):
        ## This function is a copy of AxesImage.make_image(*args, **kwargs), but with _make_image replaced by _make_image_special

        trans = self.get_transform()
        # image is created in the canvas coordinate.
        x1, x2, y1, y2 = self.get_extent()
        bbox = Bbox(np.array([[x1, y1], [x2, y2]]))
        transformed_bbox = TransformedBbox(bbox, trans)

        return self._make_image_special(
            self._A, self._density, bbox, transformed_bbox, self.axes.bbox, magnification,
            unsampled=unsampled)

    # -----------------------------------------------------------------

    def _make_image_special(self, A, density, in_bbox, out_bbox, clip_bbox, magnification=1.0,
                    unsampled=False, round_to_pixel_border=True, densities=None):
        """
        This function is a copy of _ImageBase._make_image(*args, **kwargs), but with only the A.ndim == 2 case,
        and with the transformation of the alpha channel on top of the normal behaviour with the colormap's RGBA
        """

        if A is None:
            raise RuntimeError('You must first set the image '
                               'array or the image attribute')
        if A.size == 0:
            raise RuntimeError("_make_image must get a non-empty image. "
                               "Your Artist's draw method must filter before "
                               "this method is called.")

        clipped_bbox = Bbox.intersection(out_bbox, clip_bbox)

        if clipped_bbox is None:
            return None, 0, 0, None

        out_width_base = clipped_bbox.width * magnification
        out_height_base = clipped_bbox.height * magnification

        if out_width_base == 0 or out_height_base == 0:
            return None, 0, 0, None

        if self.origin == 'upper':
            # Flip the input image using a transform.  This avoids the
            # problem with flipping the array, which results in a copy
            # when it is converted to contiguous in the C wrapper
            t0 = Affine2D().translate(0, -A.shape[0]).scale(1, -1)
        else:
            t0 = IdentityTransform()

        t0 += (
                Affine2D()
                .scale(
                    in_bbox.width / A.shape[1],
                    in_bbox.height / A.shape[0])
                .translate(in_bbox.x0, in_bbox.y0)
                + self.get_transform())

        t = (t0
             + Affine2D().translate(
                    -clipped_bbox.x0,
                    -clipped_bbox.y0)
             .scale(magnification, magnification))

        # So that the image is aligned with the edge of the axes, we want
        # to round up the output width to the next integer.  This also
        # means scaling the transform just slightly to account for the
        # extra subpixel.
        if (t.is_affine and round_to_pixel_border and
                (out_width_base % 1.0 != 0.0 or out_height_base % 1.0 != 0.0)):
            out_width = int(ceil(out_width_base))
            out_height = int(ceil(out_height_base))
            extra_width = (out_width - out_width_base) / out_width_base
            extra_height = (out_height - out_height_base) / out_height_base
            t += Affine2D().scale(1.0 + extra_width, 1.0 + extra_height)
        else:
            out_width = int(out_width_base)
            out_height = int(out_height_base)

        if not unsampled:
            #if A.ndim not in (2, 3):
            #    raise ValueError("Invalid dimensions, got {}".format(A.shape))

            # if we are a 2D array, then we are running through the
            # norm + colormap transformation.  However, in general the
            # input data is not going to match the size on the screen so we
            # have to resample to the correct number of pixels
            # need to

            # TODO slice input array first
            inp_dtype = A.dtype
            a_min = A.min()
            a_max = A.max()
            # figure out the type we should scale to.  For floats,
            # leave as is.  For integers cast to an appropriate-sized
            # float.  Small integers get smaller floats in an attempt
            # to keep the memory footprint reasonable.
            if a_min is np.ma.masked:
                # all masked, so values don't matter
                a_min, a_max = np.int32(0), np.int32(1)
            if inp_dtype.kind == 'f':
                scaled_dtype = A.dtype
            else:
                # probably an integer of some type.
                da = a_max.astype(np.float64) - a_min.astype(np.float64)
                if da > 1e8:
                    # give more breathing room if a big dynamic range
                    scaled_dtype = np.float64
                else:
                    scaled_dtype = np.float32

            # scale the input data to [.1, .9].  The Agg
            # interpolators clip to [0, 1] internally, use a
            # smaller input scale to identify which of the
            # interpolated points need to be should be flagged as
            # over / under.
            # This may introduce numeric instabilities in very broadly
            # scaled data
            A_scaled = np.empty(A.shape, dtype=scaled_dtype)
            A_scaled[:] = A
            # clip scaled data around norm if necessary.
            # This is necessary for big numbers at the edge of
            # float64's ability to represent changes.  Applying
            # a norm first would be good, but ruins the interpolation
            # of over numbers.
            self.norm.autoscale_None(A)
            dv = (np.float64(self.norm.vmax) -
                  np.float64(self.norm.vmin))
            vmid = self.norm.vmin + dv / 2
            fact = 1e7 if scaled_dtype == np.float64 else 1e4
            newmin = vmid - dv * fact
            if newmin < a_min:
                newmin = None
            else:
                a_min = np.float64(newmin)
            newmax = vmid + dv * fact
            if newmax > a_max:
                newmax = None
            else:
                a_max = np.float64(newmax)
            if newmax is not None or newmin is not None:
                A_scaled = np.clip(A_scaled, newmin, newmax)

            A_scaled -= a_min
            # a_min and a_max might be ndarray subclasses so use
            # asscalar to avoid errors
            a_min = np.asscalar(a_min.astype(scaled_dtype))
            a_max = np.asscalar(a_max.astype(scaled_dtype))

            if a_min != a_max:
                A_scaled /= ((a_max - a_min) / 0.8)
            A_scaled += 0.1
            A_resampled = np.zeros((out_height, out_width),
                                   dtype=A_scaled.dtype)
            # resample the input data to the correct resolution and shape
            _image.resample(A_scaled, A_resampled,
                            t,
                            _interpd_[self.get_interpolation()],
                            self.get_resample(), 1.0,
                            self.get_filternorm() or 0.0,
                            self.get_filterrad() or 0.0)

            #alpha = self.get_alpha()

            new_density = np.zeros((out_height, out_width),
                                 dtype=density.dtype)
            _image.resample(density, new_density,
                            t,
                            _interpd_[self.get_interpolation()],
                            self.get_resample(), 1.0,
                            self.get_filternorm() or 0.0,
                            self.get_filterrad() or 0.0)

            # we are done with A_scaled now, remove from namespace
            # to be sure!
            del A_scaled
            # un-scale the resampled data to approximately the
            # original range things that interpolated to above /
            # below the original min/max will still be above /
            # below, but possibly clipped in the case of higher order
            # interpolation + drastically changing data.
            A_resampled -= 0.1
            if a_min != a_max:
                A_resampled *= ((a_max - a_min) / 0.8)
            A_resampled += a_min
            # if using NoNorm, cast back to the original datatype
            if isinstance(self.norm, mcolors.NoNorm):
                A_resampled = A_resampled.astype(A.dtype)

            mask = np.empty(A.shape, dtype=np.float32)
            if A.mask.shape == A.shape:
                # this is the case of a nontrivial mask
                mask[:] = np.where(A.mask, np.float32(np.nan),
                                   np.float32(1))
            else:
                mask[:] = 1

            # we always have to interpolate the mask to account for
            # non-affine transformations
            out_mask = np.zeros((out_height, out_width),
                                dtype=mask.dtype)
            _image.resample(mask, out_mask,
                            t,
                            _interpd_[self.get_interpolation()],
                            True, 1,
                            self.get_filternorm() or 0.0,
                            self.get_filterrad() or 0.0)
            # we are done with the mask, delete from namespace to be sure!
            del mask
            # Agg updates the out_mask in place.  If the pixel has
            # no image data it will not be updated (and still be 0
            # as we initialized it), if input data that would go
            # into that output pixel than it will be `nan`, if all
            # the input data for a pixel is good it will be 1, and
            # if there is _some_ good data in that output pixel it
            # will be between [0, 1] (such as a rotated image).

            out_alpha = np.array(out_mask)
            out_mask = np.isnan(out_mask)
            out_alpha[out_mask] = 1
            new_alpha = alpha_from_densities(new_density, tresfrac=self.c_alpha_log_tresh, scaleval=self.c_alpha_log_enhance, logscale=self.c_alpha_logdensity)
            out_alpha *= new_alpha

            # mask and run through the norm
            output = self.norm(np.ma.masked_array(A_resampled, out_mask))

            # at this point output is either a 2D array of normed data
            # (of int or float)
            # or an RGBA array of re-sampled input
            output = self.to_rgba(output, bytes=True, norm=False)
            # output is now a correctly sized RGBA array of uint8

            # Apply alpha *after* if the input was greyscale without a mask
            if A.ndim == 2:
                # alpha = self.get_alpha()
                # if alpha is None:
                #    alpha = 1
                alpha = 1
                alpha_channel = output[:, :, 3]
                alpha_channel[:] = np.asarray(
                    np.asarray(alpha_channel, np.float32) * out_alpha * alpha,
                    np.uint8)

        else:
            if self._imcache is None:
                self._imcache = self.to_rgba(A, bytes=True, norm=(A.ndim == 2))
            output = self._imcache

            # Subset the input image to only the part that will be
            # displayed
            subset = TransformedBbox(
                clip_bbox, t0.frozen().inverted()).frozen()
            output = output[
                     int(max(subset.ymin, 0)):
                     int(min(subset.ymax + 1, output.shape[0])),
                     int(max(subset.xmin, 0)):
                     int(min(subset.xmax + 1, output.shape[1]))]

            t = Affine2D().translate(
                int(max(subset.xmin, 0)), int(max(subset.ymin, 0))) + t

        return output, clipped_bbox.x0, clipped_bbox.y0, t

    # -----------------------------------------------------------------

    def set_clim(self, vmin, vmax):
        self._density_vmin = vmin
        self._density_vmax = vmax

    # -----------------------------------------------------------------

    def set_norm(self, norm):
        if norm.vmin is not None:
            self._density_vmin = norm.vmin
        if norm.vmax is not None:
            self._density_vmax = norm.vmax
        super(GenericDensityArtist, self).set_norm(norm)

    # -----------------------------------------------------------------

    def remove(self):
        if self._timer is not None:
            self._timer.stop()
            self._timer = None
        super(GenericDensityArtist, self).remove()
        # We explicitly clean up the reference to the histogram2d function since
        # this may in some cases cause circular references.
        self._histogram2d_func = None

# -----------------------------------------------------------------

def log_scale(x, a=100):
    return np.log10(1+(a*x)) / np.log10(1+a)

# -----------------------------------------------------------------

def simple_log(x):
    return np.log10(x) / np.log10(x.max())

# -----------------------------------------------------------------

def alpha_from_densities(densities, tresfrac=0.05, scaleval=10., logscale=True):

    """
    This function ...
    :param densities:
    :param tresfrac:
    :param scaleval:
    :param logscale:
    :return:
    """

    densities = densities + (np.max(densities) * tresfrac)  # lower bound for alpha

    if logscale:
        if scaleval > 0:
            counts = densities / np.max(densities)
            alpha = log_scale(counts, scaleval)
        else: alpha = simple_log(densities)
    else: alpha = densities / np.max(densities)

    # Return alpha
    return alpha

# -----------------------------------------------------------------

class FixedDataDensityHelper:

    compute_when_pressed = True

    def __init__(self, ax, x, y, c=None, downres_factor=4):

        self._ax = ax
        self._c = None
        self._downres = False

        if downres_factor < 1 or downres_factor % 1 != 0:
            raise ValueError('downres_factor should be a strictly positive integer value')

        self._downres_factor = downres_factor
        self.set_xy(x, y)
        self.set_c(c)

    def downres(self):
        self._downres = True

    def upres(self):
        self._downres = False

    def set_xy(self, x, y):
        self._x = x
        self._y = y
        self._x_log = None
        self._y_log = None
        self._x_log_sub = None
        self._y_log_sub = None
        step = self._downres_factor ** 2
        self._x_sub = self._x[::step]
        self._y_sub = self._y[::step]

    def set_c(self, c):
        self._c = c
        step = self._downres_factor ** 2
        if self._c is None:
            self._c_sub = None
        else:
            self._c_sub = self._c[::step]

    def _update_x_log(self):
        step = self._downres_factor ** 2
        with np.errstate(invalid='ignore'):
            self._x_log = np.log10(self._x)
        self._x_log_sub = self._x_log[::step]

    def _update_y_log(self):
        step = self._downres_factor ** 2
        with np.errstate(invalid='ignore'):
            self._y_log = np.log10(self._y)
        self._y_log_sub = self._y_log[::step]

    def __call__(self, bins=None, range=None, return_density=False):

        ny, nx = bins
        (ymin, ymax), (xmin, xmax) = range

        xscale = self._ax.get_xscale()
        yscale = self._ax.get_yscale()

        if xscale == 'log':
            xmin, xmax = log10(xmin), log10(xmax)
            if self._x_log is None:
                # We do this here insead of in set_xy to save time since in
                # set_xy we don't know yet if the axes will be log or not.
                self._update_x_log()
            if self._downres:
                x = self._x_log_sub
            else:
                x = self._x_log
        elif xscale == 'linear':
            if self._downres:
                x = self._x_sub
            else:
                x = self._x
        else:  # pragma: nocover
            raise ValueError('Unexpected xscale: {0}'.format(xscale))

        if yscale == 'log':
            ymin, ymax = log10(ymin), log10(ymax)
            if self._y_log is None:
                # We do this here insead of in set_xy to save time since in
                # set_xy we don't know yet if the axes will be log or not.
                self._update_y_log()
            if self._downres:
                y = self._y_log_sub
            else:
                y = self._y_log
        elif yscale == 'linear':
            if self._downres:
                y = self._y_sub
            else:
                y = self._y
        else:  # pragma: nocover
            raise ValueError('Unexpected xscale: {0}'.format(xscale))

        if self._downres:
            nx_sub = nx // self._downres_factor
            ny_sub = ny // self._downres_factor
            bins = (ny_sub, nx_sub)
            weights = self._c_sub
        else:
            bins = (ny, nx)
            weights = self._c

        if weights is None:
            array = histogram2d(y, x, bins=bins,
                                range=((ymin, ymax), (xmin, xmax)))
            density = array / np.nanmax(array)
        else:
            array = histogram2d(y, x, bins=bins, weights=weights,
                                range=((ymin, ymax), (xmin, xmax)))
            count = histogram2d(y, x, bins=bins,
                                range=((ymin, ymax), (xmin, xmax)))
            density = count / np.nanmax(count)

            with np.errstate(invalid='ignore'):
                array /= count

        if return_density: return array, density
        else: return array

# -----------------------------------------------------------------

class ScatterDensityArtist(GenericDensityArtist):

    """
    Matplotlib artist to make a density plot of (x, y) scatter data.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        The axes to plot the artist into.
    x, y : iterable
        The data to plot.
    c : iterable
        Values to use for color-encoding. This is meant to be the same as
        the argument with the same name in :meth:`~matplotlib.axes.Axes.scatter`
        although for now only 1D iterables of values are accepted. Note that
        values are averaged inside each pixel of the density map *before*
        applying the colormap, which in some cases will be different from what
        the average color of markers would have been inside each pixel.
    dpi : int or `None`
        The number of dots per inch to include in the density map. To use
        the native resolution of the drawing device, set this to None.
    downres_factor : int
        For interactive devices, when panning, the density map will
        automatically be made at a lower resolution and including only a
        subset of the points. The new dpi of the figure when panning will
        then be dpi / downres_factor, and the number of elements in the
        arrays will be reduced by downres_factor**2.
    cmap : `matplotlib.colors.Colormap`
        The colormap to use for the density map.
    color : str or tuple
        The color to use for the density map. This can be any valid
        Matplotlib color. If specified, this takes precedence over the
        colormap.
    alpha : float
        Overall transparency of the density map.
    norm : `matplotlib.colors.Normalize`
        The normalization class for the density map.
    vmin, vmax : float or func
        The lower and upper levels used for scaling the density map. These can
        optionally be functions that take the density array and returns a single
        value (e.g. a function that returns the 5% percentile, or the minimum).
        This is useful since when zooming in/out, the optimal limits change.
    update_while_panning : bool, optional
        Whether to compute histograms on-the-fly while panning.
    kwargs
        Any additional keyword arguments are passed to AxesImage.
    """

    def __init__(self, ax, x, y, downres_factor=4, c=None, **kwargs):
        self.histogram2d_helper = FixedDataDensityHelper(ax, x, y, c=c, downres_factor=downres_factor)
        super(ScatterDensityArtist, self).__init__(ax, histogram2d_func=self.histogram2d_helper, **kwargs)

    # -----------------------------------------------------------------

    def set_xy(self, x, y):
        self.histogram2d_helper.set_xy(x, y)

    # -----------------------------------------------------------------

    def set_c(self, c):
        self.histogram2d_helper.set_c(c)

    # -----------------------------------------------------------------

    def on_press(self, event=None, force=False):
        if not force:
            if self._update_while_panning and self.histogram2d_helper._downres_factor == 1:
                return
        self.histogram2d_helper.downres()
        return super(ScatterDensityArtist, self).on_press(force=force)

    # -----------------------------------------------------------------

    def on_release(self, event=None):
        self.histogram2d_helper.upres()
        return super(ScatterDensityArtist, self).on_release()

# -----------------------------------------------------------------
