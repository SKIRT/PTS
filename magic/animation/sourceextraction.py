#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.animation.sourceextraction Contains the SourceExtractionAnimation class.

# -----------------------------------------------------------------

# Import standard modules
import io
import imageio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as plt_Rectangle
from matplotlib.patches import Ellipse as plt_Ellipse

# Import astronomical modules
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# Import the relevant PTS classes and modules
from ...core.basics.animation import Animation
from ..region.ellipse import PixelEllipseRegion
from ..basics.mask import Mask

# -----------------------------------------------------------------

class SourceExtractionAnimation(Animation):

    """
    This class ...
    """

    def __init__(self, frame):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SourceExtractionAnimation, self).__init__()

        # Make a reference to the image frame
        self.frame = frame

        # Set the number of frames per second
        self.fps = 1 # 1 frame per second because the frames are very distinct

        # Have to be set
        self.principal_shape = None
        self.mask = Mask.empty_like(self.frame)

    # -----------------------------------------------------------------

    def add_source(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        # Create buffer and figure
        buf = io.BytesIO()
        fig = plt.figure(figsize=(15, 15))

        # Add first subplot
        ax = fig.add_subplot(2, 3, 1)

        # Plot the frame (before removal)
        norm = ImageNormalize(stretch=LogStretch())
        norm_residual = ImageNormalize(stretch=LogStretch())
        min_value = np.nanmin(self.frame)
        max_value = 0.5 * (np.nanmax(self.frame) + min_value)

        ax.imshow(self.frame, origin="lower", interpolation="nearest", vmin=min_value, vmax=max_value, norm=norm)

        r = plt_Rectangle((source.x_min, source.y_min), source.xsize, source.ysize, edgecolor="yellow",
                          facecolor="none", lw=5)
        ax.add_patch(r)

        if isinstance(self.principal_shape, PixelEllipseRegion):

            ell = plt_Ellipse((self.principal_shape.center.x, self.principal_shape.center.y),
                              2.0 * self.principal_shape.radius.x, 2.0 * self.principal_shape.radius.y,
                              self.principal_shape.angle.to("deg").value, edgecolor="red", facecolor="none", lw=5)
            ax.add_patch(ell)

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Add the second subplot
        ax = fig.add_subplot(2, 3, 2)

        frame_xsize = self.frame.xsize
        frame_ysize = self.frame.ysize

        new_xmin = int(source.center.x - 0.1 * frame_xsize)
        new_xmax = int(source.center.x + 0.1 * frame_xsize)
        new_ymin = int(source.center.y - 0.1 * frame_ysize)
        new_ymax = int(source.center.y + 0.1 * frame_ysize)

        if new_xmin < 0: new_xmin = 0
        if new_ymin < 0: new_ymin = 0
        if new_xmax > frame_xsize: new_xmax = frame_xsize
        if new_ymax > frame_ysize: new_ymax = frame_ysize

        ax.imshow(self.frame[new_ymin:new_ymax, new_xmin:new_xmax], origin="lower", interpolation="nearest",
                  vmin=min_value, vmax=max_value, norm=norm)
        r = plt_Rectangle((source.x_min - new_xmin, source.y_min - new_ymin), source.xsize, source.ysize,
                          edgecolor="yellow", facecolor="none", lw=5)
        ax.add_patch(r)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Add the third subplot
        ax = fig.add_subplot(2, 3, 3)

        # Adapt the mask
        self.mask[source.y_slice, source.x_slice] += source.mask

        # Plot the mask
        # ax.imshow(self.mask[new_ymin:new_ymax, new_xmin:new_xmax], origin="lower", cmap='Greys')
        ax.imshow(self.mask, origin="lower", cmap="Greys")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Add the fourth subplot
        ax = fig.add_subplot(2, 3, 4)

        # Plot the source cutout
        # ax.imshow(source.cutout, origin="lower", interpolation="nearest", vmin=min_value, vmax=max_value, norm=norm)
        ax.imshow(np.ma.MaskedArray(source.cutout, mask=source.background_mask), origin="lower",
                  interpolation="nearest", vmin=min_value, vmax=max_value, norm=norm)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # ell = plt_Ellipse((source.center.x-source.x_min, source.center.y-source.y_min), 2.0*source.radius.x, 2.0*source.radius.y, edgecolor="yellow", facecolor="none", lw=5)
        # ax.add_patch(ell)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Add the fourth subplot
        ax = fig.add_subplot(2, 3, 5)

        # Plot the source background
        # ax.imshow(np.ma.MaskedArray(source.cutout, mask=source.mask), origin="lower", interpolation="nearest", vmin=min_value, vmax=max_value, norm=norm)
        new_cutout = source.cutout.copy()
        new_cutout[source.mask] = source.background[source.mask]
        ax.imshow(new_cutout, origin="lower", interpolation="nearest", vmin=min_value, vmax=max_value, norm=norm)
        # ell = plt_Ellipse((source.center.x-source.x_min, source.center.y-source.y_min), 2.0 * source.radius.x, 2.0 * source.radius.y, edgecolor="yellow", facecolor="none", lw=5)
        # ax.add_patch(ell)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Add the sixth subplot
        ax = fig.add_subplot(2, 3, 6)

        # Plot the residual cutout (background removed)
        ax.imshow(source.subtracted, origin="lower", interpolation="nearest", vmin=0.0, vmax=max_value,
                  norm=norm_residual)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        plt.tight_layout()

        plt.savefig(buf, format="png")
        plt.close()
        buf.seek(0)
        im = imageio.imread(buf)
        buf.close()

        self.add_frame(im)

# -----------------------------------------------------------------
