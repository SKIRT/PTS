#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.plot.error Contains the ErrorPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
from scipy import ndimage
import copy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.gridspec as gridspec
import glob
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import LogNorm
import pyfits
from collections import OrderedDict
from textwrap import wrap

from astropy.io import fits
from pyfits import PrimaryHDU, Header
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import aplpy
import wcsaxes
import matplotlib.colors as mpl_colors
import matplotlib.colorbar as mpl_colorbar

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...magic.tools.plotting import colours

# -----------------------------------------------------------------

class ErrorPlotter(object):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        :param title:
        """

        # Set the title
        self.title = title

        # Figure and grid
        self._figure = None
        self._grid = None

        # Properties
        self.style = "dark" # "dark" or "light"
        self.transparent = True
        self.format = None
        self.colormap = "viridis"
        self.vmin = None

    # -----------------------------------------------------------------

    def add_image(self, image):

        """
        This function ...
        :return:
        """



    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function ...
        :param path:
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        fig_a.hist2d(np.log10(old + yng + new), np.log10(tot), cmin=1, norm=LogNorm(), bins=200, cmap=plt.cm.BuGn_r)

# -----------------------------------------------------------------
