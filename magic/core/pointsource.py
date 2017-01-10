#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.pointsource Contains the PointSource class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .source import Source

# -----------------------------------------------------------------

class PointSource(Source):

    """
    This class...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSource, self).__init__(**kwargs)

        # Set other properties
        self.catalog = kwargs.pop("catalog", None)
        self.id = kwargs.pop("id", None)
        self.ra_error = kwargs.pop("ra_error", None)
        self.dec_error = kwargs.pop("dec_error", None)
        self.confidence = kwargs.pop("confidence", None)
        self.magnitudes = kwargs.pop("magnitudes", dict())
        self.magnitude_errors = kwargs.pop("magnitude_errors", dict())

        # PSF model
        self.psf_model = None

        # FWHM
        self.fwhm = None

        # Saturation detection
        self.saturation = None

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.psf_model is not None

    # -----------------------------------------------------------------

    @property
    def has_saturation(self):

        """
        This function ...
        :return:
        """

        return self.saturation is not None

    # -----------------------------------------------------------------

    def find_contour(self, frame, config, saturation=False):

        """
        This function ...
        :param frame:
        :param config:
        :param saturation:
        :return:
        """

        # Determine which box and mask
        if saturation:

            box = self.saturation.cutout
            mask = self.saturation.mask

            # Implementation
            self._find_contour_impl(frame.wcs, box, mask, config)

        # Call the base class implementation
        else: super(PointSource, self).find_contour(frame, config)

# -----------------------------------------------------------------
