#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from image.inpaint import replace_nans

import numpy as np

def in_paint(data, mask):

    # Fill the data with nans according to the mask
    data_ma = np.ma.array(data.astype(float), mask=mask)
    data_nans = data_ma.filled(np.NaN)

    interpolated = replace_nans(data_nans, 5, 0.5, 2, "localmean")

    return interpolated