#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.wavelengthgrids Contains the WavelengthGridGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import bisect
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..core.emissionlines import EmissionLines

# -----------------------------------------------------------------

# The names of the subgrids
subgrids = ["UV", "optical", "PAH", "dust", "extension"]

# Define the limits of the subgrids
limits = dict()
limits["UV"] = (0.02, 0.085)
limits["optical"] = (0.085, 3.)
limits["PAH"] = (3., 27.)
limits["dust"] = (27., 1000.)
limits["extension"] = (1000., 2000)

# Define the relative fineness (the number of points) of the subgrids
relpoints = dict()
relpoints["UV"] = 25./325.    # 25
relpoints["optical"] = 100./325.     # 100
relpoints["PAH"] = 125./325.  # 125
relpoints["dust"] = 50./325.  # 50
relpoints["extension"] = 25./325.  # 25

# -----------------------------------------------------------------

class WavelengthGridGenerator(FittingComponent):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridGenerator, self).__init__()

        # -- Attributes --

        self.grids = []

        self.emission_lines = None

    # -----------------------------------------------------------------

    def run(self, min_npoints, max_npoints, ngrids, fixed=None): # fixed wavelength points: I1 AND FUV !!

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Generate the grids
        self.generate(min_npoints, max_npoints, ngrids, fixed)

        # 17. Writing
        #self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Emission lines
        self.emission_lines = EmissionLines()

    # -----------------------------------------------------------------

    def generate(self, min_npoints, max_npoints, ngrids, fixed=None):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the wavelength grids ...")

        # Determine the increment
        increment = (max_npoints - min_npoints) / float(ngrids)

        # List of floats
        values = np.arange(min_npoints, max_npoints, increment)

        npoints_list = []
        for value in values:
            if int(value) == npoints_list[-1]: continue
            else: npoints_list.append(int(value))

        # Loop over the different number of points
        for npoints in npoints_list:

            # Create the grid and add it to the list
            grid = self.create_grid(npoints, fixed=fixed)
            self.grids.append(grid)

    # -----------------------------------------------------------------

    def create_grid(self, npoints, add_emission_lines=False, fixed=None):

        """
        This function ...
        :param npoints:
        :param add_emission_lines:
        :param fixed:
        :return:
        """

        # Inform the user
        with_without = " with " if add_emission_lines else " without "
        log.info("Creating a wavelength grid with " + str(npoints) + with_without + "emission lines ...")

        # A list of the wavelength points
        wavelengths = []

        # Loop over the subgrids
        for subgrid in subgrids:

            # Debugging
            log.debug("Adding the " + subgrid + " subgrid ...")

            # Determine minimum, maximum and number of wavelength points for this subgrid
            min_lambda = limits[subgrid][0]
            max_lambda = limits[subgrid][1]
            points = int(round(relpoints[subgrid] * npoints))

            # Generate and add the wavelength points
            wavelengths += make_grid(min_lambda, max_lambda, points)

        # Add the emission lines
        if add_emission_lines:

            # Add emission line grid points
            logdelta = 0.001
            for line in self.emission_lines:

                center = line.center
                left = line.left
                right = line.right

                # logcenter = np.log10(center)
                logleft = np.log10(left if left > 0 else center) - logdelta
                logright = np.log10(right if right > 0 else center) + logdelta
                newgrid = []
                for w in wavelengths:
                    logw = np.log10(w)
                    if logw < logleft or logw > logright:
                        newgrid.append(w)
                newgrid.append(center)
                if left > 0:
                    newgrid.append(left)
                if right > 0:
                    newgrid.append(right)
                wavelengths = newgrid

        # Add fixed wavelength points
        if fixed is not None:
            for wavelength in fixed: wavelengths.append(wavelength)

        # Sort the wavelength points
        wavelengths = sorted(wavelengths)

        # Create the wavelength grid
        grid = WavelengthGrid.from_wavelengths(wavelengths)

        # Return the grid
        return grid

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------

def make_grid(wmin, wmax, N):

    """
    This function returns a wavelength grid (in micron) with a given resolution (nr of points per decade)
    # in the specified range (in micron), aligned with the 10^n grid points.
    """

    result = []

    # generate wavelength points p on a logarithmic scale with lambda = 10**p micron
    #  -2 <==> 0.01
    #   4 <==> 10000
    for i in range(-2*N,4*N+1):
        p = float(i)/N
        w = 10.**p
        if wmin <= w < wmax: result.append(w)

    # Return the grid
    return result

# -----------------------------------------------------------------
