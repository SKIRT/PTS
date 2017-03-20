#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.wavelengthgrids Contains the WavelengthGridGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..basics.emissionlines import EmissionLines, EmissionLine
from ..basics.configurable import Configurable
from ..basics.table import SmartTable
from ..basics.range import QuantityRange
from ..basics.unit import parse_unit as u
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..plot.transmission import TransmissionPlotter
from ..basics.range import IntegerRange

# -----------------------------------------------------------------

# The names of the subgrids
subgrids = ["UV", "optical", "PAH", "dust", "extension"]

# Define the ranges of the subgrids
ranges = dict()
ranges["UV"] = QuantityRange(0.02, 0.085, unit="micron")
ranges["optical"] = QuantityRange(0.085, 3., unit="micron")
ranges["PAH"] = QuantityRange(3., 27., unit="micron")
ranges["dust"] = QuantityRange(27., 1000., unit="micron")
ranges["extension"] = QuantityRange(1000., 2000, unit="micron")

# Define the relative fineness (the number of points) of the subgrids
relpoints = dict()
relpoints["UV"] = 25./325.    # 25
relpoints["optical"] = 100./325.     # 100
relpoints["PAH"] = 125./325.  # 125
relpoints["dust"] = 50./325.  # 50
relpoints["extension"] = 25./325.  # 25

# -----------------------------------------------------------------

class WavelengthGridsTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(WavelengthGridsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("UV points", int, None, "number of points in UV spectrum (range: " + str(ranges["UV"]) + ")")
        self.add_column_info("Optical points", int, None, "number of points in the optical spectrum (range: " + str(ranges["optical"]) + ")")
        self.add_column_info("PAH points", int, None, "number of points in the PAH spectrum (range: " + str(ranges["PAH"]) + ")")
        self.add_column_info("Dust points", int, None, "number of points in the dust spectrum (range: " + str(ranges["dust"]) + ")")
        self.add_column_info("Extension points", int, None, "number of points in the extension spectrum (range: " + str(ranges["extension"]) + ")")
        self.add_column_info("Broad band filters", str, None, "broad band filters for which the wavelenth range was resampled")
        self.add_column_info("Narrow band filters", str, None, "narrow band filters for which the wavelength was added")
        self.add_column_info("Emission lines", int, None, "number of emission lines")
        self.add_column_info("Fixed points", int, None, "number of fixed points")
        self.add_column_info("Total points", int, None, "total number of points")

    # -----------------------------------------------------------------

    def add_grid(self, grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added):

        """
        This function ...
        :param grid:
        :param subgrid_npoints:
        :param emission_npoints:
        :param fixed_npoints:
        :param broad_resampled:
        :param narrow_added:
        :return:
        """

        # Get values
        uv_npoints = subgrid_npoints["UV"] if "UV" in subgrid_npoints else 0
        optical_npoints = subgrid_npoints["optical"] if "optical" in subgrid_npoints else 0
        pah_npoints = subgrid_npoints["PAH"] if "PAH" in subgrid_npoints else 0
        dust_npoints = subgrid_npoints["dust"] if "dust" in subgrid_npoints else 0
        extension_npoints = subgrid_npoints["extension"] if "extension" in subgrid_npoints else 0

        # Create strings from the filter lists
        broad_string = ",".join(broad_resampled)
        narrow_string = ",".join(narrow_added)

        # Add row
        self.add_row([uv_npoints, optical_npoints, pah_npoints, dust_npoints, extension_npoints, broad_string,
                      narrow_string, emission_npoints, fixed_npoints, len(grid)])

# -----------------------------------------------------------------

class WavelengthGridGenerator(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridGenerator, self).__init__(config, interactive)

        # -- Attributes --

        # Settings
        self.npoints_range = None
        self.ngrids = None
        self.fixed = None
        self.add_emission_lines = False
        self.lines = None
        self.min_wavelength = None
        self.max_wavelength = None
        self.filters = None

        # The wavelength grids
        self.grids = []

        # The wavelength grid property table
        self.table = None

        # The emission line object
        self.emission_lines = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Generate the grids
        self.generate()

        # 3. Show
        if self.config.show: self.show()

        # 4. Plot
        if self.config.plot: self.plot()

        # 5. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    @property
    def single_grid(self):

        """
        This function ...
        :return:
        """

        if len(self.grids) == 0: raise RuntimeError("No grid")
        elif len(self.grids) == 1: return self.grids[0]
        else: raise RuntimeError("More than one grid")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Set options
        self.ngrids = kwargs.pop("ngrids")
        if self.ngrids == 1: self.npoints_range = IntegerRange.infinitesimal(kwargs.pop("npoints"))
        else: self.npoints_range = kwargs.pop("npoints_range")
        self.fixed = kwargs.pop("fixed", None)
        self.add_emission_lines = kwargs.pop("add_emission_lines", False)
        self.lines = kwargs.pop("lines", None)
        self.min_wavelength = kwargs.pop("min_wavelength", None)
        self.max_wavelength = kwargs.pop("max_wavelength", None)
        self.filters = kwargs.pop("filters", None)

        # Create the emission lines instance
        if self.add_emission_lines:

            # Use all liness
            if self.lines is None: self.emission_lines = EmissionLines()

            # Use specific lines
            else: self.emission_lines = [EmissionLine.from_string(string) for string in self.lines]

        # Initialize the table
        self.table = WavelengthGridsTable()

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :param:
        """

        # Inform the user
        log.info("Generating the wavelength grids ...")

        # Loop over the different number of points
        for npoints in self.npoints_range.linear(self.ngrids):

            # Create the grid and add it to the list
            self.create_grid(npoints)

    # -----------------------------------------------------------------

    def create_grid(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        # Inform the user
        with_without = " with " if self.add_emission_lines else " without "
        log.info("Creating a wavelength grid with " + str(npoints) + " points" + with_without + "emission lines ...")

        # Generate the grid
        grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added = \
            create_one_subgrid_wavelength_grid(npoints, self.emission_lines, fixed=self.fixed,
                                               min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength,
                                               filters=self.filters, min_wavelengths_in_filter=self.config.min_wavelengths_in_filter,
                                               min_wavelengths_in_fwhm=self.config.min_wavelengths_in_fwhm)

        # Debugging
        log.debug("Generated a wavelength grid with:")
        log.debug("")
        log.debug(" - number of points: " + str(len(grid)))
        log.debug(" - number of points in subgrids: " + str(subgrid_npoints))
        log.debug(" - number of emission points: " + str(emission_npoints))
        log.debug(" - number of fixed points: " + str(fixed_npoints))
        log.debug(" - filters for which extra sampling was performed: " + str(broad_resampled))
        log.debug(" - narrow band filters for which wavelength was added: " + str(narrow_added))
        log.debug("")

        # Add the grid
        self.grids.append(grid)

        # Add entry to the table
        self.table.add_grid(grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the wavelength grids ...")

        # Loop over the grids
        for grid in self.grids:

            print("Wavelength grid with " + str(len(grid)) + " wavelength points:")
            print("")
            print(grid)
            print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        #plotter = TransmissionPlotter()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the grids
        self.write_grids()

        # Write the table
        self.write_table()

    # -----------------------------------------------------------------

    def write_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the grids ...")

        pass

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the grids table ...")

        pass

# -----------------------------------------------------------------

def create_one_subgrid_wavelength_grid(npoints, emission_lines=None, fixed=None, min_wavelength=None, max_wavelength=None,
                                       filters=None, min_wavelengths_in_filter=5, min_wavelengths_in_fwhm=3):

    """
    This function ...
    :param npoints:
    :param emission_lines:
    :param fixed:
    :param min_wavelength:
    :param max_wavelength:
    :param filters:
    :param min_wavelengths_in_filter:
    :param min_wavelengths_in_fwhm:
    :return:
    """

    # Debugging
    log.debug("Creating wavelength grid with " + str(npoints) + " points ...")

    # A list of the wavelength points
    wavelengths = []

    # Keep track of the number of points per subgrid
    subgrid_npoints = dict()

    # Loop over the subgrids
    for subgrid in subgrids:

        # Debugging
        log.debug("Adding the " + subgrid + " subgrid ...")

        # Determine minimum, maximum
        min_lambda = ranges[subgrid].min
        max_lambda = ranges[subgrid].max

        # Skip subgrids out of range
        if min_wavelength is not None and max_lambda < min_wavelength: continue
        if max_wavelength is not None and min_lambda > max_wavelength: continue

        # Determine the normal number of wavelength points for this subgrid
        points = int(round(relpoints[subgrid] * npoints))

        # Correct the number of wavelengths based on the given min and max wavelength
        #if min_wavelength is not None: min_lambda = max(min_lambda, min_wavelength)
        #if max_wavelength is not None: max_lambda = min(max_lambda, max_wavelength)

        # Generate and add the wavelength points
        wavelengths_subgrid_original = make_grid(min_lambda, max_lambda, points)

        # Filter based on given boundaries
        wavelengths_subgrid = []
        for wav in wavelengths_subgrid_original:
            if min_wavelength is not None and wav < min_wavelength: continue
            if max_wavelength is not None and wav > max_wavelength: continue
            wavelengths_subgrid.append(wav)

        # Set the number of points for this subgrid
        subgrid_npoints[subgrid] = len(wavelengths_subgrid)

        # Add the wavelength points
        wavelengths += wavelengths_subgrid

    # Loop over the filters
    #broad_resampled = 0
    #narrow_added = 0
    broad_resampled = []
    narrow_added = []
    if filters is not None:

        # Debugging
        log.debug("Adding wavelengths for sampling filter bandpasses ...")

        # Loop over the filters
        for fltr in filters:

            # Debugging
            log.debug("Adding wavelength(s) for the " + str(fltr) + " filter ...")

            # Broad band filter: make sure there are at least 10 wavelength points in the range min_wavelength > max_wavelength
            if isinstance(fltr, BroadBandFilter):

                min_wavelength = fltr.min
                max_wavelength = fltr.max

                # Check that at least 5 wavelength points sample the range of the filter
                # Get the indices of the wavelengths that fall within this range in the current list of wavelengths
                current_indices = [i for i in range(len(wavelengths)) if min_wavelength < wavelengths[i] < max_wavelength]

                # Check if there at least 10
                if len(current_indices) >= min_wavelengths_in_filter: continue

                # Otherwise, delete the current wavelengths and add 10 new ones
                for index in sorted(current_indices, reverse=True): del wavelengths[index]

                # One more filter for which we have resampled
                broad_resampled.append(str(fltr))

                # Generate new wavelengths for sampling the filter range on a logarithmic grid
                new_wavelengths = fltr.range.log(min_wavelengths_in_filter, as_list=True)

                # Add the new wavelengths
                wavelengths += new_wavelengths

                # Sort the wavelength points
                wavelengths = sorted(wavelengths)

                # If FWHM of the filter is defined
                if fltr.fwhm is not None:

                    min_wavelength_fwhm = fltr.mean - fltr.fwhm
                    max_wavelength_fwhm = fltr.mean + fltr.fwhm

                    # Check that at least 3 wavelength points sample the inner range of the filter
                    current_indices = [i for i in range(len(wavelengths)) if min_wavelength_fwhm < wavelengths[i] < max_wavelength_fwhm]

                    # Check if there are at least 3
                    if len(current_indices) >= min_wavelengths_in_fwhm: continue

                    # Otherwise, delete the current wavelengths and add 3 new ones
                    for index in sorted(current_indices, reverse=True): del wavelengths[index]

                    # Generate new wavelengths for sampling the inner filter range on a logarithmic grid
                    new_wavelengths = fltr.fwhm_range.log(min_wavelengths_in_fwhm, as_list=True)

                    # Add the new wavelengths
                    wavelengths += new_wavelengths

                    # Sorting is done at the end

            # For a narrow band filter, add the exact wavelength of the filter to the wavelength grid
            elif isinstance(fltr, NarrowBandFilter):

                # Add the wavelength
                wavelengths.append(fltr.wavelength)

                # One more filter for which we have added a wavelength
                #narrow_added += 1
                narrow_added.append(str(fltr))

            # Unrecognized filter
            else: raise ValueError("Unrecognized filter object: " + str(fltr))

    # Add the emission lines
    emission_npoints = 0
    if emission_lines is not None:

        # Debugging
        log.debug("Adding the emission lines ...")

        # Add the mission lines
        wavelengths = add_emission_lines(wavelengths, emission_lines, min_wavelength, max_wavelength)
        emission_npoints = len(emission_lines)

    # Add fixed wavelength points
    fixed_npoints = 0
    if fixed is not None:

        # Debugging
        log.debug("Adding the fixed points to the grid ...")

        fixed_npoints = len(fixed)
        for wavelength in fixed:

            # Debugging
            log.debug("Adding fixed wavelength " + str(wavelength) + " ...")

            if min_wavelength is not None and wavelength < min_wavelength: continue
            if max_wavelength is not None and wavelength > max_wavelength: continue
            wavelengths.append(wavelength)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid and some information about the subgrids
    return grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added

# -----------------------------------------------------------------

def create_one_logarithmic_wavelength_grid(wrange, npoints, emission_lines=None, fixed=None):

    """
    This function ...
    :param wrange:
    :param npoints:
    :param emission_lines:
    :param fixed:
    :return:
    """

    # Verify the grid parameters
    if npoints < 2: raise ValueError("the number of points in the grid should be at least 2")
    if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")

    # Calculate log of boundaries
    logmin = np.log10(float(wrange.min))
    logmax = np.log10(float(wrange.max))

    # Calculate the grid points
    wavelengths = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)

    # Add the emission lines
    emission_npoints = 0
    if emission_lines is not None:

        # Add the mission lines
        wavelengths = add_emission_lines(wavelengths, emission_lines)
        emission_npoints = len(emission_lines)

    # Add fixed wavelength points
    fixed_npoints = 0
    if fixed is not None:
        fixed_npoints = len(fixed)
        for wavelength in fixed: wavelengths.append(wavelength)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid
    return grid, emission_npoints, fixed_npoints

# -----------------------------------------------------------------

def create_one_nested_log_wavelength_grid(wrange, npoints, wrange_zoom, npoints_zoom, emission_lines=None, fixed=None):

    """
    This function ...
    :param wrange:
    :param npoints:
    :param wrange_zoom:
    :param npoints_zoom:
    :param emission_lines:
    :param fixed:
    :return:
    """

    # Verify the grid parameters
    if npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
    if npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
    if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")
    if (wrange_zoom.min <= wrange.min
        or wrange_zoom.max <= wrange_zoom.min
        or wrange.max <= wrange_zoom.max):
        raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

    logmin = np.log10(float(wrange.min))
    logmax = np.log10(float(wrange.max))
    logmin_zoom = np.log10(float(wrange_zoom.min))
    logmax_zoom = np.log10(float(wrange_zoom.max))

    # Build the high- and low-resolution grids independently
    base_grid = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)
    zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=npoints_zoom, endpoint=True, base=10)

    # Merge the two grids
    wavelengths = []

    # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
    for wavelength in base_grid:
        if wavelength < wrange_zoom.min: wavelengths.append(wavelength)

    # Add the wavelengths of the high-resolution grid
    for wavelength in zoom_grid: wavelengths.append(wavelength)

    # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
    for wavelength in base_grid:
        if wavelength > wrange_zoom.max: wavelengths.append(wavelength)

    # Add the emission lines
    emission_npoints = 0
    if emission_lines is not None:

        # Add the mission lines
        wavelengths = add_emission_lines(wavelengths, emission_lines)
        emission_npoints = len(emission_lines)

    # Add fixed wavelength points
    fixed_npoints = 0
    if fixed is not None:
        fixed_npoints = len(fixed)
        for wavelength in fixed: wavelengths.append(wavelength)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid
    return grid, emission_npoints, fixed_npoints

# -----------------------------------------------------------------

def add_emission_lines(wavelengths, emission_lines, min_wavelength=None, max_wavelength=None):

    """
    This function ...
    :param wavelengths:
    :param emission_lines:
    :param min_wavelength:
    :param max_wavelength:
    :return:
    """

    # Add emission line grid points
    logdelta = 0.001
    for line in emission_lines:

        center_micron = line.center.to("micron").value
        left_micron = line.left.to("micron").value
        right_micron = line.right.to("micron").value

        if min_wavelength is not None and line.center < min_wavelength: continue
        if max_wavelength is not None and line.center > max_wavelength: continue

        # logcenter = np.log10(center)
        logleft = np.log10(left_micron if left_micron > 0 else center_micron) - logdelta
        logright = np.log10(right_micron if right_micron > 0 else center_micron) + logdelta
        newgrid = []
        for w in wavelengths:
            logw = np.log10(w.to("micron").value)
            if logw < logleft or logw > logright:
                newgrid.append(w)
        newgrid.append(line.center)
        if left_micron > 0:
            newgrid.append(line.left)
        if right_micron > 0:
            newgrid.append(line.right)
        wavelengths = newgrid

    # Return the new wavelength list
    return wavelengths

# -----------------------------------------------------------------

def make_grid(wmin, wmax, N):

    """
    This function returns a wavelength grid (in micron) with a given resolution (nr of points per decade)
    # in the specified range (in micron), aligned with the 10^n grid points.
    """

    result = []

    wmin = wmin.to("micron").value
    wmax = wmax.to("micron").value

    # generate wavelength points p on a logarithmic scale with lambda = 10**p micron
    #  -2 <==> 0.01
    #   4 <==> 10000
    for i in range(-2*N,4*N+1):
        p = float(i)/N
        w = 10.**p
        if wmin <= w < wmax: result.append(w * u("micron"))

    # Return the grid
    return result

# -----------------------------------------------------------------
