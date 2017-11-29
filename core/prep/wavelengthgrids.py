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
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..basics.emissionlines import EmissionLines, EmissionLine
from ..basics.configurable import Configurable
from ..basics.table import SmartTable
from ..basics.range import QuantityRange
from ..units.parsing import parse_unit as u
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..basics.range import IntegerRange
from ..tools import sequences
from ..basics.containers import DefaultOrderedDict
from ..tools.stringify import stringify_list_fancy

# -----------------------------------------------------------------

# The names of the subgrids
#subgrids = ["UV", "optical", "PAH", "dust", "extension"]
subgrids = ["EUV", "stellar", "aromatic", "thermal", "microwave"]

# Define the ranges of the subgrids
ranges = OrderedDict()
ranges["EUV"] = QuantityRange(0.02, 0.085, unit="micron")
ranges["stellar"] = QuantityRange(0.085, 3., unit="micron")
ranges["aromatic"] = QuantityRange(3., 27., unit="micron")
ranges["thermal"] = QuantityRange(27., 1000., unit="micron")
ranges["microwave"] = QuantityRange(1000., 2000, unit="micron")

# Define the relative fineness (the number of points) of the subgrids
relpoints = OrderedDict()
relpoints["EUV"] = 25./325.           # 25
relpoints["stellar"] = 100./325.     # 100
relpoints["aromatic"] = 125./325.         # 125
relpoints["thermal"] = 50./325.         # 50
relpoints["microwave"] = 25./325.    # 25

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
        self.add_column_info("Broad band filters", str, None, "broad band filters for which the wavelength range was resampled")
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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridGenerator, self).__init__(*args, **kwargs)

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
        self.adjust_to = None

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
        self.adjust_to = kwargs.pop("adjust_to", None)

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
        grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added, replaced, new = \
            create_one_subgrid_wavelength_grid(npoints, self.emission_lines, fixed=self.fixed,
                                               min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength,
                                               filters=self.filters, min_wavelengths_in_filter=self.config.min_wavelengths_in_filter,
                                               min_wavelengths_in_fwhm=self.config.min_wavelengths_in_fwhm, adjust_to=self.adjust_to)

        has_replaced = len(replaced) > 0
        has_new = len(new) > 0

        # Debugging
        log.debug("Generated a wavelength grid with:")
        log.debug("")
        log.debug(" - number of points: " + str(len(grid)))
        log.debug(" - number of points in subgrids: ")
        for subgrid in subgrid_npoints: log.debug("     * " + subgrid + ": " + str(subgrid_npoints[subgrid]))
        log.debug(" - number of emission points: " + str(emission_npoints))
        log.debug(" - number of fixed points: " + str(fixed_npoints))
        log.debug(" - filters for which extra sampling was performed: " + str(broad_resampled))
        log.debug(" - narrow band filters for which wavelength was added: " + str(narrow_added))
        if has_replaced: log.debug(" - replaced wavelengths:")
        for old, new in replaced: log.debug("    * " + str(old) + " -> " + str(new))
        if has_new: log.debug(" - new wavelengths:")
        for line in stringify_list_fancy(new)[1].split("\n"): log.debug("    " + line)
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

def get_subgrid_wavelengths(npoints, min_wavelength=None, max_wavelength=None):

    """
    This function ...
    :param min_wavelength:
    :param max_wavelength:
    :param npoints:
    :return:
    """

    subgrid_wavelengths = OrderedDict()

    # Loop over the subgrids
    for subgrid in subgrids:

        # Debugging
        log.debug("Checking " + subgrid + " subgrid wavelengths ...")

        # Determine minimum, maximum
        min_lambda = ranges[subgrid].min
        max_lambda = ranges[subgrid].max

        # Skip subgrids out of range
        if min_wavelength is not None and max_lambda < min_wavelength: continue
        if max_wavelength is not None and min_lambda > max_wavelength: continue

        # Determine the normal number of wavelength points for this subgrid
        points = int(round(relpoints[subgrid] * npoints))

        # Generate and add the wavelength points
        wavelengths_subgrid_original = make_grid(min_lambda, max_lambda, points)

        # Filter based on given boundaries
        wavelengths_subgrid = []
        for wav in wavelengths_subgrid_original:
            if min_wavelength is not None and wav < min_wavelength: continue
            if max_wavelength is not None and wav > max_wavelength: continue
            wavelengths_subgrid.append(wav)

        # Add the sequence of wavelengths
        subgrid_wavelengths[subgrid] = wavelengths_subgrid

    # Return
    return subgrid_wavelengths

# -----------------------------------------------------------------

def resample_filter_wavelengths(wavelengths, filters, min_wavelengths_in_filter=5, min_wavelengths_in_fwhm=3):

    """
    This function ...
    :param wavelengths:
    :param filters:
    :param min_wavelengths_in_filter:
    :param min_wavelengths_in_fwhm:
    :return:
    """

    # Initialize dictionary for the wavelengths for each filter
    filter_wavelengths = DefaultOrderedDict(list)

    # Debugging
    log.debug("Adding wavelengths for sampling filter bandpasses ...")

    # Loop over the filters
    for fltr in filters:

        # Debugging
        log.debug("Adding wavelength(s) for the " + str(fltr) + " filter ...")

        # Broad band filter: make sure there are at least 10 wavelength points in the range min_wavelength > max_wavelength
        if isinstance(fltr, BroadBandFilter):

            # Get filter minimum and maximum wavelength
            min_wavelength = fltr.min
            max_wavelength = fltr.max

            # Check that at least 5 wavelength points sample the range of the filter
            # Get the indices of the wavelengths that fall within this range in the current list of wavelengths
            current_indices = [i for i in range(len(wavelengths)) if min_wavelength < wavelengths[i] < max_wavelength]

            # Check if there at least ..
            if len(current_indices) >= min_wavelengths_in_filter: continue

            # Otherwise, delete the current wavelengths and add 10 new ones
            for index in sorted(current_indices, reverse=True): del wavelengths[index]

            # One more filter for which we have resampled
            #broad_resampled.append(str(fltr))

            # Generate new wavelengths for sampling the filter range on a logarithmic grid
            new_wavelengths = fltr.range.log(min_wavelengths_in_filter, as_list=True)

            # Add the new wavelengths
            wavelengths += new_wavelengths

            # Sort the wavelength points
            # SORTING IS DONE AT THE END
            #wavelengths = sorted(wavelengths)

            # Add the wavelengths
            filter_wavelengths[fltr].extend(new_wavelengths)

            # If FWHM of the filter is defined
            if fltr.fwhm is not None:

                min_wavelength_fwhm = fltr.mean - fltr.fwhm
                max_wavelength_fwhm = fltr.mean + fltr.fwhm

                # Check that at least 3 wavelength points sample the inner range of the filter
                current_indices = [i for i in range(len(wavelengths)) if min_wavelength_fwhm < wavelengths[i] < max_wavelength_fwhm]
                current_indices_new = [i for i in range(len(new_wavelengths)) if min_wavelength_fwhm < new_wavelengths[i] < max_wavelength_fwhm]

                # Check if there are at least 3
                if len(current_indices) >= min_wavelengths_in_fwhm: continue

                # Otherwise, delete the current wavelengths and add 3 new ones
                for index in sorted(current_indices, reverse=True): del wavelengths[index]

                # Generate new wavelengths for sampling the inner filter range on a logarithmic grid
                new_wavelengths = fltr.fwhm_range.log(min_wavelengths_in_fwhm, as_list=True)

                # Add the new wavelengths
                wavelengths += new_wavelengths

                # SORTING IS DONE AT THE END

                # Add to dictionary
                for index in sorted(current_indices_new, reverse=True): del filter_wavelengths[fltr][index]
                filter_wavelengths[fltr].extend(new_wavelengths)

        # For a narrow band filter, add the exact wavelength of the filter to the wavelength grid
        elif isinstance(fltr, NarrowBandFilter):

            # Add the wavelength
            wavelengths.append(fltr.wavelength)

            # One more filter for which we have added a wavelength
            #narrow_added.append(str(fltr))

            # Add the wavelength
            filter_wavelengths[fltr].append(fltr.wavelength)

        # Unrecognized filter
        else: raise ValueError("Unrecognized filter object: " + str(fltr))

    # Return the filter wavelengths
    return filter_wavelengths

# -----------------------------------------------------------------

def adjust_to_wavelengths(wavelengths, adjust_to):

    """
    This function ...
    :param wavelengths:
    :param adjust_to:
    :return:
    """

    # Debugging
    log.debug("Adjusting wavelength grid to include exact wavelengths ...")

    # Replace dict
    replace_dict = defaultdict(list)

    # Original wavelengths with replaced
    replaced = []

    # New wavelengths
    new = []

    # Loop over the wavelengths
    for wavelength in adjust_to:

        # Debugging
        log.debug("Adjusting the wavelength grid for " + str(wavelength) + "...")

        # Find closest
        index = sequences.find_closest_index(wavelengths, wavelength)

        # Add wavelength
        replace_dict[index].append(wavelength)

    # Loop over the indices for which a replacement will be performed
    for index in replace_dict:

        # Only one replacement
        if len(replace_dict[index]) == 1:

            # Get the new wavelength
            wavelength = replace_dict[index][0]

            # Set
            replaced.append((wavelengths[index], wavelength))

            # Replace the wavelength
            wavelengths[index] = wavelength

        # Multiple wavelengths for one index
        else:

            # Get the wavelengths
            adjust_wavelengths = replace_dict[index]

            # Get the closest wavelength to the previous wavelength
            previous_wavelength = wavelengths[index]

            # Find closest
            adjust_index = sequences.find_closest_index(adjust_wavelengths, previous_wavelength)

            # Set
            replaced.append((wavelengths[index], adjust_wavelengths[adjust_index]))

            # Replace closest
            wavelengths[index] = adjust_wavelengths[adjust_index]

            # Add other
            for j in range(len(adjust_wavelengths)):

                # Already added
                if j == adjust_index: continue

                # Add the wavelength
                wavelengths.append(adjust_wavelengths[j])

                # Add to new
                new.append(adjust_wavelengths[j])

    # Return
    return replaced, new

# -----------------------------------------------------------------

def create_one_subgrid_wavelength_grid(npoints, emission_lines=None, fixed=None, min_wavelength=None, max_wavelength=None,
                                       filters=None, min_wavelengths_in_filter=5, min_wavelengths_in_fwhm=3, adjust_to=None,
                                       return_elements=False):

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
    :param adjust_to:
    :param return_elements:
    :return:
    """

    # Debugging
    log.debug("Creating wavelength grid with " + str(npoints) + " points ...")

    # A list of the wavelength points
    wavelengths = []

    # Get subgrid wavelength sequences
    subgrid_wavelengths = get_subgrid_wavelengths(npoints, min_wavelength=min_wavelength, max_wavelength=max_wavelength)

    # Loop over the subgrids
    for subgrid in subgrid_wavelengths:

        # Debugging
        log.debug("Adding the " + subgrid + " subgrid ...")

        # Add the wavelength points
        wavelengths += subgrid_wavelengths[subgrid]

    # Loop over the filters
    if filters is not None: filter_wavelengths = resample_filter_wavelengths(wavelengths, filters, min_wavelengths_in_filter=min_wavelengths_in_filter, min_wavelengths_in_fwhm=min_wavelengths_in_fwhm)
    else: filter_wavelengths = dict()

    # Adjust to passed wavelengths
    if adjust_to is not None: replaced, new = adjust_to_wavelengths(wavelengths, adjust_to)
    else: replaced = new = []

    # Add the emission lines
    if emission_lines is not None:

        # Debugging
        log.debug("Adding the emission lines ...")

        # Add the mission lines
        wavelengths, line_wavelengths = added_emission_lines(wavelengths, emission_lines, min_wavelength, max_wavelength, return_added=True)

    # No emission lines
    else: line_wavelengths = dict()

    # Add fixed wavelength points
    new_fixed = []
    if fixed is not None:

        # Debugging
        log.debug("Adding the fixed points to the grid ...")

        # Loop over the wavelengths
        for wavelength in fixed:

            # Check if not already there (e.g. from adjust_to)
            if wavelength in wavelengths: continue

            # Debugging
            log.debug("Adding fixed wavelength " + str(wavelength) + " ...")

            if min_wavelength is not None and wavelength < min_wavelength: continue
            if max_wavelength is not None and wavelength > max_wavelength: continue
            wavelengths.append(wavelength)

            # Add to fixed
            new_fixed.append(wavelength)

    #print(wavelengths)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid and some information about the subgrids
    if return_elements: return grid, subgrid_wavelengths, filter_wavelengths, replaced, new, line_wavelengths, new_fixed
    else:

        # Keep track of the number of points per subgrid
        subgrid_npoints = OrderedDict()
        for subgrid in subgrid_wavelengths:
            subgrid_npoints[subgrid] = len(subgrid_wavelengths[subgrid])

        # Set broad resampled nwavelengths
        broad_resampled = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, BroadBandFilter): continue
            broad_resampled.append(fltr)

        # Set narrow band filter added nwavelengths
        narrow_added = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, NarrowBandFilter): continue
            narrow_added.append(fltr)

        # Set emission npoints
        emission_npoints = 0
        for line_identifier in line_wavelengths:
            emission_npoints += len(line_wavelengths[line_identifier])

        # Set fixed npoints
        fixed_npoints = len(new_fixed)

        # Return
        return grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added, replaced, new

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
        wavelengths = added_emission_lines(wavelengths, emission_lines)
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
        wavelengths = added_emission_lines(wavelengths, emission_lines)
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

def added_emission_lines(wavelengths, emission_lines, min_wavelength=None, max_wavelength=None, return_added=False):

    """
    This function ...
    :param wavelengths:
    :param emission_lines:
    :param min_wavelength:
    :param max_wavelength:
    :param return_added:
    :return:
    """

    # Initialize dictionary to contain the wavelengths added to the grid for each emission line
    line_wavelengths = DefaultOrderedDict(list)

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
            if logw < logleft or logw > logright: newgrid.append(w)

        newgrid.append(line.center)
        line_wavelengths[line.identifier].append(line.center)

        if left_micron > 0:
            newgrid.append(line.left)
            line_wavelengths[line.identifier].append(line.left)

        if right_micron > 0:
            newgrid.append(line.right)
            line_wavelengths[line.identifier].append(line.right)

        wavelengths = newgrid

    # Return the new wavelength list
    if return_added: return wavelengths, line_wavelengths
    else: return wavelengths

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
