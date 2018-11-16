#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.distribution Contains the Distribution and Distribution2D classes.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from scipy.interpolate import InterpolatedUnivariateSpline

# Import astronomical modules
from astropy.modeling import models, fitting
from astropy.io import fits

# Import the relevant PTS classes and modules
from ..tools.utils import lazyproperty
from .curve import Curve
from .log import log
from .table import SmartTable

# -----------------------------------------------------------------

nan_value = float("nan")

# -----------------------------------------------------------------

class Distribution(Curve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        if "name" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:

            # Get x properties
            name = kwargs.pop("name")
            unit = kwargs.pop("unit", None)
            description = kwargs.pop("description", None)

            # Get y properties
            y_name = kwargs.pop("y_name", "Frequency")
            y_description = kwargs.pop("y_description", None)

            # Set x properties
            x_unit = unit
            x_name = name
            x_description = description

            # Add settings to dict
            kwargs["x_unit"] = x_unit
            kwargs["x_name"] = x_name
            kwargs["y_name"] = y_name
            kwargs["x_description"] = x_description
            kwargs["y_description"] = y_description

        # Call the constructor of the base class
        super(Distribution, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def name(self):
        return self.x_name

    # -----------------------------------------------------------------

    @classmethod
    def from_data(cls, name, data, **kwargs):

        """
        This function ...
        :param name:
        :param data:
        :param kwargs:
        :return:
        """

        from ..tools import numbers

        # Get the data
        if isinstance(data, np.ndarray): pass
        elif hasattr(data, "data"): data = data.data
        else: data = np.asarray(data)

        # Get a mask of the finite values
        #finite = np.logical_not(np.logical_or(np.isnan(data), np.isinf(data)))
        finite = np.isfinite(data)

        # Get the finite values
        values = data[finite]

        # Get min or max value
        min_value = kwargs.pop("min_value", None)
        max_value = kwargs.pop("max_value", None)
        if min_value is not None and max_value is not None: mask = (values >= min_value) * (values <= max_value)
        elif min_value is not None: mask = values >= min_value
        elif max_value is not None: mask = values <= max_value
        else: mask = None

        # Mask
        if mask is not None: values = values[mask]

        # Logscale?
        logarithmic = kwargs.get("logarithmic", False)

        # Sigma clip
        sigma_clip = kwargs.pop("sigma_clip", False)
        if sigma_clip:
            from .log import log
            noriginal = len(values)
            sigma_level = kwargs.pop("sigma_level", 3.)
            values, nmasked = numbers.sigma_clip(values, sigma_level=sigma_level, return_nmasked=True, logarithmic=logarithmic)
            log.debug("Masked " + str(nmasked) + " out of " + str(noriginal) + " values by sigma-clipping")

        # Create the distribution
        return cls.from_values(name, values, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_values(cls, name, values, nbins=20, weights=None, unit=None, logarithmic=False, density=True,
                    sigma_clip=False, sigma_level=3, clip_above=None, clip_below=None, ignore_value=None):

        """
        This function ...
        :param name:
        :param values:
        :param nbins:
        :param weights:
        :param unit:
        :param logarithmic:
        :param density:
        :param sigma_clip:
        :param sigma_level:
        :param clip_above:
        :param clip_below:
        :param ignore_value:
        :return:
        """

        from ..tools import sequences
        from ..tools import numbers

        # Sigma-clip?
        if sigma_clip: values = numbers.sigma_clip(values, sigma_level=sigma_level)

        # Clip below or above?
        if clip_above is not None: values = sequences.clip_above(values, clip_above)
        if clip_below is not None: values = sequences.clip_below(values, clip_below)

        # Ignore value?
        if ignore_value is not None: values = sequences.removed_item(values, ignore_value)

        # Check whether has units
        if sequences.have_units(values):
            if unit is not None: values = sequences.without_units(values, unit=unit)
            else:
                unit = sequences.get_unit(values)
                values = sequences.without_units(values)

        # Check size
        nvalues = len(values)
        if weights is not None and len(weights) != nvalues: raise ValueError("Number of weight values (" + str(len(weights)) + ") should be equal to the number of values (" + str(nvalues) + ")")

        # Create histogram
        if logarithmic:
            logvalues = np.log10(np.array(values))
            logvalues = logvalues[np.isfinite(logvalues)]
            counts, edges = np.histogram(logvalues, bins=nbins, density=density, weights=weights)
            edges = 10**edges
        else: counts, edges = np.histogram(values, bins=nbins, density=density, weights=weights)

        # Determine the centers
        centers = []
        nedges = len(edges)
        for i in range(nedges - 1):
            if logarithmic: center = numbers.geometric_mean(edges[i], edges[i+1])
            else: center = numbers.arithmetic_mean(edges[i], edges[i+1])
            centers.append(center)

        # Create and return
        if density: return cls.from_probabilities(name, counts, centers, unit=unit)
        else: return cls.from_counts(name, counts, centers, unit=unit)

    # -----------------------------------------------------------------

    @classmethod
    def by_rank(cls, name, frequencies, y_name="Frequency", unit=None, sort=False):

        """
        This function ...
        :param name:
        :param frequencies:
        :param y_name:
        :param unit:
        :param sort:
        :return:
        """

        # Get the ranks as values
        nfrequencies = len(frequencies)
        values = np.array(range(nfrequencies), dtype=float)

        # Return
        return cls.from_columns(name, values, frequencies, y_name=y_name, unit=unit, sort=sort)

    # -----------------------------------------------------------------

    @classmethod
    def from_columns(cls, name, values, frequencies, y_name="Frequency", unit=None, sort=False):

        """
        This function ...
        :param name:
        :param values:
        :param frequencies:
        :param y_name:
        :param unit:
        :param sort:
        :return:
        """

        # Create distribution
        distr = cls(name=name, y_name=y_name, unit=unit)

        # Check if sorted
        from ..tools import sequences
        if not sequences.is_sorted(values):
            if sort:
                # Create sorted lists
                indices = sequences.argsort(values)
                values = [values[index] for index in indices]
                frequencies = [frequencies[index] for index in indices]
            else: raise ValueError("Values are not sorted")

        # Set data
        #print(distr.unit)
        #print(name)
        #print(len(distr[name]))
        #print(len(values))

        # BY DOING IT LIKE THIS, THE UNIT OF THE COLUMNS GET LOST (COLUMNS ARE REPLACED BY THE NEW ARRAYS)
        # ALSO: WHEN ACCESSING ANY COLUMN BEFORE THESE LINES, (E.G. print(distr[name]) or self.get_column_unit(name) ...)
        # THESE LINES FAIL BECAUSE THEY ARE THEN INITIALIZED WITH LENGTH 0 AND SETTING ARRAY TO EXISTING COLUMN GIVES INCORRECT SHAPE ERROR
        #distr[name] = np.array(values)
        #distr[y_name] = np.array(frequencies)

        # Add the rows
        for value, frequency in zip(values, frequencies):
            super(Distribution, distr).add_row([value, frequency]) # because the add_row function is prohibited in the Distribution class (because of its lazy features)

        # Return the distribution
        return distr

    # -----------------------------------------------------------------

    @classmethod
    def from_probabilities(cls, name, probabilities, values, unit=None, sort=False):

        """
        This function ...
        :param probabilities:
        :param values:
        :param name:
        :param unit:
        :param sort:
        :return:
        """

        # Create
        return cls.from_columns(name, values, probabilities, y_name="Probability", unit=unit, sort=sort)

    # -----------------------------------------------------------------

    @classmethod
    def from_counts(cls, name, counts, values, unit=None, sort=False):

        """
        This function ...
        :param name:
        :param counts:
        :param values:
        :param unit:
        :param sort:
        :return:
        """

        return cls.from_columns(name, values, counts, y_name="Counts", unit=unit, sort=sort)

    # -----------------------------------------------------------------

    def add_row(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise RuntimeError("Cannot add rows to a distribution object")

    # -----------------------------------------------------------------

    def remove_row(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise RuntimeError("Cannot remove rows from a distribution object")

    # -----------------------------------------------------------------

    def normalized(self, value=1.0, method="max"):

        """
        This function ...
        :param value:
        :param method:
        :return:
        """

        counts = self.frequencies

        # By setting highest bin width to value
        if method == "max":

            max_count = np.max(counts)
            factor = value / max_count
            normalized = counts * factor
            y_name = "Probability"

        # Plain sum of all bin counts
        elif method == "sum":

            total_count = np.sum(counts)
            factor = value / total_count
            normalized = counts * factor
            y_name = "Normalized count"

        # Integral: incorporating bin widths
        elif method == "integral": raise ValueError("Not implemented yet")

        # Invalid
        else: raise ValueError("Invalid option for 'method'")

        # Create new distribution
        return self.__class__.from_columns(self.name, self.values, normalized, y_name=y_name, sort=False)

    # -----------------------------------------------------------------

    @property
    def value_name(self):
        return self.x_name

    # -----------------------------------------------------------------

    @lazyproperty
    def values(self):
        return self.get_x(asarray=True, unit=self.x_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def frequencies(self):
        return self.get_y(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def unit(self):
        return self.column_unit(self.x_name)

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.has_column_unit(self.x_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_value(self):
        from ..tools import numbers
        return numbers.weighed_arithmetic_mean(self.values, weights=self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def mean(self):
        if self.has_unit: return self.mean_value * self.unit
        else: return self.mean_value

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_value(self):
        from ..tools import numbers
        return numbers.weighed_standard_deviation(self.values, weights=self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev(self):
        if self.has_unit: return self.stddev_value * self.unit
        else: return self.stddev_value

    # -----------------------------------------------------------------

    @lazyproperty
    def geometric_mean_value(self):
        from ..tools import numbers
        return numbers.weighed_geometric_mean(self.values, weights=self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def geometric_mean(self):
        if self.has_unit: return self.geometric_mean_value * self.unit
        else: return self.geometric_mean_value

    # -----------------------------------------------------------------

    @lazyproperty
    def median_value(self):
        from ..tools import numbers
        return numbers.median(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def median(self):
        if self.has_unit: return self.median_value * self.unit
        else: return self.median_value

    # -----------------------------------------------------------------

    @lazyproperty
    def percentile_16_value(self):
        return find_percentile_16(self.values, self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def percentile_16(self):
        if self.has_unit: return self.percentile_16_value * self.unit
        else: return self.percentile_16_value

    # -----------------------------------------------------------------

    @lazyproperty
    def percentile_84_value(self):
        return find_percentile_84(self.values, self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def percentile_84(self):
        if self.has_unit: return self.percentile_84_value * self.unit
        else: return self.percentile_84_value

    # -----------------------------------------------------------------

    @lazyproperty
    def fwhm_value(self):
        from ..tools import numbers
        return numbers.fwhm(self.values, self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def fwhm(self):
        if self.has_unit: return self.fwhm_value * self.unit
        else: return self.fwhm_value

    # -----------------------------------------------------------------

    @property
    def nvalues(self):
        return len(self.values)

    # -----------------------------------------------------------------

    @property
    def has_single_value(self):
        return self.nvalues == 1

    # -----------------------------------------------------------------

    @property
    def nbins(self):
        return self.nvalues

    # -----------------------------------------------------------------

    @property
    def has_single_bin(self):
        return self.nbins == 1

    # -----------------------------------------------------------------

    @lazyproperty
    def bin_widths(self):
        widths = []
        for i in range(self.nedges - 1): widths.append(self.edges[i+1] - self.edges[i])
        return widths

    # -----------------------------------------------------------------

    @lazyproperty
    def bin_widths_log(self):
        widths = []
        for i in range(self.nedges - 1): widths.append(self.edges_log[i+1] - self.edges[i])
        return widths

    # -----------------------------------------------------------------

    def get_edges_for_bin(self, index, logscale=False):

        """
        This function ...
        :param index:
        :param logscale:
        :return:
        """

        if logscale: return self.edges_log[index], self.edges_log[index+1]
        else: return self.edges[index], self.edges[index+1]

    # -----------------------------------------------------------------

    @lazyproperty
    def bin_width(self):

        """
        This function ...
        :return:
        """

        # Check
        if not all_close(self.bin_widths): raise RuntimeError("Bin widths not equal")

        # Calculate width
        width = (self.max_value - self.min_value) / self.nbins
        assert np.isclose(width, self.bin_widths[0])

        # Return
        return width

    # -----------------------------------------------------------------

    @lazyproperty
    def bin_width_log(self):

        """
        This function ...
        :return:
        """

        # Check
        #if not all_close(self.bin_widths_log): raise RuntimeError("Bin widths in log space not equal")
        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def most_frequent_index(self):

        """
        This function ...
        :return:
        """

        return np.argmax(self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def most_frequent_value(self):

        """
        This function ...
        :return:
        """

        return self.values[self.most_frequent_index]

    # -----------------------------------------------------------------

    @lazyproperty
    def most_frequent(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.most_frequent_value * self.unit
        else: return self.most_frequent_value

    # -----------------------------------------------------------------

    @lazyproperty
    def least_frequent_index(self):

        """
        This function ...
        :return:
        """

        return np.argmin(self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def least_frequent_value(self):

        """
        This function ...
        :return:
        """

        return self.values[self.least_frequent_index]

    # -----------------------------------------------------------------

    @lazyproperty
    def least_frequent(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.least_frequent_value * self.unit
        else: return self.least_frequent_value

    # -----------------------------------------------------------------

    @lazyproperty
    def least_frequent_value_non_zero(self):

        """
        This function ...
        :return:
        """

        ma = np.ma.masked_equal(self.frequencies, 0.0, copy=False)
        index = np.argmin(ma)
        return self.values[index]

    # -----------------------------------------------------------------

    @lazyproperty
    def least_frequent_non_zero(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.least_frequent_value_non_zero * self.unit
        else: return self.least_frequent_value_non_zero

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_indices(self):

        """
        This function ...
        :return:
        """

        return list(sorted(range(self.nvalues), key=lambda index: self[self.y_name][index]))

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_indices_reversed(self):

        """
        This function ...
        :return:
        """

        return list(reversed(self.sorted_indices))

    # -----------------------------------------------------------------

    def get_leading_indices(self, fraction, return_fraction=False):

        """
        This function ...
        :param fraction:
        :param return_fraction:
        :return:
        """

        # Get the total count/frequency/prob
        total = np.sum(self[self.y_name])

        indices = []
        current_fraction = 0.0

        # Loop over the indices of the rows, from most frequent (highest prob) to lowest
        for index in self.sorted_indices_reversed:

            # Get the count/frequency/prob
            frequency = self[self.y_name][index]
            relative = float(frequency) / total

            # Add index
            indices.append(index)

            # Add to current fraction
            current_fraction += relative

            # Break the loop if we have enough
            if current_fraction >= fraction: break

        # Return the indices
        if return_fraction: return indices, current_fraction
        else: return indices

    # -----------------------------------------------------------------

    def get_leading_values(self, fraction, add_unit=True, return_fraction=False):

        """
        This function ...
        :param fraction:
        :param add_unit:
        :param return_fraction:
        :return:
        """

        # Get the indices
        indices, total_fraction = self.get_leading_indices(fraction, return_fraction=True)

        # Get the values, in the same order
        values = []
        for index in indices:
            value = self.get_value(self.value_name, index, add_unit=add_unit)
            values.append(value)

        # Return the values
        if return_fraction: return values, total_fraction
        else: return values

    # -----------------------------------------------------------------

    def get_leading_edges(self, fraction, logscale=False, add_unit=True, return_fraction=False):

        """
        This function ...
        :param fraction:
        :param logscale:
        :param add_unit:
        :param return_fraction:
        :return:
        """

        # Get the indices
        indices, total_fraction = self.get_leading_indices(fraction, return_fraction=True)

        min_edge = None
        max_edge = None

        # Loop over the indices, get edges
        for index in indices:

            # Get edges
            edges = self.get_edges_for_bin(index, logscale=logscale)

            # Adapt edges
            if min_edge is None or edges[0] < min_edge: min_edge = edges[0]
            if max_edge is None or edges[1] > max_edge: max_edge = edges[1]

        # Add unit?
        if add_unit and self.has_unit:
            min_edge = min_edge * self.unit
            max_edge = max_edge * self.unit

        # Return
        if return_fraction: return min_edge, max_edge, total_fraction
        else: return min_edge, max_edge

    # -----------------------------------------------------------------

    def get_leading_values_and_edges(self, fraction, add_unit=True, logscale=False, return_fraction=False):

        """
        This function ...
        :param fraction:
        :param add_unit:
        :param logscale:
        :param return_fraction:
        :return:
        """

        # Get the indices
        indices, total_fraction = self.get_leading_indices(fraction, return_fraction=True)

        # Get the values, in the same order
        values = []
        for index in indices:
            value = self.get_value(self.value_name, index, add_unit=add_unit)
            values.append(value)

        min_edge = None
        max_edge = None

        # Loop over the indices, get edges
        for index in indices:

            # Get edges
            edges = self.get_edges_for_bin(index, logscale=logscale)

            # Adapt edges
            if min_edge is None or edges[0] < min_edge: min_edge = edges[0]
            if max_edge is None or edges[1] > max_edge: max_edge = edges[1]

        # Add unit?
        if add_unit and self.has_unit:
            min_edge = min_edge * self.unit
            max_edge = max_edge * self.unit

        # Return
        if return_fraction: return values, min_edge, max_edge, total_fraction
        else: return values, min_edge, max_edge

    # -----------------------------------------------------------------

    @property
    def first_value(self):

        """
        This function ...
        :return:
        """

        return self.values[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def first(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.first_value * self.unit
        else: return self.first_value

    # -----------------------------------------------------------------

    @property
    def min_value(self):

        """
        This function ...
        :return:
        """

        return self.first_value

    # -----------------------------------------------------------------

    @property
    def second_value(self):

        """
        This function ...
        :return:
        """

        return self.values[1]

    # -----------------------------------------------------------------

    @lazyproperty
    def second(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.second_value * self.unit
        else: return self.second_value

    # -----------------------------------------------------------------

    @property
    def last_value(self):

        """
        This function ...
        :return:
        """

        return self.values[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def last(self):

        """
        Thisf unction ...
        :return:
        """

        if self.has_unit: return self.last_value * self.unit
        else: return self.last_value

    # -----------------------------------------------------------------

    @property
    def max_value(self):

        """
        This function ...
        :return:
        """

        return self.last_value

    # -----------------------------------------------------------------

    @property
    def second_last_value(self):

        """
        This function ...
        :return:
        """

        return self.values[-2]

    # -----------------------------------------------------------------

    @lazyproperty
    def second_last(self):

        """
        This function ...
        :return:
        """

        if self.has_unit: return self.second_last_value * self.unit
        else: return self.second_last_value

    # -----------------------------------------------------------------

    @property
    def min_frequency(self):

        """
        This function ...
        :return:
        """

        return np.min(self.frequencies)

    # -----------------------------------------------------------------

    @property
    def max_frequency(self):

        """
        This function ...
        :return:
        """

        return np.max(self.frequencies)

    # -----------------------------------------------------------------

    @property
    def min_frequency_nonzero(self):

        """
        This function ...
        :return:
        """

        ma = np.ma.masked_equal(self.frequencies, 0.0, copy=False)
        return np.min(ma)

    # -----------------------------------------------------------------

    @lazyproperty
    def local_maxima(self):

        """
        This function ...
        :return:
        """

        return get_local_maxima(self.values, self.frequencies)

    # -----------------------------------------------------------------

    @lazyproperty
    def local_minima(self):

        """
        This function ...
        :return:
        """

        return get_local_minima(self.values, self.frequencies)

    # -----------------------------------------------------------------

    def less_frequent_than(self, frequency, including=False):

        """
        This function ...
        :param frequency:
        :param including:
        :return:
        """

        return get_less_frequent_than(self.values, self.frequencies, frequency, including=including)

    # -----------------------------------------------------------------

    def more_frequent_than(self, frequency, including=False):

        """
        This function ...
        :param frequency:
        :param including:
        :return:
        """

        return get_more_frequent_than(self.values, self.frequencies, frequency, including=including)

    # -----------------------------------------------------------------

    @lazyproperty
    def smooth(self):

        """
        This function ...
        :return:
        """

        order = 2
        s = InterpolatedUnivariateSpline(self.values, self.frequencies, k=order)

        # Return the spline curve
        return s

    # -----------------------------------------------------------------

    @lazyproperty
    def smooth_log(self):

        """
        This function ...
        :return:
        """

        not_zero = self.frequencies != 0

        centers = self.values[not_zero]
        counts = self.frequencies[not_zero]

        order = 2
        s = InterpolatedUnivariateSpline(centers, np.log10(counts), k=order)

        # Return the spline curve
        return s

    # -----------------------------------------------------------------

    def smooth_values(self, x_min=None, x_max=None, npoints=200):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param npoints:
        :return:
        """

        if x_min is None: x_min = self.min_value
        if x_max is None: x_max = self.max_value

        x_smooth = np.linspace(x_min, x_max, npoints)

        s = self.smooth
        y_smooth = s(x_smooth)

        # Return
        return x_smooth, y_smooth

    # -----------------------------------------------------------------

    def smooth_values_log(self, x_min=None, x_max=None, npoints=200):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param npoints:
        :return:
        """

        if x_min is None: x_min = self.min_value
        if x_max is None: x_max = self.max_value

        x_smooth = np.linspace(x_min, x_max, npoints)

        s = self.smooth_log
        y_smooth_log = s(x_smooth)

        # Return
        return x_smooth, 10.**y_smooth_log

    # -----------------------------------------------------------------

    @lazyproperty
    def local_maxima_smooth(self):

        """
        This function ...
        :return:
        """

        x_smooth, y_smooth = self.smooth_values()
        return get_local_maxima(x_smooth, y_smooth)

    # -----------------------------------------------------------------

    @lazyproperty
    def local_minima_smooth(self):

        """
        This function ...
        :return:
        """

        x_smooth, y_smooth = self.smooth_values()
        return get_local_minima(x_smooth, y_smooth)

    # -----------------------------------------------------------------

    def cumulative_smooth(self, x_min, x_max, npoints=200):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param npoints:
        :return:
        """

        if self._cum_smooth is None:

            x_smooth, y_smooth = self.smooth_values(x_min, x_max, npoints)

            # Set negative values to zero
            y_smooth[y_smooth < 0.0] = 0.0

            # Normalize y by calculating the integral
            #total = simps(y_smooth, x_smooth)

            # NO, by calculating the sum (why?)
            total = np.sum(y_smooth)

            # Now, y should be normalized within x_min:x_max
            y_smooth /= total

            # Return the cumulative distribution
            #return x_smooth, np.cumsum(y_smooth)

            self._cum_smooth = (x_smooth, np.cumsum(y_smooth))

        return self._cum_smooth

    # -----------------------------------------------------------------

    def cumulative_log_smooth(self, x_min, x_max, npoints=200):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param npoints:
        :return:
        """

        x_smooth, y_smooth = self.smooth_values_log(x_min, x_max, npoints)

        # Normalize y by calculating the integral
        #total = simps(y_smooth, x_smooth)

        # NO, by calculating the sum (why?)
        total = np.sum(y_smooth)

        # Now, y should be normalized within x_min:x_max
        y_smooth /= total

        # Return the cumulative distribution
        return x_smooth, np.cumsum(y_smooth)

    # -----------------------------------------------------------------

    def random(self, x_min, x_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :return:
        """

        # Draw a random uniform variate between 0 and 1
        uniform = np.random.uniform(low=0.0, high=1.0)

        x, y_cumulative = self.cumulative_smooth(x_min, x_max)

        from ..tools import nr
        index = nr.locate_clip(y_cumulative, uniform)

        return x[index]

    # -----------------------------------------------------------------

    @lazyproperty
    def edges(self):

        """
        This function ...
        :return:
        """

        from ..tools import numbers

        # Add placeholder for first edge
        edges = [None]

        # Add next edges
        for i in range(self.nvalues - 1):
            edge = numbers.arithmetic_mean(self.values[i], self.values[i+1])
            edges.append(edge)

        # Calculate first and last edge widths
        if self.has_single_bin: first_edge_width = 0.1 * self.values[0]
        else: first_edge_width = self.second_value - self.first_value
        first_edge_half_width = 0.5 * first_edge_width
        if self.has_single_bin: last_edge_width = 0.1 * self.values[0]
        else: last_edge_width = self.last_value - self.second_last_value
        last_edge_half_width = 0.5 * last_edge_width

        # Set first edge position
        first_edge = self.first_value - first_edge_half_width
        edges[0] = first_edge

        # Add last edge position
        last_edge = self.last_value + last_edge_half_width
        edges.append(last_edge)

        # Return the edges
        return edges

    # -----------------------------------------------------------------

    @property
    def nedges(self):
        return len(self.edges)

    # -----------------------------------------------------------------

    @property
    def first_edge(self):
        return self.edges[0]

    # -----------------------------------------------------------------

    @property
    def second_edge(self):
        return self.edges[1]

    # -----------------------------------------------------------------

    @property
    def last_edge(self):
        return self.edges[-1]

    # -----------------------------------------------------------------

    @property
    def second_last_edge(self):
        return self.edges[-2]

    # -----------------------------------------------------------------

    @property
    def min_edge(self):
        return self.first_edge

    # -----------------------------------------------------------------

    @property
    def max_edge(self):
        return self.last_edge

    # -----------------------------------------------------------------

    @lazyproperty
    def edges_log(self):

        """
        This function ...
        :return:
        """

        from ..tools import numbers

        # Add placeholder for first edge
        edges = [None]

        # Add next edges
        for i in range(self.nvalues - 1):
            edge = numbers.geometric_mean(self.values[i], self.values[i + 1])
            edges.append(edge)

        # Calculate first and last edge widths (factors)
        if self.has_single_bin: first_edge_factor = 2.
        else: first_edge_factor = self.second_value / self.first_value
        first_edge_sqrt_factor = np.sqrt(first_edge_factor)
        if self.has_single_bin: last_edge_factor = 2.
        else: last_edge_factor = self.last_value / self.second_last_value
        last_edge_sqrt_factor = np.sqrt(last_edge_factor)

        # Set first edge position
        first_edge = self.first_value / first_edge_sqrt_factor
        edges[0] = first_edge

        # Add last edge position
        last_edge = self.last_value * last_edge_sqrt_factor
        edges.append(last_edge)

        # Return the edges
        return edges

    # -----------------------------------------------------------------

    @property
    def first_edge_log(self):

        """
        This function ...
        :return:
        """

        return self.edges_log[0]

    # -----------------------------------------------------------------

    @property
    def second_edge_log(self):

        """
        This function ...
        :return:
        """

        return self.edges_log[1]

    # -----------------------------------------------------------------

    @property
    def last_edge_log(self):

        """
        Thisn function ...
        :return:
        """

        return self.edges_log[-1]

    # -----------------------------------------------------------------

    @property
    def second_last_edge_log(self):

        """
        Thisn function ...
        :return:
        """

        return self.edges_log[-2]

    # -----------------------------------------------------------------

    @property
    def min_edge_log(self):

        """
        This function ...
        :return:
        """

        return self.first_edge_log

    # -----------------------------------------------------------------

    @property
    def max_edge_log(self):

        """
        This function ...
        :return:
        """

        return self.last_edge_log

    # -----------------------------------------------------------------

    def fit_gaussian(self):

        """
        This function ...
        :return:
        """

        # TODO: this is not working, why ??

        gaussian_init = models.Gaussian1D()
        #fitter = fitting.SLSQPLSQFitter()
        fitter = fitting.LevMarLSQFitter()

        #print(self.centers)
        #print(self.counts)

        #plt.figure()
        #plt.plot(self.centers, self.counts, 'ko')
        #plt.show()

        gaussian_fit = fitter(gaussian_init, self.values, self.frequencies)

        return gaussian_fit

# -----------------------------------------------------------------

def all_equal(array):

    """
    This function ...
    :param array:
    :return:
    """

    first = array[0]

    for i in range(1, len(array)):
        if array[i] != first: return False

    return True

# -----------------------------------------------------------------

def all_close(array):

    """
    This function ...
    :param array:
    :return:
    """

    first = array[0]

    for i in range(1, len(array)):
        if not np.isclose(array[i], first, rtol=1e-4): return False

    return True

# -----------------------------------------------------------------

def find_percentiles(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    if len(values) > 1: return find_percentile_16(values, probabilities), find_percentile_50(values, probabilities), find_percentile_84(values, probabilities)
    else: return None, None, None

# -----------------------------------------------------------------

def find_percentile_16(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    return find_percentile(values, probabilities, 15.86)

# -----------------------------------------------------------------

def find_percentile_50(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    return find_percentile(values, probabilities, 50.)

# -----------------------------------------------------------------

def find_percentile_84(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    return find_percentile(values, probabilities, 84.14)

# -----------------------------------------------------------------

def find_percentile(values, probabilities, percentile):

    """
    This function ...
    :param values:
    :param probabilities:
    :param percentile:
    :return:
    """

    npoints = 10000
    interpfunc = interp1d(values, probabilities, kind='linear')

    parRange = np.linspace(min(values), max(values), npoints)
    interProb = interpfunc(parRange)

    cumInteg = np.zeros(npoints-1)

    for i in range(1,npoints-1):
        cumInteg[i] = cumInteg[i-1] + (0.5*(interProb[i+1] + interProb[i]) * (parRange[i+1] - parRange[i]))

    cumInteg = cumInteg / cumInteg[-1]
    idx = (np.abs(cumInteg-percentile/100.)).argmin()

    return parRange[idx]

# -----------------------------------------------------------------

def get_local_maxima(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    m = argrelextrema(y, np.greater)[0].tolist()

    # Find the index of the absolute maximum (should also be included, is not for example when it is at the edge)
    index = np.argmax(y)
    if index not in m: m.append(index)

    x_maxima = [x[i] for i in m]
    y_maxima = [y[i] for i in m]

    # Return
    return x_maxima, y_maxima

# -----------------------------------------------------------------

def get_local_minima(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    m = argrelextrema(y, np.less)[0].tolist()

    # Find the indx of the absolute minimum (should also be included, is not for example when it is at the edge)
    index = np.argmin(y)
    if index not in m: m.append(index)

    x_minima = [x[i] for i in m]
    y_minima = [y[i] for i in m]

    # Return
    return x_minima, y_minima

# -----------------------------------------------------------------

def get_more_frequent_than(x, y, y0, including=False):

    """
    This function ...
    :param x:
    :param y:
    :param y0:
    :param including:
    :return:
    """

    # Get mask
    if including: where = y >= y0
    else: where = y > y0

    # Return
    return x[where], y[where]

# -----------------------------------------------------------------

def get_less_frequent_than(x, y, y0, including=False):

    """
    Thisf unction ...
    :param x:
    :param y:
    :param y0:
    :param including:
    :return:
    """

    # Get mask
    if including: where = y <= y0
    else: where = y < y0

    # Return
    return x[where], y[where]

# -----------------------------------------------------------------

default_x_name = "x"
default_y_name = "y"

# -----------------------------------------------------------------

class Distribution2D(object):

    """
    This class ...
    """

    def __init__(self, name, x, y, counts, x_name=default_x_name, y_name=default_y_name, x_unit=None,
                 y_unit=None, description=None):

        """
        The constructor ...
        :param name:
        :param x:
        :param y:
        :param counts:
        :param x_name:
        :param y_name:
        :param x_unit:
        :param y_unit:
        :param description:
        """

        from ..units.parsing import parse_unit

        # The name and description
        self.name = name
        self.description = description
        self.x_name = x_name
        self.y_name = y_name

        # Set the data
        self.x = x # X EDGES
        self.y = y # Y EDGES
        self.counts = counts

        # Set the units
        self.x_unit = parse_unit(x_unit) if x_unit is not None else None
        self.y_unit = parse_unit(y_unit) if y_unit is not None else None

        # Check lenghts
        self.check()

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..units.parsing import parse_unit as u
        from ..tools.parsing import real_list

        # Debugging
        log.debug("Loading the 2D histogram in FITS format from '" + path + "' ...")

        # Get data and header
        counts = fits.getdata(path)
        header = fits.getheader(path)

        # Add properties to the header
        name = header["NAME"]
        x_name = header["XNAME"]
        y_name = header["YNAME"]
        x_unit = u(header["XUNIT"]) if "XUNIT" in header else None
        y_unit = u(header["YUNIT"]) if "YUNIT" in header else None
        x = np.array(real_list(header["X"]))
        y = np.array(real_list(header["Y"]))
        description = header["DESCR"] if "DESCR" in header else None

        # Create distribution
        distribution = cls(name, x, y, counts, x_name=x_name, y_name=y_name, x_unit=x_unit, y_unit=y_unit, description=description)

        # Return
        return distribution

    # -----------------------------------------------------------------

    @property
    def has_description(self):
        return self.description is not None

    # -----------------------------------------------------------------

    @property
    def has_x_unit(self):
        return self.x_unit is not None

    # -----------------------------------------------------------------

    @property
    def has_y_unit(self):
        return self.y_unit is not None

    # -----------------------------------------------------------------

    @property
    def nx(self):
        return len(self.x) # NUMBER OF X EDGES

    # -----------------------------------------------------------------

    @property
    def ny(self):
        return len(self.y) # NUMBER OF Y EDGES

    # -----------------------------------------------------------------

    @property
    def nxbins(self):
        return self.nx - 1

    # -----------------------------------------------------------------

    @property
    def nybins(self):
        return self.ny - 1

    # -----------------------------------------------------------------

    @lazyproperty
    def nbins(self):
        if self.nxbins == self.nybins: return self.nxbins
        else: return (self.nxbins, self.nybins,)

    # -----------------------------------------------------------------

    @property
    def counts_nx(self):
        return self.counts.shape[1]

    # -----------------------------------------------------------------

    @property
    def counts_ny(self):
        return self.counts.shape[0]

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        from ..tools import types, sequences

        # Check types
        if not types.is_real_array(self.x): raise ValueError("x data must be an array of real values")
        if not types.is_real_array(self.y): raise ValueError("y data must be an array of real values")
        if not types.is_real_or_integer_array(self.counts): raise ValueError("counts must be an array of real or integer values")

        # Check
        sizes_x = [self.nxbins, self.counts_nx]
        if not sequences.all_equal(sizes_x): raise ValueError("Sizes of x data arrays are not compatible")
        sizes_y = [self.nybins, self.counts_ny]
        if not sequences.all_equal(sizes_y): raise ValueError("Sizes of y data arrays are not compatible")

    # -----------------------------------------------------------------

    @classmethod
    def from_values(cls, name, x, y, weights=None, nbins=200, x_name=default_x_name, y_name=default_y_name,
                    x_unit=None, y_unit=None, description=None, frequency=False):

        """
        This function ...
        :param name:
        :param x:
        :param y:
        :param weights:
        :param nbins:
        :param x_name:
        :param y_name:
        :param x_unit:
        :param y_unit:
        :param description:
        :param frequency:
        :return:
        """

        # Determine the bin step
        #binstep = max(x) / float(nbins)
        #rBins_F, FBins_r = getRadBins(x, y, 1, weights)
        #rBins_F, FBins_r = getRadBins(x, y, binstep, weights)
        #rBins_F[rBins_F > 25] = np.nan
        #rBins_F = None
        #FBins_r = None
        #print("rBins_F", rBins_F)
        #print("FBins_r", FBins_r)

        # Calculate the 2D histogram
        prob, xedges, yedges = np.histogram2d(x, y, bins=nbins, normed=frequency, weights=weights)

        # H needs to be rotated and flipped
        prob = np.rot90(prob)
        prob = np.flipud(prob)

        # Mask zeros
        #zeroes = np.equal(prob, 0)
        #valid = np.logical_not(zeroes)
        #print(prob.shape, len(prob), len(valid), len(xedges), len(yedges))
        #Hmasked = np.ma.masked_where(prob == 0, prob)  # Mask pixels with a value of zero
        #prob = np.ma.masked_where(prob == 0, prob)

        # Clip zeroes
        #xedges = xedges[valid]
        #yedges = yedges[valid]
        #prob = prob[valid]

        # Create the distribution
        # name, x, y, counts, x_name="x", y_name="y", x_unit=None, y_unit=None, description=None
        distribution = cls(name, xedges, yedges, prob, x_name, y_name, x_unit=x_unit, y_unit=y_unit, description=description)

        # Return
        return distribution

    # -----------------------------------------------------------------

    @lazyproperty
    def x_centers(self):

        """
        This function ...
        :return:
        """

        return (self.x[:-1] + self.x[1:]) / 2.

    # -----------------------------------------------------------------

    @lazyproperty
    def y_centers(self):

        """
        Thisf unction ...
        :return:
        """

        return (self.y[:-1] + self.y[1:]) / 2.

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.saveto_fits(path)

    # -----------------------------------------------------------------

    def saveto_fits(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..tools.stringify import tostr

        # Check path
        if not path.endswith("fits"): raise ValueError("Not a valid path for saving as FITS file")

        # Create header with information
        header = fits.Header()

        # Add properties to the header
        header["NAME"] = self.name
        header["XNAME"] = self.x_name
        header["YNAME"] = self.y_name
        if self.has_x_unit: header["XUNIT"] = tostr(self.x_unit)
        if self.has_y_unit: header["YUNIT"] = tostr(self.y_unit)
        header["X"] = tostr(list(self.x))
        header["Y"] = tostr(list(self.y))
        if self.has_description: header["DESCR"] = self.description

        # Create the HDU
        hdu = fits.PrimaryHDU(self.counts, header)

        # Write the HDU to a FITS file
        hdu.writeto(path, clobber=True)

        # Inform the user that the file has been created
        log.debug("File '" + path + "' created")

        # Update path
        self.path = path

    # -----------------------------------------------------------------

    def saveto_table(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Saving the 2D histogram in table format to '" + path + "' ...")

        # Set the columns
        columns = [self.x, self.y]
        names = [self.x_name, self.y_name]

        # Loop over the columns of the counts
        for i in self.counts_nx:
            column = self.counts_nx[i]
            column = np.append(column, nan_value) # to match size
            columns.append(column)
            names.append(self.name + str(i))

        # Capitalize column names
        names = [name.capitalize() for name in names]

        # Set units
        units = [self.x_unit, self.y_unit, None]

        # Debugging
        log.debug("Creating table ...")

        # Create table
        table = SmartTable.from_columns(*columns, names=names, units=units, as_columns=True)

        # Debugging
        log.debug("Writing table ...")

        # Set the description
        if self.has_description: table.meta["description"] = self.description

        # Save the table
        table.saveto(path)

        # Update the path
        self.path = path

# -----------------------------------------------------------------

def getRadBins(xarr, yarr, binstep, weights):

    """
    This function ...
    :param xarr:
    :param yarr:
    :param binstep:
    :param weights:
    :return:
    """

    idx = np.argsort(xarr)

    sortx = xarr[idx]
    sorty = yarr[idx]
    sortweights = weights[idx]

    #avx = np.array([])
    #avy = np.array([])
    avx = []
    avy = []

    i = 0
    counter = 0

    while i*binstep < sortx[-1]:

        n=0
        sumy = 0
        sumweight = 0

        while sortx[counter] < (i+1)*binstep and counter < len(sortx)-1:

            sumy += sorty[counter]*sortweights[counter]
            sumweight += sortweights[counter]
            n += 1
            counter += 1

        if sumweight != 0:

            #avx = np.append(avx,i*binstep+0.5)
            #avy = np.append(avy,sumy / sumweight)

            avx.append(i*binstep+0.5)
            avy.append(sumy / sumweight)

        i += 1

    #return avx, avy
    avx = np.array(avx)
    avy = np.array(avy)
    return avx, avy

# -----------------------------------------------------------------
