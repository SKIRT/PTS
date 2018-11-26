#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.datacube Contains the DataCube class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict, OrderedDict
import numpy as np

# Import astronomical modules
from astropy.units import Unit
from astropy.convolution import Gaussian2DKernel

# Import the relevant PTS classes and modules
from .image import Image
from .frame import Frame
from ...core.data.sed import SED, ObservedSED
from ...core.basics.curve import WavelengthCurve
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.log import log
from ..basics.mask import Mask, MaskBase
from ...core.basics.errorbar import ErrorBar
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...core.filter.broad import BroadBandFilter
from ..basics.vector import Pixel
from ...core.tools.parallelization import ParallelTarget
from ...core.tools import types
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr, get_list_string_max_nvalues
from ...core.tools import sequences, numbers
from ...core.units.unit import PhotometricUnit
from ..region.region import PixelRegion, SkyRegion
from .frame import nan_value

# -----------------------------------------------------------------

parallel_filter_convolution_dirname = "datacube-parallel-filter-convolution"

# -----------------------------------------------------------------

def get_parameters_path(output_path=None):

    """
    This function ...
    :param output_path:
    :return:
    """

    # Determine output path
    if output_path is None: output_path = fs.cwd()

    # Get parameters path
    paths = fs.files_in_path(output_path, extension="xml", endswith="parameters")

    # Return single
    return sequences.get_single(paths, method="none")

# -----------------------------------------------------------------

def load_all_skirt_datacube_paths(output_path=None):

    """
    This function ...
    :param output_path:
    :return:
    """

    # Determine output path
    if output_path is None: output_path = fs.cwd()

    # Get datacube paths
    datacube_paths = fs.files_in_path(output_path, extension="fits")

    # Determine prefix
    prefix = None
    for path in datacube_paths:
        filename = fs.strip_extension(fs.name(path))
        if prefix is None: prefix = filename.split("_")[0]
        elif prefix != filename.split("_")[0]: raise IOError("Not all datacubes have the same simulation prefix")
    if prefix is None: raise IOError("No datacubes were found")

    # Arrange the datacubes per instrument and per contribution
    datacube_paths_instruments = defaultdict(dict)
    for path in datacube_paths:
        filename = fs.strip_extension(fs.name(path))
        instrument = filename.split(prefix + "_")[1].split("_")[0]
        contr = filename.split("_")[-1]
        datacube_paths_instruments[instrument][contr] = path

    # Return the paths
    return datacube_paths_instruments

# -----------------------------------------------------------------

def load_skirt_datacube_paths(output_path=None, contribution="total"):

    """
    This function ...
    :param output_path:
    :param contribution:
    :return:
    """

    # All contributions
    all_paths = load_all_skirt_datacube_paths(output_path=output_path)

    # Per instrument, for the specific contribution
    datacube_paths = dict()

    # Loop over the instruments
    for instrument_name in all_paths:

        # Get the path
        path = all_paths[instrument_name][contribution]

        # Add the path to the dictionary
        datacube_paths[instrument_name] = path

    # Return the paths
    return datacube_paths

# -----------------------------------------------------------------

def load_skirt_sed_paths(output_path=None):

    """
    This function ...
    :param output_path:
    :return:
    """

    # Determine output path
    if output_path is None: output_path = fs.cwd()

    # Get SED paths
    sed_paths = fs.files_in_path(output_path, extension="dat", endswith="_sed")

    # Determine prefix
    prefix = None
    for path in sed_paths:
        filename = fs.strip_extension(fs.name(path))
        if prefix is None: prefix = filename.split("_")[0]
        elif prefix != filename.split("_")[0]: raise IOError("Not all SED files have the same simulation prefix")
    if prefix is None: raise IOError("No SED files were found")

    # Arrange the paths per instrument and per contribution
    sed_paths_instruments = dict()
    for path in sed_paths:
        filename = fs.strip_extension(fs.name(path))
        instrument = filename.split(prefix + "_")[1].split("_sed")[0]
        sed_paths_instruments[instrument] = path

    # Return the SED file paths
    return sed_paths_instruments

# -----------------------------------------------------------------

def load_skirt_datacubes(output_path=None, contribution="total", wavelength_range=None):

    """
    This function ...
    :param output_path:
    :param contribution:
    :param wavelength_range:
    :return:
    """

    # Load paths
    datacube_paths_instruments = load_all_skirt_datacube_paths(output_path=output_path)
    sed_paths_instruments = load_skirt_sed_paths(output_path=output_path)

    # Initialize dictionary to contain the datacubes
    datacubes = dict()

    # Loop over the instruments
    for instrument_name in datacube_paths_instruments:

        # Get the datacube path
        datacube_path = datacube_paths_instruments[instrument_name][contribution]

        # Get the SED path
        sed_path = sed_paths_instruments[instrument_name]

        # Load datacube
        datacube = DataCube.from_file_and_sed_file(datacube_path, sed_path, wavelength_range=wavelength_range)

        # Add to dictionary
        datacubes[instrument_name] = datacube

    # Return the dictionary of datacubes
    return datacubes

# -----------------------------------------------------------------

class DataCube(Image):

    """
    This class...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    def __init__(self, name="untitled", distance=None, pixelscale=None, psf_filter=None, fwhm=None):

        """
        The constructor ...
        :param name:
        :param distance:
        :param pixelscale:
        :param psf_filter:
        :param fwhm:
        """

        # Call the constructor of the base class
        super(DataCube, self).__init__(name)

        # The wavelength grid
        self.wavelength_grid = None

        # Set properties
        if distance is not None: self.distance = distance
        if pixelscale is not None: self.pixelscale = pixelscale
        if psf_filter is not None: self.psf_filter = psf_filter
        if fwhm is not None: self.fwhm = fwhm

    # -----------------------------------------------------------------

    @classmethod
    def from_array(cls, array, wavelength_grid, axis=0, wcs=None, unit=None):

        """
        This function ...
        :param array:
        :param wavelength_grid:
        :param axis:
        :param wcs:
        :param unit:
        :return:
        """

        # Get the number of frames
        nframes = array.shape[axis]

        # Initialize list of frames
        frames = []

        # Loop over the frames
        for index in range(nframes):

            # Get frame data
            if axis == 0: data = array[index, :, :]
            elif axis == 1: data = array[:, index, :]
            elif axis == 2: data = array[:, : index]
            else: raise ValueError("'axis' parameter should be integer 0-2")

            # Create the frame
            frame = Frame(data, wcs=wcs, unit=unit)

            # Add the frame
            frames.append(frame)

        # Create and return
        return cls.from_frames(frames, wavelength_grid=wavelength_grid)

    # -----------------------------------------------------------------

    @classmethod
    def from_sed_and_map(cls, sed, frame):

        """
        This fnuction
        :param sed: SED to scale the frame to
        :param frame: map of surface density distribution
        :return:
        """

        # Initialize list of frames
        frames = []

        # Loop over the photometry points
        for index in range(sed.npoints):

            # Get the photometry
            photometry = sed.get_photometry(index)

            # Create scaled frame
            scaled_frame = frame.normalized(to=photometry, silent=True)

            # Add the frame
            frames.append(scaled_frame)

        # Create the wavelength grid
        wavelength_grid = sed.wavelength_grid()

        # Create and return
        return cls.from_frames(frames, wavelength_grid=wavelength_grid)

    # -----------------------------------------------------------------

    @classmethod
    def from_sed_and_map_file(cls, sed_path, frame_path):

        """
        This function ...
        :param sed_path:
        :param frame_path:
        :return:
        """

        from ...core.data.sed import load_sed

        # Load the SED (can be ObservedSED, and can be from SKIRT output)
        sed = load_sed(sed_path)

        # Load the map
        frame = Frame.from_file(frame_path)

        # Create and return
        return cls.from_sed_and_map(sed, frame)

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_output(cls, instrument_name, output_path=None, contribution="total", wavelength_range=None):

        """
        This function ...
        :param instrument_name:
        :param output_path:
        :param contribution:
        :param wavelength_range:
        :return:
        """

        # Load paths
        datacube_paths_instruments = load_all_skirt_datacube_paths(output_path=output_path)
        sed_paths_instruments = load_skirt_sed_paths(output_path=output_path)

        # Get the desired datacube path
        datacube_path = datacube_paths_instruments[instrument_name][contribution]

        # Get the SED path for the instrument
        sed_path = sed_paths_instruments[instrument_name]

        # Get the parameters path
        parameters_path = get_parameters_path(output_path=output_path)
        if parameters_path is not None:
            from ...core.simulation.skifile import SkiFile
            skifile = SkiFile(parameters_path)
            distance = skifile.get_instrument_distance(instrument_name)
        else: distance = None

        # Return
        return cls.from_file_and_sed_file(datacube_path, sed_path, wavelength_range=wavelength_range, distance=distance)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_sed_file(cls, image_path, sed_path, wavelength_range=None, distance=None):

        """
        This function ...
        :param image_path:
        :param sed_path:
        :param wavelength_range:
        :param distance:
        :return:
        """

        from ...core.data.sed import load_sed

        # Load the SED (can be ObservedSED, and can be from SKIRT output)
        sed = load_sed(sed_path)

        # Create
        return cls.from_file_and_sed(image_path, sed, wavelength_range=wavelength_range, distance=distance)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_sed(cls, image_path, sed, wavelength_range=None, distance=None):

        """
        This function ...
        :param image_path:
        :param sed:
        :param wavelength_range:
        :param distance:
        :return:
        """

        # Get the wavelength grid
        wavelength_grid = WavelengthGrid.from_sed(sed)

        # Create the datacube
        return cls.from_file(image_path, wavelength_grid, wavelength_range=wavelength_range, distance=distance)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_wavelength_grid_file(cls, image_path, wavelengths_path, wavelength_range=None, distance=None):

        """
        This function ...
        :param image_path:
        :param wavelengths_path:
        :param wavelength_range:
        :param distance:
        :return:
        """

        # Get the wavelength grid
        wavelength_grid = WavelengthGrid.from_file(wavelengths_path)

        # Create the datacube
        return cls.from_file(image_path, wavelength_grid, wavelength_range=wavelength_range, distance=distance)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, image_path, wavelength_grid=None, wavelength_range=None, distance=None, wcs=None):

        """
        This function ...
        :param image_path:
        :param wavelength_grid:
        :param wavelength_range:
        :param distance:
        :param wcs:
        :return:
        """

        # Get the indices for the given wavelength range
        if wavelength_range is not None:
            if wavelength_grid is None: raise ValueError("Specifying wavelength range without passing the wavelength grid explicitly is not (yet) supported")
            indices = wavelength_grid.wavelength_indices(wavelength_range=wavelength_range)
        else: indices = None

        # Call the corresponding base class function
        datacube = super(DataCube, cls).from_file(image_path, always_call_first_primary=False, no_filter=True,
                                                  density=True, density_strict=True, indices=indices,
                                                  absolute_index_names=False) # IMPORTANT: ASSUME THAT DATACUBES ARE ALWAYS DEFINED IN SPECTRAL DENSITY UNITS!

        # Set the distance if specified
        if distance is not None: datacube.distance = distance

        # Set coordinate system if specified
        if wcs is not None: datacube.wcs = wcs

        # Get the wavelength grid from the metadata if needed
        if wavelength_grid is None:

            from ...core.tools.parsing import real_list, unit

            # Check
            if "wvlngths" not in datacube.metadata: raise ValueError("Wavelengths not specified in header")
            if "wavunit" not in datacube.metadata: raise ValueError("Wavelength unit not specified in header")
            wavelengths = real_list(datacube.metadata.pop("wvlngths"))
            wavelength_unit = unit(datacube.metadata.pop("wavunit"))

            # Create the wavelength grid
            wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit=wavelength_unit)

        # Slice the wavelength grid
        if indices is not None: wavelength_grid = wavelength_grid[indices]

        # Check wavelength grid size
        assert len(wavelength_grid) == datacube.nframes

        # Set the wavelength grid
        datacube.wavelength_grid = wavelength_grid

        # Loop over the frames
        for i in range(datacube.nframes):

            # Frame name
            frame_name = "frame" + str(i)

            # Set the wavelength of the frame
            datacube.frames[frame_name].wavelength = datacube.wavelength_grid[i]

        # Return the datacube instance
        return datacube

    # -----------------------------------------------------------------

    @classmethod
    def from_files(cls, paths, **kwargs):

        """
        This function ...
        :param paths: paths of frames
        :param kwargs:
        :return:
        """

        frames = [Frame.from_file(path) for path in paths]
        return cls.from_frames(frames, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_frames(cls, frames, wavelength_grid=None, is_sorted=False, **kwargs):

        """
        This function ...
        :param frames:
        :param wavelength_grid:
        :param is_sorted:
        :param kwargs:
        :return:
        """

        # Create a datacube instance
        datacube = cls(**kwargs)

        # The indices of the frames, sorted on wavelength
        nframes = len(frames)
        if is_sorted: sorted_indices = range(nframes)
        else: sorted_indices = sorted(range(nframes), key=lambda i: frames[i].wavelength_micron)

        # The list of wavelengths
        wavelengths = []

        # Add the frames
        nframes = 0
        for index in sorted_indices:

            # Add the frame
            frame_name = "frame" + str(nframes)
            datacube.add_frame(frames[index], frame_name, silent=True)

            # Add the wavelength
            if wavelength_grid is None: wavelengths.append(frames[index].wavelength_micron)

            # Increment the number of frames
            nframes += 1

        # Create the wavelength grid
        if wavelength_grid is None: wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit="micron")

        # Set the wavelength grid
        datacube.wavelength_grid = wavelength_grid

        # Return the datacube
        return datacube

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.wavelength_grid.wavelengths(unit=unit, asarray=asarray, add_unit=add_unit)

    # -----------------------------------------------------------------

    def wavelengths_string(self, unit=None, add_unit=True, delimiter=", "):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param delimiter:
        :return:
        """

        # Set unit
        if unit is None: unit = self.wavelength_unit

        # Get scalar values
        wavelengths = self.wavelengths(unit=unit, add_unit=False)
        wavelength_strings = [tostr(wavelength, scientific=False) for wavelength in wavelengths]

        # Create string
        string = delimiter.join(wavelength_strings)

        # Add unit?
        if add_unit: string += " " + tostr(unit)

        # Return the string
        return string

    # -----------------------------------------------------------------

    def wavelength_deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.wavelength_grid.deltas(unit=unit, asarray=asarray, add_unit=add_unit)

    # -----------------------------------------------------------------

    def wavelength_indices(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param include_min:
        :param include_max:
        :return:
        """

        return self.wavelength_grid.wavelength_indices(min_wavelength, max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------

    def get_indices(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):
        return self.wavelength_indices(min_wavelength=min_wavelength, max_wavelength=max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------

    def splice(self, min_wavelength=None, max_wavelength=None, copy=False):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param copy: copy the frames
        :return:
        """

        # Get the indices
        #indices = self.get_indices(x_min=x_min, x_max=x_max)
        wavelength_grid, indices = self.wavelength_grid.get_subgrid(min_wavelength=min_wavelength,
                                                                    max_wavelength=max_wavelength, return_indices=True,
                                                                    include_min=True, include_max=True)

        # Get the frames
        if copy: frames = [self.frames[index].copy() for index in indices]
        else: frames = [self.frames[index] for index in indices]

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=wavelength_grid, distance=self.distance)

    # -----------------------------------------------------------------

    def truncate(self, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        # Get the indices
        wavelength_grid, indices = self.wavelength_grid.get_subgrid(min_wavelength=min_wavelength,
                                                                    max_wavelength=max_wavelength, return_indices=True,
                                                                    include_min=True, include_max=True)

        # Determine the frame names
        #print(indices)
        frame_names = ["frame" + str(index) for index in indices]
        #print(frame_names)

        # Remove the frames
        self.remove_frames_except(*frame_names, silent=True)

        # Set the new wavelength grid
        self.wavelength_grid = wavelength_grid

        # Rename the frames
        for index in range(self.nframes):
            current_name = self.frame_names[index]
            new_name = "frame" + str(index)
            self.rename_frame(current_name, new_name, keep_position=True, silent=True)

    # -----------------------------------------------------------------

    def truncated(self, min_wavelength=None, max_wavelength=None, copy=False):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param copy:
        :return:
        """

        # Get the indices
        wavelength_grid, indices = self.wavelength_grid.get_subgrid(min_wavelength=min_wavelength,
                                                                    max_wavelength=max_wavelength, return_indices=True,
                                                                    include_min=True, include_max=True)

        # Determine the frame names
        frame_names = ["frame" + str(index) for index in indices]

        # Create new datacube
        new = self.__class__(name=self.name)

        # Add the frames
        for index, frame_name in enumerate(frame_names):
            frame = self.frames[frame_name]
            new_name = "frame" + str(index)
            new.add_frame(frame, new_name, copy=copy, silent=True)

        # Set the wavelength grid
        new.wavelength_grid = wavelength_grid

        # Return the new cube
        return new

    # -----------------------------------------------------------------

    def flatten_above(self, wavelength, flatten_value=0., include=True):

        """
        This function ...
        :param wavelength:
        :param flatten_value:
        :param include:
        :return:
        """

        from ...core.units.parsing import parse_quantity

        # Check whether wavelength is passed
        if hasattr(wavelength, "unit"):
            if not types.is_length_quantity(wavelength): raise ValueError("Not a wavelength: '" + tostr(wavelength) + "'")
        elif types.is_string_type(wavelength): wavelength = parse_quantity(wavelength)
        else: raise ValueError("Invalid argument: must be length quantity or string")

        # Get the indices
        indices = self.get_indices(wavelength, None, include_min=include)

        # Set flatten value with unit
        if self.has_unit:
            if hasattr(flatten_value, "unit"): pass # OK
            elif flatten_value == 0.: pass # OK
            else: raise ValueError("Unit of the flatten value is not defined")
        # Unit of the datacube is not defined
        elif hasattr(flatten_value, "unit") or types.is_string_type(flatten_value): raise ValueError("Unit of the flatten value is defined, but column unit is not")

        # Flatten the frames
        for index in indices:

            # Determine the flatten value for this frame
            if hasattr(flatten_value, "unit"): value = flatten_value.to(self.unit, wavelength=self.get_wavelength(index), distance=self.distance, pixelscale=self.pixelscale).value
            else: value = flatten_value

            # Set the value
            self.frames[index].fill(value)

    # -----------------------------------------------------------------

    def flatten_below(self, wavelength, flatten_value=0., include=True):

        """
        This function ...
        :param wavelength:
        :param flatten_value:
        :param include:
        :return:
        """

        from ...core.units.parsing import parse_quantity

        # Check whether wavelength is passed
        if hasattr(wavelength, "unit"):
            if not types.is_length_quantity(wavelength): raise ValueError("Not a wavelength: '" + tostr(wavelength) + "'")
        elif types.is_string_type(wavelength): wavelength = parse_quantity(wavelength)
        else: raise ValueError("Invalid argument: must be length quantity or string")

        # Get the indices
        indices = self.get_indices(None, wavelength, include_max=include)

        # Set flatten value with unit
        if self.has_unit:
            if hasattr(flatten_value, "unit"): pass # OK
            elif flatten_value == 0.: pass # OK
            else: raise ValueError("Unit of the flatten value is not defined")
        # Unit of the datacube is not defined
        elif hasattr(flatten_value, "unit") or types.is_string_type(flatten_value): raise ValueError("Unit of the flatten value is defined, but column unit is not")

        # Flatten the frames
        for index in indices:

            # Determine the flatten value for this frame
            if hasattr(flatten_value, "unit"): value = flatten_value.to(self.unit, wavelength=self.get_wavelength(index), distance=self.distance, pixelscale=self.pixelscale).value
            else: value = flatten_value

            # Set the value
            self.frames[index].fill(value)

    # -----------------------------------------------------------------

    def flattened_above(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_above(value, flatten_value=flatten_value, include=include)
        return new

    # -----------------------------------------------------------------

    def flattened_below(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_below(value, flatten_value=flatten_value, include=include)
        return new

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number or frame
        if types.is_real_or_integer(other) or isinstance(other, Frame): frames = [self.frames[index] + other for index in range(self.nframes)]

        # Quantity
        elif types.is_quantity(other):

            # Get the wavelengths
            wavelengths = self.wavelengths()

            # Create frames
            frames = []
            for index in range(self.nframes):
                frame = self.frames[index]
                new_frame = frame + other.to(frame.unit, wavelength=wavelengths[index]).value
                frames.append(new_frame)

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubs must have an equal number of frames")

            # Create new frames
            # Adding frames is implemented in frame class!
            frames = [self.frames[index] + other.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Return new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    __radd__ = __add__

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number or frame
        if types.is_real_or_integer(other) or isinstance(other, Frame):

            # Add to the frames
            for index in range(self.nframes): self.frames[index] += other

        # Quantity
        elif types.is_quantity(other):

            # Get the wavelengths
            wavelengths = self.wavelengths()

            # Add to the frames
            for index in range(self.nframes):
                self.frames[index] += other.to(self.frames[index].unit, wavelength=wavelengths[index]).value

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubs must have an equal number of frames")

            # Add to the frames
            for index in range(self.nframes): self.frames[index] += other.frames[index]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Return self
        return self

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number or frame
        if types.is_real_or_integer(other) or isinstance(other, Frame): frames = [self.frames[index] - other for index in range(self.nframes)]

        # Quantity
        elif types.is_quantity(other):

            # Get the wavelengths
            wavelengths = self.wavelengths()

            # Create frames
            frames = []
            for index in range(self.nframes):
                new_frame = self.frames[index] - other.to(self.frames[index].unit, wavelength=wavelengths[index]).value
                frames.append(new_frame)

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Create new frames
            # Subtracting frames is implemented in Frame class!
            frames = [self.frames[index] - other.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    def __rsub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number or frame
        if types.is_real_or_integer(other) or isinstance(other, Frame): frames = [other - self.frames[index] for index in range(self.nframes)]

        # Quantity
        elif types.is_quantity(other):

            # Get the wavelengths
            wavelengths = self.wavelengths()

            # Create frames
            frames = []
            for index in range(self.nframes):
                new_frame = other.to(self.frames[index].unit, wavelength=wavelengths[index]).value - self.frames[index]
                frames.append(new_frame)

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Create new frames
            frames = [other.frames[index] - self.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    def __isub__(self, other):

        """
        Thisf unction ...
        :param other:
        :return:
        """

        # Regular number or frame
        if types.is_real_or_integer(other) or isinstance(other, Frame):

            # Subtract the frames
            for index in range(self.nframes): self.frames[index] -= other

        # Quantity
        elif types.is_quantity(other):

            # Get the wavelengths
            wavelengths = self.wavelengths()

            # Subtract frames
            for index in range(self.nframes):
                self.frames[index] -= other.to(self.frames[index].unit, wavelength=wavelengths[index]).value

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Subtract the frames
            for index in range(self.nframes): self.frames[index] -= other.frames[index]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Return self
        return self

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number, quantity, or frame
        if types.is_real_or_integer(other) or types.is_quantity(other) or isinstance(other, Frame): frames = [self.frames[index] * other for index in range(self.nframes)]

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Create new frames
            # Multiplying frames is implemented in Frame class!
            frames = [self.frames[index] * other.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    __rmul__ = __mul__

    # -----------------------------------------------------------------

    def __imul__(self, other):
        
        """
        This function ...
        :param other: 
        :return: 
        """

        # Regular number, quantity, or frame
        if types.is_real_or_integer(other) or types.is_quantity(other) or isinstance(other, Frame):

            # Multiply the frames
            for index in range(self.nframes): self.frames[index] *= other

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Multiply the frames
            for index in range(self.nframes): self.frames[index] *= other.frames[index]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Return self
        return self

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number, quantity, or frame
        if types.is_real_or_integer(other) or types.is_quantity(other) or isinstance(other, Frame): frames = [self.frames[index] / other for index in range(self.nframes)]

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Create new frames
            # Dividing frames is implemented in Frame class!
            frames = [self.frames[index] / other.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    def __rdiv__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        #return self.inverse() * other

        # Regular number, quantity or frame
        if types.is_real_or_integer(other) or types.is_quantity(other) or isinstance(other, Frame): frames = [other / self.frames[index] for index in range(self.nframes)]

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Create new frames
            frames = [other.frames[index] / self.frames[index] for index in range(self.nframes)]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Create new datacube
        return self.__class__.from_frames(frames, wavelength_grid=self.wavelength_grid.copy(), distance=self.distance)

    # -----------------------------------------------------------------

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    # -----------------------------------------------------------------

    def __idiv__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Regular number, quantity, or Frame
        if types.is_real_or_integer(other) or types.is_quantity(other) or isinstance(other, Frame):

            # Divide all frames
            for index in range(self.nframes): self.frames[index] /= other

        # Datacube
        elif isinstance(other, DataCube):

            # Check the number of frames
            if self.nframes != other.nframes: raise ValueError("Datacubes must have an equal number of frames")

            # Divide all frames
            for index in range(self.nframes): self.frames[index] /= other.frames[index]

        # Invalid
        else: raise ValueError("Invalid argument")

        # Return self
        return self

    # -----------------------------------------------------------------

    __itruediv__ = __idiv__

    # -----------------------------------------------------------------

    def replace_negatives_by_nans(self):

        """
        This function ...
        :return:
        """

        for frame_name in self.frame_names: self.frames[frame_name].replace_negatives_by_nans()

    # -----------------------------------------------------------------

    def replace_infs_by_nans(self):

        """
        This fnuction ...
        :return:
        """

        for frame_name in self.frame_names: self.frames[frame_name].replace_infs_by_nans()

    # -----------------------------------------------------------------

    def replace_by_nans_where_greater_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for frame_name in self.frame_names: self.frames[frame_name].replace_by_nans_where_greater_than(value)

    # -----------------------------------------------------------------

    def cutoff_greater(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for frame_name in self.frame_names: self.frames[frame_name].cutoff_greater(value)

    # -----------------------------------------------------------------

    def get_mask(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        if isinstance(region_or_mask, PixelRegion): mask = region_or_mask.to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, SkyRegion): mask = region_or_mask.to_pixel(self.wcs).to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, Mask): mask = region_or_mask
        else: raise ValueError("Argument must be region or mask")

        # Return
        return mask

    # -----------------------------------------------------------------

    def create_kernel(self, sigma=None, smoothing_factor=None):

        """
        This function ...
        :param sigma:
        :param smoothing_factor:
        :return:
        """

        # Determine sigma
        if sigma is None:

            # Check whether we have the necessary information
            if self.fwhm is None: raise ValueError("FWHM of the frame should be defined or sigma should be passed")
            if self.pixelscale is None: raise ValueError("Pixelscale of the frame is not defined")

            # Get the sigma in pixels
            sigma = self.sigma_pix

            # Smoothing factor
            if smoothing_factor is not None:

                if smoothing_factor < 1: raise ValueError("Smoothing factor cannot be smaller than one")
                log.debug("Original sigma of the frame resolution is " + tostr(sigma) + " pixels")
                log.debug("Interpolated regions will be smoother by a factor of " + str(smoothing_factor))
                sigma = sigma * smoothing_factor

        # Smoothing factor is passed but also sigma
        elif smoothing_factor is not None: raise ValueError("Smoothing factor cannot be passed when sigma is passed: multiply the specified sigma with the desired smoothing factor")

        # Debugging
        log.debug("Creating a kernel with a sigma of " + tostr(sigma) + " pixels ...")

        # We smooth with a Gaussian kernel with stddev passed by the user
        # Create the kernel
        kernel = Gaussian2DKernel(stddev=sigma)

        # Return
        return kernel

    # -----------------------------------------------------------------

    def interpolate(self, region_or_mask, sigma=None, max_iterations=10, plot=False, not_converge="keep",
                    min_max_in=None, smoothing_factor=None, replace_nans=None):

        """
        Thisfunction ...
        :param region_or_mask:
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :param smoothing_factor:
        :param replace_nans:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Create kernel
        kernel = self.create_kernel(sigma=sigma, smoothing_factor=smoothing_factor)

        # Interpolate each frame
        nframes = self.nframes
        for index, frame_name in enumerate(self.frame_names):

            # Get the frame
            frame = self.frames[frame_name]

            # Get a mask of the original NaN pixels
            original_nans = frame.nans

            # Set originally NaN pixels to something else? zero? -> CAN AFFECT THE INTERPOLATION OF NEIGHBOURING PIXELS
            if replace_nans: frame[original_nans] = replace_nans

            # Set nans at masked pixels
            original_values = frame[mask]
            frame[mask] = nan_value

            # Debugging
            log.debug("Interpolating frame " + str(index+1) + " of " + str(nframes) + " ...")

            # Interpolate the nans
            try: frame.interpolate_nans_with_kernel(kernel, max_iterations=max_iterations, plot=plot, not_converge=not_converge, min_max_in=min_max_in)
            except RuntimeError as e:

                # Reset the original values (e.g. infs)
                frame[mask] = original_values

                # Set original Nans Back to Nan
                frame[original_nans] = nan_value

                # Reraise the error
                raise e

            # Set original NaNs back to NaN
            frame[original_nans] = nan_value

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def interpolate_nans(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None,
                     smoothing_factor=None, error_on_max=True, return_nans=False):

        """
        This function ...
        :param sigma: 
        :param max_iterations: 
        :param plot: 
        :param not_converge: 
        :param min_max_in: 
        :param smoothing_factor:
        :param error_on_max:
        :param return_nans:
        :return: 
        """

        # Create kernel
        kernel = self.create_kernel(sigma=sigma, smoothing_factor=smoothing_factor)

        # Initialize image for NaNs
        nans_image = Image("nans") if return_nans else None

        # Interpolate each frame
        nframes = self.nframes
        for index, frame_name in enumerate(self.frame_names):
            log.debug("Interpolating frame " + str(index+1) + " of " + str(nframes) + " ...")
            nans = self.frames[frame_name].interpolate_nans_with_kernel(kernel, plot=plot, max_iterations=max_iterations,
                                                                    not_converge=not_converge, min_max_in=min_max_in,
                                                                    error_on_max=error_on_max)
            if return_nans: nans_image.add_mask(nans, frame_name.replace("frame", "mask"))

        # Return nans?
        if return_nans: return nans_image

    # -----------------------------------------------------------------

    def interpolate_largest_nans(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None,
                                 smoothing_factor=None, error_on_max=True, return_masks=False, replace_nans=None):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :param smoothing_factor:
        :param error_on_max:
        :param return_masks:
        :param replace_nans:
        :return:
        """

        # Create kernel
        kernel = self.create_kernel(sigma=sigma, smoothing_factor=smoothing_factor)

        # Initialize image for masks
        masks_image = Image("largest_nans") if return_masks else None

        # Interpolate each frame
        nframes = self.nframes
        for index, frame_name in enumerate(self.frame_names):

            # Get the frame
            frame = self.frames[frame_name]

            # Get a mask of the original NaN pixels
            original_nans = frame.nans

            # Get the mask
            largest_nans = original_nans.largest()
            other_nans = largest_nans.inverse()
            mask = largest_nans

            # Set originally NaN pixels to something else? zero? -> CAN AFFECT THE INTERPOLATION OF NEIGHBOURING PIXELS
            if replace_nans is not None: frame[original_nans] = replace_nans

            # Set nans at masked pixels
            #original_values = frame[mask]
            frame[mask] = nan_value

            # Debugging
            log.debug("Interpolating frame " + str(index+1) + " of " + str(nframes) + " ...")

            # Interpolate
            self.frames[frame_name].interpolate_nans_with_kernel(kernel, plot=plot, max_iterations=max_iterations, not_converge=not_converge,
                                                                 min_max_in=min_max_in, error_on_max=error_on_max)
            if return_masks: masks_image.add_mask(mask, frame_name.replace("frame", "mask"))

            # Set other NaNs back to NaN
            frame[other_nans] = nan_value

        # Return masks?
        if return_masks: return masks_image

    # -----------------------------------------------------------------

    def interpolate_not_largest_nans(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None,
                                     smoothing_factor=None, error_on_max=True, return_masks=False, replace_nans=None):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :param smoothing_factor:
        :param error_on_max:
        :param return_masks:
        :param replace_nans:
        :return:
        """

        # Create kernel
        kernel = self.create_kernel(sigma=sigma, smoothing_factor=smoothing_factor)

        # Initialize image for masks
        masks_image = Image("not_largest_nans") if return_masks else None

        # Interpolate each frame
        nframes = self.nframes
        for index, frame_name in enumerate(self.frame_names):

            # Get the frame
            frame = self.frames[frame_name]

            # Get a mask of the original NaN pixels
            original_nans = frame.nans

            # Get the mask
            largest_nans = original_nans.largest()
            #mask = largest_nans.inverse()
            mask = original_nans * largest_nans.inverse()

            # Set originally NaN pixels to something else? zero? -> CAN AFFECT THE INTERPOLATION OF NEIGHBOURING PIXELS
            if replace_nans is not None: frame[original_nans] = replace_nans

            # Set nans at masked pixels
            #original_values = frame[mask]
            frame[mask] = nan_value

            # Debugging
            log.debug("Interpolating frame " + str(index+1) + " of " + str(nframes) + " ...")

            # Interpolate
            self.frames[frame_name].interpolate_nans_with_kernel(kernel, plot=plot, max_iterations=max_iterations, not_converge=not_converge, min_max_in=min_max_in, error_on_max=error_on_max)
            if return_masks: masks_image.add_mask(mask, frame_name.replace("frame", "mask"))

            # Set largest NaNs back to NaN
            frame[largest_nans] = nan_value

        # Return masks?
        if return_masks: return masks_image

    # -----------------------------------------------------------------

    def interpolate_nans_special(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None,
                                     smoothing_factor=None, error_on_max=True, return_masks=False, replace_nans=None,
                                    dilation_rank=2, dilation_connectivity=2):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :param smoothing_factor:
        :param error_on_max:
        :param return_masks:
        :param replace_nans:
        :param dilation_rank:
        :param dilation_connectivity:
        :return:
        """

        # Create kernel
        kernel = self.create_kernel(sigma=sigma, smoothing_factor=smoothing_factor)

        # Initialize image for masks
        masks_image = Image("nans_special") if return_masks else None

        # Interpolate each frame
        nframes = self.nframes
        for index, frame_name in enumerate(self.frame_names):

            # Get the frame
            frame = self.frames[frame_name]

            # Get a mask of the original NaN pixels
            original_nans = frame.nans

            # Debugging
            log.debug("Preparing mask for frame " + str(index+1) + " of " + str(nframes) + " ...")

            # Get not-NaN mask and dilate, then fill holes
            not_nans = original_nans.inverse()
            not_nans.dilate_rc(dilation_rank, dilation_connectivity)
            not_nans.fill_holes()
            mask = original_nans * not_nans
            original_not_mask_nans = original_nans * mask.inverse() # ORIGINAL AND NOT MASKED

            # Set originally NaN pixels to something else? zero? -> CAN AFFECT THE INTERPOLATION OF NEIGHBOURING PIXELS
            if replace_nans is not None: frame[original_nans] = replace_nans

            # Set nans at masked pixels, so that they will be interpolated
            frame[mask] = nan_value

            # Debugging
            log.debug("Interpolating frame " + str(index + 1) + " of " + str(nframes) + " ...")

            # Interpolate
            self.frames[frame_name].interpolate_nans_with_kernel(kernel, plot=plot, max_iterations=max_iterations, not_converge=not_converge, min_max_in=min_max_in, error_on_max=error_on_max)
            if return_masks: masks_image.add_mask(mask, frame_name.replace("frame", "mask"))

            # Set original NaNs back to NaN
            frame[original_not_mask_nans] = nan_value

        # Return masks?
        if return_masks: return masks_image

    # -----------------------------------------------------------------

    def integrate(self):

        """
        This function ...
        :return:
        """

        # Check whether (wavelength) spectral density
        if not self.unit.is_spectral_density: raise ValueError("Datacube is not in spectral density units")

        # Is implemented: data of frames is obtained in 'corresponding_wavelength_density_unit'
        #if not self.unit.is_wavelength_density: raise NotImplementedError("Integrating over frequencies is not yet implemented")
        #if self.is_wavelength_density:

        # Get the unit for the spectral photometry
        unit = self.corresponding_wavelength_density_unit
        wavelength_unit = unit.wavelength_unit
        bolometric_unit = unit.corresponding_bolometric_unit

        # Get the wavelength deltas
        deltas = self.wavelength_deltas(unit=wavelength_unit, asarray=True)

        # Get a list of the frame data in the correct units
        data_list = self.get_data(unit=unit)

        # Get arrays
        arrays = [delta*data for delta, data in zip(deltas, data_list)]

        # Calculate the integral
        frame = Frame(np.sum(arrays, axis=0), wcs=self.wcs, distance=self.distance, unit=bolometric_unit,
                      psf_filter=self.psf_filter, fwhm=self.fwhm, pixelscale=self.pixelscale)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def mean_wavelengths(self, unit="micron"):

        """
        This function ...
        :param unit:
        :return:
        """

        # Get array of wavelengths
        wavelengths = self.wavelengths(asarray=True, unit=unit)

        # Calculate data of mean wavelength per pixel
        data = numbers.weighed_arithmetic_mean_numpy(wavelengths, weights=self.asarray(axis=2))

        # Return the frame
        return Frame(data, unit=unit, wcs=self.wcs, distance=self.distance, fwhm=self.fwhm)

    # -----------------------------------------------------------------

    def get_frames(self, copy=False, unit=None):

        """
        This function ...
        :param copy:
        :param unit:
        :return:
        """

        #return self.frames.as_list(copy=copy)

        frames = []

        # Loop over the frames
        for name in self.frames:

            # Get the frame in the desired unit
            if unit is not None: frame = self.frames[name].converted_to(unit, distance=self.distance, silent=True)
            elif copy: frame = self.frames[name].copy()
            else: frame = self.frames[name]

            # Add the frame
            frames.append(frame)

        # Return the list of frames
        return frames

    # -----------------------------------------------------------------

    def get_data(self, copy=False, unit=None):

        """
        This function ...
        :param copy:
        :param unit:
        :return:
        """

        return [frame.data for frame in self.get_frames(copy=copy, unit=unit)]

    # -----------------------------------------------------------------

    def asarray(self, axis=0):

        """
        This function ...
        :param axis:
        :return:
        """

        # Get a list that contains the data frames
        data_list = self.get_data()

        # Stack the frames into a 3D numpy array
        #if axis == 3: return np.dstack(data_list)
        #elif axis == 2: return np.hstack(data_list)
        #elif axis == 1: return np.vstack(data_list)
        #elif axis == 0: return np.stack(data_list)
        #else: raise ValueError("'axis' parameter should be integer 0-3")

        if axis == 0: return np.stack(data_list)
        elif axis == 2: return np.dstack(data_list)
        else: raise ValueError("'axis' parameter should be 0 or 2")

    # -----------------------------------------------------------------

    def get_frame_index_for_wavelength(self, wavelength, return_wavelength=False):

        """
        This function ...
        :param wavelength:
        :param return_wavelength:
        :return:
        """

        return self.wavelength_grid.closest_wavelength_index(wavelength, return_wavelength=return_wavelength)

    # -----------------------------------------------------------------

    def get_frame_index_for_filter(self, fltr, return_wavelength=False):

        """
        This function ...
        :param fltr:
        :param return_wavelength:
        :return:
        """

        return self.get_frame_index_for_wavelength(fltr.wavelength, return_wavelength=return_wavelength)

    # -----------------------------------------------------------------

    def get_frame_name_for_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        index = self.get_frame_index_for_wavelength(wavelength)
        return self.frames.keys()[index]

    # -----------------------------------------------------------------

    def get_frame_for_wavelength(self, wavelength, copy=False, interpolate=False):

        """
        This function ...
        :param wavelength:
        :param copy:
        :param interpolate:
        :return:
        """

        # Interpolation is not yet supported
        if interpolate: raise NotImplementedError("Interpolation is not yet implemented")

        # Get frame index
        index = self.get_frame_index_for_wavelength(wavelength)

        if copy: frame = self.frames[index].copy()
        else: frame = self.frames[index]

        # Set properties
        frame.unit = self.unit
        frame.wcs = self.wcs
        frame.psf_filter = self.psf_filter
        frame.fwhm = self.fwhm
        frame.pixelscale = self.pixelscale
        frame.distance = self.distance

        # Return
        return frame

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        # If the slicing item is a mask
        if isinstance(item, MaskBase) or isinstance(item, Mask):

            # Create a 3D Numpy array containing
            stack = []
            for frame_name in self.frames:
                stack.append(self.frames[frame_name][item])
            #return np.array(stack)
            return stack # return the list of frame slices

        # If the slicing item is a pixel (x,y)
        if isinstance(item, Pixel):

            # Create a 1D Numpy array
            stack = np.zeros(self.nframes)

            # Set the values
            for index, frame_name in enumerate(self.frames.keys()):
                stack[index] = self.frames[item]

            # Return the numpy array
            return stack

        # Not implemented
        elif isinstance(item, slice): raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    def is_identical(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Check the wavelength grid
        if not np.all(self.wavelengths(asarray=True, unit="micron") == other.wavelengths(asarray=True, unit="micron")): return False

        # Call the implementation of the base class Image
        return super(DataCube, self).is_identical(other)

    # -----------------------------------------------------------------

    def local_sed(self, region, min_wavelength=None, max_wavelength=None, errorcube=None):

        """
        This function ...
        :param region:
        :param min_wavelength:
        :param max_wavelength:
        :param errorcube:
        :return:
        """

        # Determine the unit for the SED
        unit = self.corresponding_non_angular_or_intrinsic_area_unit

        # Determine the conversion factor
        conversion_factor = self.unit.conversion_factor(unit, distance=self.distance, pixelscale=self.pixelscale)
        error_conversion_factor = errorcube.unit.conversion_factor(unit, distance=self.distance, pixelscale=self.pixelscale) if errorcube is not None else None

        # Initialize the SED
        sed = ObservedSED(photometry_unit=unit)

        # Create a mask from the region (or shape)
        mask = region.to_mask(self.xsize, self.ysize)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixels that belong to the region
            flux = np.sum(self.frames[frame_name][mask]) * conversion_factor * unit

            # Get error
            if errorcube is not None:

                # Get the error in the region (quadratic sum)
                error = np.sqrt(np.sum(errorcube.frames[frame_name][mask]**2)) * error_conversion_factor * unit
                errorbar = ErrorBar(error)

                # Add an entry to the SED
                sed.add_point(self.frames[frame_name].filter, flux, errorbar)

            # Add an entry to the SED
            else: sed.add_point(self.frames[frame_name].filter, flux)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def pixel_sed(self, x, y, min_wavelength=None, max_wavelength=None, errorcube=None):

        """
        This function ...
        :param x:
        :param y:
        :param min_wavelength:
        :param max_wavelength:
        :param errorcube:
        :return:
        """

        # Determine the unit for the SED
        unit = self.corresponding_non_angular_or_intrinsic_area_unit

        # Determine the conversion factor
        conversion_factor = self.unit.conversion_factor(unit, distance=self.distance, pixelscale=self.pixelscale)
        error_conversion_factor = errorcube.unit.conversion_factor(unit, distance=self.distance, pixelscale=self.pixelscale) if errorcube is not None else None

        # Initialize the SED
        sed = ObservedSED(photometry_unit=unit)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixel
            flux = self.frames[frame_name][y, x] * conversion_factor * unit

            # Get error
            if errorcube is not None:

                # Get the error in the pixel
                error = errorcube.frames[frame_name][y, x] * error_conversion_factor * unit
                errorbar = ErrorBar(error)

                # Add an entry to the SED
                sed.add_point(self.frames[frame_name].filter, flux, errorbar)

            # Add an entry to the SED
            else: sed.add_point(self.frames[frame_name].filter, flux)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def global_sed(self, mask=None, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param mask:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        # Determine the mask
        if types.is_string_type(mask): inverse_mask = self.masks[mask].inverse()
        elif isinstance(mask, Mask): inverse_mask = mask.inverse()
        elif mask is None: inverse_mask = None
        else: raise ValueError("Mask must be string or Mask (or None) instead of " + str(type(mask)))

        # Determine the unit for the SED
        unit = self.corresponding_non_angular_or_intrinsic_area_unit

        # Determine the conversion factor
        #print(self.distance)
        #print(self.pixelscale)
        conversion_factor = self.unit.conversion_factor(unit, distance=self.distance, pixelscale=self.pixelscale)

        # Initialize the SED
        sed = SED(photometry_unit=unit)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Get the wavelength
            wavelength = self.wavelength_grid[index]

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Calculate the total flux
            if inverse_mask is not None: total_flux = np.sum(self.frames[frame_name][inverse_mask]) * conversion_factor * unit
            else: total_flux = self.frames[frame_name].sum() * conversion_factor * unit

            # Add an entry to the SED
            sed.add_point(wavelength, total_flux)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def global_curve(self, name, mask=None, min_wavelength=None, max_wavelength=None, measure="sum", description=None):

        """
        Thisf unction ...
        :param name: name for the y variable
        :param mask:
        :param min_wavelength:
        :param max_wavelength:
        :param measure:
        :param description:
        :return:
        """

        # Determine the mask
        if types.is_string_type(mask): inverse_mask = self.masks[mask].inverse()
        elif isinstance(mask, Mask): inverse_mask = mask.inverse()
        elif mask is None: inverse_mask = None
        else: raise ValueError("Mask must be string or Mask (or None) instead of " + str(type(mask)))

        # Initialize the curve
        curve = WavelengthCurve(y_name=name, y_unit=self.unit, y_description=description)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Get the wavelength
            wavelength = self.wavelength_grid[index]

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Calculate the value
            if measure == "sum":
                if inverse_mask is not None: y_value = np.sum(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].sum()
            elif measure == "mean":
                if inverse_mask is not None: y_value = np.nanmean(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].mean()
            elif measure == "median":
                if inverse_mask is not None: y_value = np.nanmedian(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].median()
            elif measure == "stddev":
                if inverse_mask is not None: y_value = np.nanstd(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].stddev()
            elif measure == "max":
                if inverse_mask is not None: y_value = np.nanmax(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].max
            elif measure == "min":
                if inverse_mask is not None: y_value = np.nanmin(self.frames[frame_name][inverse_mask])
                else: y_value = self.frames[frame_name].min
            else: raise ValueError("Invalid value for 'measure'")

            # Add an entry to the curve
            curve.add_point(wavelength, y_value)

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    def convolve_with_filter(self, fltr, return_wavelength_grid=False):

        """
        This function ...
        :param fltr:
        :param return_wavelength_grid:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with the " + str(fltr) + " filter ...")

        # Convert the datacube to a numpy array where wavelength is the third dimension (index=2)
        array = self.asarray(axis=2)

        # Get the wavelengths of the datacube
        wavelengths = self.wavelengths(asarray=True)

        # Peform the convolution
        frame, wavelength_grid = _do_one_filter_convolution(fltr, wavelengths, array, self.unit, self.wcs)

        # Return
        if return_wavelength_grid: return frame, wavelength_grid
        else: return frame

    # -----------------------------------------------------------------

    def frame_for_filter(self, fltr, convolve=False):

        """
        This function ...
        :param fltr:
        :param convolve:
        :return:
        """

        # Inform the user
        log.info("Getting frame for the " + str(fltr) + " filter ...")

        # Needs spectral convolution?
        if needs_spectral_convolution(fltr, convolve):

            # Debugging
            log.debug("The frame for the " + str(fltr) + " filter will be calculated by convolving spectrally")

            # Create the frame
            frame = self.convolve_with_filter(fltr)

        # No spectral convolution
        else:

            # Debugging
            log.debug("Getting the frame for the " + str(fltr) + " filter ...")

            # Get the index of the wavelength closest to that of the filter
            index, wavelength = self.get_frame_index_for_wavelength(fltr.wavelength, return_wavelength=True)

            # Get a copy of the frame
            frame = self.frames[index].copy()

            # Set the wavelength
            frame.wavelength = wavelength

        # Set filter to be sure
        frame.filter = fltr

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def get_wavelength(self, index):
        return self.wavelength_grid[index]

    # -----------------------------------------------------------------

    def get_wavelengths(self, unit=None, asarray=False, add_unit=True, min_wavelength=None, max_wavelength=None, inclusive=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param min_wavelength:
        :param max_wavelength:
        :param inclusive:
        :return:
        """

        return self.wavelength_grid.wavelengths(unit=unit, asarray=asarray, add_unit=add_unit, min_wavelength=min_wavelength, max_wavelength=max_wavelength, inclusive=inclusive)

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True, convert=None, add_footprint=False, silent=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :param add_footprint:
        """

        # Call the implementation of the base class
        return super(DataCube, self).rebin(reference_wcs, exact=exact, parallel=parallel, convert=convert, add_footprint=add_footprint, silent=silent)

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs, exact=False, parallel=True, convert=None, add_footprint=False, silent=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :param add_footprint:
        :param silent:
        :return:
        """

        # Call the implementation of the base class
        return super(DataCube, self).rebinned(reference_wcs, exact=exact, parallel=parallel, convert=convert, add_footprint=add_footprint, silent=silent)

    # -----------------------------------------------------------------

    def _check_sampling_for_filter_convolution(self, fltr, wavelengths=None, ignore_bad=False,
                                     min_npoints=8, min_npoints_fwhm=5, skip_ignored_bad_convolution=True):

        """
        This function ...
        :param fltr:
        :param wavelengths:
        :return:
        """

        from ...core.misc.fluxes import WavelengthGridError

        # Get the wavelengths
        if wavelengths is None: wavelengths = self.get_wavelengths(unit="micron", add_unit=True)

        # Get the wavelength indices in the ranges
        indices_in_minmax = [i for i in range(len(wavelengths)) if wavelengths[i] in fltr.range.to("micron")]
        indices_in_fwhm = [i for i in range(len(wavelengths)) if wavelengths[i] in fltr.fwhm_range.to("micron")]

        # Get the number of wavelengths in the ranges
        nwavelengths_in_minmax = len(indices_in_minmax)
        nwavelengths_in_fwhm = len(indices_in_fwhm)

        # Too little wavelengths in range
        if nwavelengths_in_minmax < min_npoints:

            # Warning message
            message = "Too few wavelengths within the filter wavelength range (" + str(fltr.min.to("micron").value) + " to " + str(fltr.max.to("micron").value) + " micron) for convolution (" + str(nwavelengths_in_minmax) + ")"

            # Ignore: don't give error
            if ignore_bad:

                # Warning
                log.warning(message)

                # Skip filter?
                if skip_ignored_bad_convolution:

                    log.warning("Skipping the '" + str(fltr) + "' filter ...")
                    return False # CHECK FAILED

                # Use filter
                else:
                    log.warning("Spectral convolution will still be attempted ...")
                    return True # CHECK OK

            # Error message
            else: raise WavelengthGridError(message, filter=fltr)

        # Too little wavelengths in FWHM range
        elif nwavelengths_in_fwhm < min_npoints_fwhm:

            # Warning message
            message = "Too few wavelengths within the filter FWHM wavelength range (" + str(fltr.fwhm_min.to("micron").value) + " to " + str(fltr.fwhm_max.to("micron").value) + " micron) for convolution (" + str(nwavelengths_in_fwhm) + ")"

            # Ignore: don't give error
            if ignore_bad:

                # Warning
                log.warning(message)

                # Skip filter?
                if skip_ignored_bad_convolution:

                    log.warning("Skipping the '" + str(fltr) + "' filter ...")
                    return False # CHECK FAILED

                # Use filter
                else:
                    log.warning("Spectral convolution will still be attempted ...")
                    return True # CHECK OK

            # Error message
            else: raise WavelengthGridError(message, filter=fltr)

        # OK
        else:
            log.debug("Enough wavelengths within the filter range")
            return True # CHECK OK

    # -----------------------------------------------------------------

    def _check_sampling_for_filter_closest(self, fltr, wavelength=None, ignore_bad=False, skip_ignored_bad_closest=True):

        """
        This function ...
        :param fltr:
        :param wavelength:
        :param ignore_bad:
        :param skip_ignored_bad_closest:
        :return:
        """

        from ...core.misc.fluxes import WavelengthGridError

        # Get wavelength?
        if wavelength is None: wavelength = self.get_wavelength(self.get_frame_index_for_filter(fltr))

        # Check grid wavelength in FWHM
        in_fwhm = wavelength in fltr.fwhm_range

        # Check grid wavelength in inner range
        in_inner = wavelength in fltr.inner_range

        # Not in FWHM?
        if not in_fwhm:

            # Warning message
            message = "Wavelength (" + tostr(wavelength) + ") not in the FWHM range (" + tostr(fltr.fwhm_range) + ") of the filter"

            # Ignore bad sampling: don't give error
            if ignore_bad:

                # Warning
                log.warning(message)

                # Skip filter?
                if skip_ignored_bad_closest:

                    log.warning("Skipping the '" + str(fltr) + "' filter ...")
                    return False # CHECK FAILED

                # Use filter
                else:
                    log.warning("The wavelength '" + tostr(wavelength) + "' will still be used to represent the '" + str(fltr) + "' filter ...")
                    return True # CHECK OK

            # Error
            else: raise WavelengthGridError(message, filter=fltr)

        # Not in inner range
        elif not in_inner:

            # Warning message
            message = "Wavelength (" + tostr(wavelength) + ") not in the inner range (" + tostr(fltr.inner_range) + ") of the filter"

            # Ignore bad sampling: don't give error
            if ignore_bad:

                # Warning
                log.warning(message)

                # Skip filter?
                if skip_ignored_bad_closest:
                    log.warning("Skipping the '" + str(fltr) + "' filter ...")
                    #frames.append(None)
                    #continue
                    return False # CHECK FAILED

                # Use filter
                else:
                    log.warning("The wavelength '" + tostr(wavelength) + "' will still be used to represent the '" + str(fltr) + "' filter ...")
                    return True # CHECK OK

            # Error
            else: raise WavelengthGridError(message, filter=fltr)

        # OK
        else:
            log.debug("Wavelength found close to the filter (" + tostr(wavelength) + ")")
            return True # CHECK OK

    # -----------------------------------------------------------------

    def _initialize_frame_for_filter(self, fltr, convolve, used_wavelength_indices=None, check=True, ignore_bad=False,
                                     min_npoints=8, min_npoints_fwhm=5, skip_ignored_bad_convolution=True,
                                     skip_ignored_bad_closest=True, wavelengths=None):

        """
        This function ...
        :param fltr:
        :param convolve:
        :param used_wavelength_indices:
        :param check:
        :param ignore_bad:
        :param min_npoints:
        :param min_npoints_fwhm:
        :param skip_ignored_bad_convolution:
        :param skip_ignored_bad_closest:
        :param wavelengths:
        :return:
        """

        # Needs spectral convolution?
        if needs_spectral_convolution(fltr, convolve):

            # Debugging
            log.debug("The frame for the " + str(fltr) + " filter will be calculated by convolving spectrally")

            # Check
            if check and not self._check_sampling_for_filter_convolution(fltr, wavelengths=wavelengths, ignore_bad=ignore_bad,
                                                                         min_npoints=min_npoints, min_npoints_fwhm=min_npoints_fwhm,
                                                                         skip_ignored_bad_convolution=skip_ignored_bad_convolution):
                # Return no frame, and also to_convolve = False
                return None, False

            # Return no frame, but to_convolve = True
            return None, True

        # No spectral convolution
        else:

            # Debugging
            log.debug("Getting the frame for the " + str(fltr) + " filter ...")

            # Get the index of the wavelength closest to that of the filter
            index = self.get_frame_index_for_filter(fltr)

            # Get the wavelength
            wavelength = self.get_wavelength(index)

            # Check the difference between the filter wavelength and the actual grid wavelength
            if check and not self._check_sampling_for_filter_closest(fltr, wavelength=wavelength, ignore_bad=ignore_bad,
                                                                     skip_ignored_bad_closest=skip_ignored_bad_closest):

                # Return no frame, and also to_convolve = False
                return None, False

            # Check the wavelength index
            if used_wavelength_indices is not None:

                # Already in the dictionary: frame of datacube already used for other filter
                if index in used_wavelength_indices:
                    filters = used_wavelength_indices[index]
                    filter_names = [str(f) for f in filters]
                    log.warning("The frame for the wavelength '" + str(wavelength) + "' has already been used to create the " + ", ".join(filter_names) + " frame(s)")

                # Add the filter for the wavelength index
                used_wavelength_indices[index].append(fltr)

            # Make a copy of the frame
            frame = self.frames[index].copy()

            # Set the filter
            frame.filter = fltr

            # Set the exact frame wavelength
            frame.wavelength = wavelength

            # Return the frame, and to_convolve = False
            return frame, False

    # -----------------------------------------------------------------

    def _initialize_frames_for_filters(self, filters, convolve=False, check=True, ignore_bad=False, min_npoints=8,
                                       min_npoints_fwhm=5, skip_ignored_bad_convolution=True, skip_ignored_bad_closest=True):

        """
        This function ...
        :param filters:
        :param convolve:
        :param check:
        :param ignore_bad:
        :param min_npoints:
        :param min_npoints_fwhm:
        :param skip_ignored_bad_convolution:
        :param skip_ignored_bad_closest:
        :return:
        """

        # Initialize
        frames = []
        for_convolution = []
        used_wavelength_indices = defaultdict(list)

        # Get the array of wavelengths
        wavelengths = self.get_wavelengths(unit="micron", add_unit=True)

        # Loop over the filters
        for fltr in filters:

            # Create initialized frame
            frame, to_convolve = self._initialize_frame_for_filter(fltr, convolve, used_wavelength_indices=used_wavelength_indices,
                                                                   check=check, ignore_bad=ignore_bad, min_npoints=min_npoints,
                                                                   min_npoints_fwhm=min_npoints_fwhm, skip_ignored_bad_convolution=skip_ignored_bad_convolution,
                                                                   skip_ignored_bad_closest=skip_ignored_bad_closest, wavelengths=wavelengths)

            # Add the frame
            frames.append(frame)

            # Add to list for convolution
            if to_convolve: for_convolution.append(fltr)

        # Return
        return frames, for_convolution, used_wavelength_indices

    # -----------------------------------------------------------------

    def _create_convolved_frames(self, filters, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
        :return:
        """

        # Initialize a list to contain the convolved frames
        convolved_frames = []

        # Calculate convolved frames
        nconvolution = len(filters)
        needs_convolution = nconvolution > 0
        if needs_convolution:

            # Debugging
            log.debug(str(nconvolution) + " filters require spectral convolution (" + ", ".join(str(fltr) for fltr in filters) + ")")

            # Make the frames by convolution
            convolved_frames, wavelengths_for_filters = self.convolve_with_filters(filters, nprocesses=nprocesses, check_previous_sessions=check_previous_sessions, return_wavelengths=True)

            # Show which wavelengths are used to create filter frames
            log.debug("Used the following wavelengths for the spectral convolution for the other filters (in micron):")
            log.debug("")
            for fltr in wavelengths_for_filters:
                filter_name = str(fltr)
                #wavelength_strings = [str(wavelength.to("micron").value) for wavelength in wavelengths_for_filters[fltr]]
                string = get_list_string_max_nvalues([wavelength.to("micron").value for wavelength in wavelengths_for_filters[fltr]], 10)  # max 10 values
                log.debug(" - " + fmt.bold + filter_name + fmt.reset_bold + ": " + string + fmt.bold + " (" + str(len(wavelengths_for_filters[fltr])) + ")" + fmt.reset)
            log.debug("")

            # Set the filters just to be sure
            for fltr, frame in zip(filters, convolved_frames): frame.filter = fltr

        # No spectral convolution
        else: log.debug("Spectral convolution will be used for none of the filters")

        # Return
        return convolved_frames

    # -----------------------------------------------------------------

    def frames_for_filters(self, filters, convolve=False, nprocesses=8, check_previous_sessions=False, as_dict=False,
                           check=True, ignore_bad=False, min_npoints=8, min_npoints_fwhm=5,
                           skip_ignored_bad_convolution=True, skip_ignored_bad_closest=True):

        """
        This function ...
        :param filters:
        :param convolve: boolean or sequence of filters for which spectral convolution should be used
        :param nprocesses:
        :param check_previous_sessions:
        :param as_dict:
        :param check:
        :param ignore_bad:
        :param min_npoints:
        :param min_npoints_fwhm:
        :param skip_ignored_bad_convolution:
        :param skip_ignored_bad_closest:
        :return:
        """

        # Inform the user
        log.info("Getting frames for " + str(len(filters)) + " different filters ...")

        # Limit the number of processes to the number of filters
        nprocesses = min(nprocesses, len(filters))

        # Initialize
        frames, for_convolution, used_wavelength_indices = self._initialize_frames_for_filters(filters, convolve=convolve,
                                                                                               check=check,
                                                                                               ignore_bad=ignore_bad,
                                                                                               min_npoints=min_npoints,
                                                                                               min_npoints_fwhm=min_npoints_fwhm,
                                                                                               skip_ignored_bad_convolution=skip_ignored_bad_convolution,
                                                                                               skip_ignored_bad_closest=skip_ignored_bad_closest)

        # Show which wavelengths are used to create filter frames
        if len(used_wavelength_indices) > 0:
            log.debug("Used the following wavelengths of the datacubes to create frames without spectral convolution:")
            log.debug("")
            for index in used_wavelength_indices:
                wavelength = self.get_wavelength(index)
                wavelength_micron = wavelength.to("micron").value
                used_filters = used_wavelength_indices[index]
                filter_names = [str(f) for f in used_filters]
                nfilters = len(filter_names)
                if nfilters == 1: log.debug(" - " + str(wavelength_micron) + " micron: " + filter_names[0])
                else: log.debug(" - " + str(wavelength_micron) + " micron: " + fmt.bold + ", ".join(filter_names) + fmt.reset)
            log.debug("")

        # Create convolved frames
        convolved_frames = self._create_convolved_frames(for_convolution, nprocesses=nprocesses,
                                                         check_previous_sessions=check_previous_sessions)

        # Add the convolved frames to the list of frames
        for fltr, frame in zip(for_convolution, convolved_frames): frames[filters.index(fltr)] = frame

        # Return the list of frames
        if as_dict: return {fltr: frame for fltr, frame in zip(filters, frames)}
        else: return frames

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters, nprocesses=8, check_previous_sessions=False, return_wavelengths=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
        :param return_wavelengths:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with " + str(len(filters)) + " different filters ...")

        # Limit the number of processes to the number of filters
        nprocesses = min(nprocesses, len(filters))

        # PARALLEL EXECUTION
        if nprocesses > 1: return self.convolve_with_filters_parallel(filters, nprocesses=nprocesses, check_previous_sessions=check_previous_sessions, return_wavelengths=return_wavelengths)

        # SERIAL EXECUTION
        else: return self.convolve_with_filters_serial(filters, return_wavelengths=return_wavelengths)

    # -----------------------------------------------------------------

    def convolve_with_filters_serial(self, filters, return_wavelengths=False):

        """
        Thisj function ...
        :param filters:
        :param return_wavelengths:
        :return:
        """

        # Debugging
        log.debug("Convolving the datacube with " + str(len(filters)) + " different filters on one process ...")

        # Initialize list to contain the output frames per filter
        nfilters = len(filters)
        frames = [None] * nfilters

        # Debugging
        log.debug("Converting the datacube into a single 3D array ...")

        # Convert the datacube to a numpy array where wavelength is the third dimension (index=2)
        array = self.asarray(axis=2)

        # Get the array of wavelengths
        wavelengths = self.wavelengths(asarray=True, unit="micron")

        # Wavelengths used for each filter
        wavelengths_for_filters = OrderedDict()

        # Loop over the filters
        for index in range(nfilters):

            # Get the current filter
            fltr = filters[index]

            # Do the filter convolution to produce one frame
            frame, filter_wavelengths = _do_one_filter_convolution(fltr, wavelengths, array, self.unit, self.wcs)

            # Add the frame to the list
            frames[index] = frame

            # Add the wavelengths
            filter_wavelengths = [value * Unit("micron") for value in filter_wavelengths]
            wavelengths_for_filters[fltr] = filter_wavelengths

        # Return the list of resulting frames
        if return_wavelengths: return frames, wavelengths_for_filters
        else: return frames

    # -----------------------------------------------------------------

    def find_previous_filter_convolution(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Searching for the (partial) results of a previous filter convolution session with this datacube ...")

        # Loop over the previous 'datacube-parallel-filter-convolution' sessions in reversed order (last times first)
        for session_path in introspection.find_temp_dirs(startswith=parallel_filter_convolution_dirname):

            # Determine the path to the datacube
            temp_datacube_path = fs.join(session_path, "datacube.fits")

            # Determine the path to the wavelength grid
            temp_wavelengthgrid_path = fs.join(session_path, "wavelengthgrid.dat")

            # Check whether the datacube file is equal by comparing hash, and the same for the wavelength grid file
            if (self.path is not None and fs.equal_files_hash(temp_datacube_path, self.path)) and (self.wavelength_grid.path is not None and fs.equal_files_hash(temp_wavelengthgrid_path, self.wavelength_grid.path)):

                # Match!
                return session_path, temp_datacube_path, temp_wavelengthgrid_path

            # Check by opening the datacube
            else:

                # Open the datacube
                test_datacube = DataCube.from_file_and_wavelength_grid_file(temp_datacube_path, temp_wavelengthgrid_path)

                # Match!
                if self.is_identical(test_datacube): return session_path, temp_datacube_path, temp_wavelengthgrid_path

        # No match found
        return None, None, None

    # -----------------------------------------------------------------

    def convolve_with_filters_parallel(self, filters, nprocesses=8, check_previous_sessions=False, return_wavelengths=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
        :param return_wavelengths:
        :return:
        """

        # Debugging
        log.debug("Convolving the datacube with " + str(len(filters)) + " different filters with " + str(nprocesses) + " parallel processes ...")

        # Initialize list to contain the output frames per filter
        nfilters = len(filters)
        frames = [None] * nfilters

        # Find (intermediate) results from previous filter convolution
        if check_previous_sessions: temp_dir_path, temp_datacube_path, temp_wavelengthgrid_path = self.find_previous_filter_convolution()
        else: temp_dir_path = temp_datacube_path = temp_wavelengthgrid_path = None

        present_frames = None

        # Not found?
        if temp_dir_path is None:

            # Save the datacube to a temporary directory
            temp_dir_path = introspection.create_temp_dir(time.unique_name(parallel_filter_convolution_dirname))

            # Save the datacube
            temp_datacube_path = fs.join(temp_dir_path, "datacube.fits")
            self.saveto(temp_datacube_path)

            # Save the wavelength grid
            temp_wavelengthgrid_path = fs.join(temp_dir_path, "wavelengthgrid.dat")
            self.wavelength_grid.saveto(temp_wavelengthgrid_path)

        # Found
        else:

            # Success
            log.success("The results of a previous filter convolution session with this datacube were found in the '" + fs.name(temp_dir_path) + "' temporary directory")

            # Look which frames are already created in the directory
            result_paths = fs.files_in_path(temp_dir_path, exact_not_name=["datacube", "wavelengthgrid"], extension="fits", sort=int)

            # Create a dictionary of the paths of the already created frames
            present_frames = dict()
            for path in result_paths:
                name = fs.strip_extension(fs.name(path))
                index = int(name)
                present_frames[index] = path

        # Create process pool
        #pool = Pool(processes=nprocesses)

        # Get string for the unit of the datacube
        unitstring = str(self.unit)

        # Parallel execution
        with ParallelTarget(_do_one_filter_convolution_from_file, nprocesses) as target:

            # EXECUTE THE LOOP IN PARALLEL
            for index in range(nfilters):

                # Check whether already present
                if present_frames is not None and index in present_frames:
                    log.success("The convolved frame for the '" + str(filters[index]) + "' filter is already created in a previous session: skipping calculation ...")
                    continue

                # Debugging
                log.debug("Convolving the datacube to create the '" + str(filters[index]) + "' frame [index " + str(index) + "] ...")

                # Get filtername
                fltrname = str(filters[index])

                # Determine path for resulting frame
                result_path = fs.join(temp_dir_path, str(index) + ".fits")

                # Get the current filter
                #pool.apply_async(_do_one_filter_convolution_from_file, args=(temp_datacube_path, temp_wavelengthgrid_path, result_path, unitstring, fltrname,))  # All simple types (strings)

                # Call the target function
                target(temp_datacube_path, temp_wavelengthgrid_path, result_path, unitstring, fltrname)

        # CLOSE AND JOIN THE PROCESS POOL
        #pool.close()
        #pool.join()

        # Check whether they have all been created
        retry = []
        for index in range(nfilters):

            # Determine the path of the resulting frame
            result_path = fs.join(temp_dir_path, str(index) + ".fits")

            # File exists -> OK
            if fs.is_file(result_path): continue

            # Get filter name
            fltrname = str(filters[index])

            # Give warning
            log.warning("The frame for the '" + fltrname + "' filter has not been created")

            # Add to retry list
            retry.append(index)

        #
        if len(retry) > 0:
            log.warning("Frames " + ", ".join(str(index) for index in retry) + " are missing from the output")
            log.warning("corresponding to the following filters: " + ", ".join(str(filters[index]) for index in retry))

        # RETRY SPECIFIC FRAMES
        with ParallelTarget(_do_one_filter_convolution_from_file, nprocesses) as target:

            # Loop over the retry indices
            for index in retry:

                # Debugging
                log.debug("Peforming convolution of the datacube to create the '" + str(filters[index]) + "' frame again [index " + str(index) + "]...")

                # Get filtername
                fltrname = str(filters[index])

                # Determine path for resulting frame
                result_path = fs.join(temp_dir_path, str(index) + ".fits")

                # Call the target function
                target(temp_datacube_path, temp_wavelengthgrid_path, result_path, unitstring, fltrname)

        # Load the resulting frames
        for index in range(nfilters):

            # Determine path of resulting frame
            result_path = fs.join(temp_dir_path, str(index) + ".fits")

            if not fs.is_file(result_path): raise RuntimeError("Something went wrong: frame " + str(index) + " for the '" + str(filters[index]) + "' filter is missing from the output")

            # Inform the user
            log.debug("Loading the frame for filter " + str(filters[index]) + " from '" + result_path + "' ...")

            # Load the frame and set it in the list
            frames[index] = Frame.from_file(result_path)

        # Return the list of resulting frames
        if return_wavelengths:

            wavelengths_for_filters = OrderedDict()

            # Loop over the filters, set the wavelength grid used for convolution
            for fltr in filters:

                # Get the array of wavelengths
                wa = self.wavelengths(asarray=True, unit="micron")
                wb = fltr._Wavelengths

                # create a combined wavelength grid, restricted to the overlapping interval
                w1 = wa[(wa >= wb[0]) & (wa <= wb[-1])]
                w2 = wb[(wb >= wa[0]) & (wb <= wa[-1])]
                w = np.unique(np.hstack((w1, w2)))
                filter_wavelengths = w

                # Add the list of wavelengths
                filter_wavelengths = [value * Unit("micron") for value in filter_wavelengths]
                wavelengths_for_filters[fltr] = filter_wavelengths

            # Return
            return frames, wavelengths_for_filters

        # Return the list of resulting frames
        else: return frames

    # -----------------------------------------------------------------

    def check_unit(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking the units of the frames ...")

        # Get the unit
        unit = sequences.find_first_not_none([self.frames[index].unit for index in range(self.nframes)], return_none=True)
        if unit is None:
            log.warning("Datacube has no unit")
            return

        # Loop over the frames, check the unit
        for index in range(self.nframes):

            # Get the frame unit
            frame_unit = self.frames[index].unit

            # Check
            if frame_unit is None: self.frames[index].unit = unit
            elif frame_unit != unit: raise ValueError("Frame " + str(index+1) + " has a different unit: '" + tostr(frame_unit) + "' instead of '" + tostr(unit) + "'")

    # -----------------------------------------------------------------

    def convert_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                   brightness_strict=False, silent=True):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param silent:
        :return:
        """

        # Parse "to unit": VERY IMPORTANT, BECAUSE DOING SELF.UNIT = TO_UNIT WILL OTHERWISE REPARSE AND WILL BE OFTEN INCORRECT!! (NO DENSITY OR BRIGHTNESS INFO)
        to_unit = PhotometricUnit(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict,
                                  density_strict=density_strict)

        # Already in the correct unit
        if to_unit == self.unit:
            log.debug("Datacube is already in the desired unit")
            return 1.

        # Get the wavelengths
        wavelengths = self.wavelengths(unit="micron", add_unit=True)

        # Set distance
        if distance is None: distance = self.distance

        # Convert the frames
        for index in range(self.nframes):

            # Debugging
            if not silent: log.debug("Converting frame " + str(index+1) + " ...")

            # Get the wavelength
            wavelength = wavelengths[index]

            # Convert the frame
            self.frames[index].convert_to(to_unit, distance=distance, wavelength=wavelength, silent=silent)

    # -----------------------------------------------------------------

    def converted_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                     brightness_strict=False, silent=True):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param silent:
        :return:
        """

        # Create new
        new = self.__class__(name=self.name)

        # Get the wavelengths
        wavelengths = self.wavelengths(unit="micron", add_unit=True)

        # Set distance
        if distance is None: distance = self.distance

        # Add the frames
        for index in range(self.nframes):

            # Get the wavelength
            wavelength = wavelengths[index]

            # Create converted frame, passing the wavelength for (potential) use in the conversion
            frame = self.frames[index].converted_to(to_unit, distance=distance, density=density, brightness=brightness,
                                                    density_strict=density_strict, brightness_strict=brightness_strict,
                                                    wavelength=wavelength, silent=silent)

            # Add the frame
            frame_name = "frame" + str(index)
            new.add_frame(frame, frame_name, silent=True)

        # Create and set the wavelength grid
        new.wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit="micron")

        # Return the new datacube
        return new

    # -----------------------------------------------------------------

    def convert_by_factor(self, factor, new_unit, silent=True):

        """
        This function ...
        :param factor:
        :param new_unit:
        :param silent:
        :return:
        """

        # Loop over the frames
        for i in range(self.nframes):

            # Debugging
            if not silent: log.debug("Converting frame " + str(i + 1) + " ...")

            # Convert
            self.frames[i].convert_by_factor(factor, new_unit)

    # -----------------------------------------------------------------

    def converted_by_factor(self, factor, new_unit):

        """
        This function ...
        :param factor:
        :param new_unit:
        :return:
        """

        # Create new
        new = self.__class__(name=self.name)

        # Get the wavelengths
        wavelengths = self.wavelengths(unit="micron", add_unit=True)

        # Add the frames
        nframes = 0
        for index in range(self.nframes):

            # Create converted frame
            frame = self.frames[index].converted_by_factor(factor, new_unit)

            # Add the frame
            frame_name = "frame" + str(nframes)
            new.add_frame(frame, frame_name, silent=True)

            # Increment the number of frames
            nframes += 1

        # Create and set the wavelength grid
        new.wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit="micron")

        # Return the new datacube
        return new

    # -----------------------------------------------------------------

    def get_conversion_factor(self, new_unit):

        """
        This function ...
        :param new_unit:
        :return:
        """

        # Get the conversion factor to the new unit
        factor = self.unit.conversion_factor(new_unit, distance=self.distance, pixelscale=self.pixelscale)

        # Debugging
        log.debug("The conversion factor is " + str(factor))

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    @property
    def is_brightness(self):
        return self.unit.is_brightness

    # -----------------------------------------------------------------

    @property
    def corresponding_brightness_unit(self):
        return self.unit.corresponding_brightness_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_brightness_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already brightness
        if self.is_brightness: return

        # Inform the user
        log.info("Converting the datacube to the corresponding brightness unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the unit
        unit = self.corresponding_brightness_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert
        self.convert_by_factor(factor, unit)

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_brightness_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already brightness
        if self.is_brightness: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding brightness unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the unit
        unit = self.corresponding_brightness_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Return converted
        return self.converted_by_factor(factor, unit)

    # -----------------------------------------------------------------

    @property
    def is_surface_brightness(self):
        return self.unit.is_surface_brightness

    # -----------------------------------------------------------------

    @property
    def is_intrinsic_brightness(self):
        return self.unit.is_intrinsic_brightness

    # -----------------------------------------------------------------

    @property
    def corresponding_angular_area_unit(self):
        return self.unit.corresponding_angular_area_unit

    # -----------------------------------------------------------------

    @property
    def corresponding_intrinsic_area_unit(self):
        return self.unit.corresponding_intrinsic_area_unit

    # -----------------------------------------------------------------

    @property
    def corresponding_non_brightness_unit(self):
        return self.unit.corresponding_non_brightness_unit

    # -----------------------------------------------------------------

    @property
    def is_per_angular_or_intrinsic_area(self):
        return self.unit.is_per_angular_or_intrinsic_area

    # -----------------------------------------------------------------

    @property
    def corresponding_angular_or_intrinsic_area_unit(self):
        return self.unit.corresponding_angular_or_intrinsic_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already
        if self.is_per_angular_or_intrinsic_area: return

        # Inform the user
        log.info("Converting the datacube to the corresponding angular or intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the unit
        unit = self.corresponding_angular_or_intrinsic_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert
        self.convert_by_factor(factor, unit)

        # Return the conversion factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already
        if self.is_per_angular_or_intrinsic_area: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding angular or intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the unit
        unit = self.corresponding_angular_or_intrinsic_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Return converted
        return self.converted_by_factor(factor, unit)

    # -----------------------------------------------------------------

    @property
    def corresponding_non_angular_or_intrinsic_area_unit(self):
        return self.unit.corresponding_non_angular_or_intrinsic_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_non_angular_or_intrinsic_area_unit(self):

        """
        Thisnfunction ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already
        if not self.is_per_angular_or_intrinsic_area: return

        # Inform the user
        log.info("Converting the datacube to the corresponding non- angular or intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the corresponding non angular or intrinsic area unit
        unit = self.corresponding_non_angular_or_intrinsic_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert
        self.convert_by_factor(factor, unit)

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_non_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        # Already
        if not self.is_per_angular_or_intrinsic_area: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding non- angular or intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the unit
        unit = self.corresponding_non_angular_or_intrinsic_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Return converted
        return self.converted_by_factor(factor, unit)

    # -----------------------------------------------------------------

    @property
    def is_per_angular_area(self):
        return self.unit.is_per_angular_area

    # -----------------------------------------------------------------

    def convert_to_corresponding_angular_area_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already per angular area
        if self.is_per_angular_area: return

        # Inform the user
        log.info("Converting the datacube to the corresponding angular area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the new unit
        unit = self.corresponding_angular_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert by the factor
        self.convert_by_factor(factor, unit)

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_angular_area_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already per ngular are
        if self.is_per_angular_area: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding angular area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the new unit
        unit = self.corresponding_angular_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Return converted
        return self.converted_by_factor(factor, unit)

    # -----------------------------------------------------------------

    @property
    def is_per_intrinsic_area(self):
        return self.unit.is_per_intrinsic_area

    # -----------------------------------------------------------------

    def convert_to_corresponding_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already per intrinsic area
        if self.is_per_intrinsic_area: return

        # Inform the user
        log.info("Converting the datacube to the corresponding intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the corresponding intrinsic area unit
        unit = self.corresponding_intrinsic_area_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert
        self.convert_by_factor(factor, unit)

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_intrinsic_area_unit(self):

        """
        Thins function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already per intrinsic area: return a copy
        if self.is_per_intrinsic_area: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding intrinsic area unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Create
        return self.converted_to(self.corresponding_intrinsic_area_unit)

    # -----------------------------------------------------------------

    def convert_to_corresponding_non_brightness_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already not a brightness
        if not self.is_brightness: return

        # Inform the user
        log.info("Converting the datacube to the corresponding non-brightness unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the new unit
        unit = self.corresponding_non_brightness_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Convert
        self.convert_by_factor(factor, unit)

        # Return the factor
        return factor

    # -----------------------------------------------------------------

    def converted_to_corresponding_non_brightness_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already not a brightness: return a copy
        if not self.is_brightness: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding non-brightness unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the new unit
        unit = self.corresponding_non_brightness_unit

        # Get the conversion factor
        factor = self.get_conversion_factor(unit)

        # Return converted
        return self.converted_by_factor(factor, unit)

    # -----------------------------------------------------------------

    @property
    def is_wavelength_density(self):
        return self.unit.is_wavelength_density

    # -----------------------------------------------------------------

    @property
    def corresponding_wavelength_density_unit(self):
        return self.unit.corresponding_wavelength_density_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_wavelength_density_unit(self, distance=None, silent=False, wavelength_unit=None):

        """
        This function ...
        :param distance:
        :param silent:
        :param wavelength_unit:
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already wavelength density
        if self.is_wavelength_density: return

        # Inform the user
        log.info("Converting the datacube to the corresponding wavelength density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the list of wavelengths
        wavelengths = self.wavelength_grid.wavelengths(unit="micron", add_unit=True)

        # Set the distance
        if distance is None: distance = self.distance

        # Keep list of conversion factors
        factors = []

        # Loop over the frames
        for i in range(self.nframes):

            # Get the wavelength
            wavelength = wavelengths[i]

            # Convert the frame
            factor = self.frames[i].convert_to_corresponding_wavelength_density_unit(wavelength=wavelength, distance=distance, wavelength_unit=wavelength_unit)
            factors.append(factor)

            # Debugging
            if not silent: log.debug("Conversion factor for frame " + str(i + 1) + " (wavelength = " + tostr(self.get_wavelength(i)) + "): " + str(factor))

        # Return the factors
        return factors

    # -----------------------------------------------------------------

    def converted_to_corresponding_wavelength_density_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already wavelength density
        if self.is_wavelength_density: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding wavelength density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Create
        return self.converted_to(self.corresponding_wavelength_density_unit)

    # -----------------------------------------------------------------

    @property
    def is_frequency_density(self):
        return self.unit.is_frequency_unit

    # -----------------------------------------------------------------

    @property
    def corresponding_frequency_density_unit(self):
        return self.unit.corresponding_frequency_density_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_frequency_density_unit(self, silent=False, frequency_unit=None):

        """
        This function ...
        :param silent:
        :param frequency_unit:
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already frequency density
        if self.is_frequency_density: return

        # Inform the user
        log.info("Converting the datacube to the corresponding frequency density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the list of wavelengths
        wavelengths = self.wavelength_grid.wavelengths(unit="micron", add_unit=True)

        # Keep list of factors
        factors = []

        # Loop over the frames
        for i in range(self.nframes):

            # Get the wavelength
            wavelength = wavelengths[i]

            # Convert the frame
            factor = self.frames[i].convert_to_corresponding_frequency_density_unit(wavelength=wavelength, frequency_unit=frequency_unit)
            factors.append(factor)

            # Debugging
            if not silent: log.debug("Conversion factor for frame " + str(i + 1) + " (wavelength = " + tostr(self.get_wavelength(i)) + "): " + str(factor))

        # Return the factors
        return factors

    # -----------------------------------------------------------------

    def converted_to_corresponding_frequency_density_unit(self):

        """
        This function ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already frequency density
        if self.is_frequency_density: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding frequency density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Create
        return self.converted_to(self.corresponding_frequency_density_unit)

    # -----------------------------------------------------------------

    @property
    def is_neutral_density(self):
        return self.unit.is_neutral_density

    # -----------------------------------------------------------------

    @property
    def corresponding_neutral_density_unit(self):
        return self.unit.corresponding_neutral_density_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_neutral_density_unit(self, silent=False):

        """
        Thisf unction ...
        :param silent:
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # ALready neutral density
        if self.is_neutral_density: return

        # Inform the user
        log.info("Converting the datacube to the corresponding neutral density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Get the list of wavelengths
        wavelengths = self.wavelength_grid.wavelengths(unit="micron", add_unit=True)

        # Keep list of conversion factors
        factors = []

        # Loop over the frames
        for i in range(self.nframes):

            # Get the wavelength
            wavelength = wavelengths[i]

            # Convert the frame
            factor = self.frames[i].convert_to_corresponding_neutral_density_unit(wavelength=wavelength)
            factors.append(factor)

            # Debugging
            if not silent: log.debug("Conversion factor for frame " + str(i + 1) + " (wavelength = " + tostr(self.get_wavelength(i)) + "): " + str(factor))

        # Return the factors
        return factors

    # -----------------------------------------------------------------

    def converted_to_corresponding_neutral_density_unit(self):

        """
        Thisf unction ...
        :return:
        """

        # Check whether unit is defined
        if not self.has_unit: raise ValueError("Unit of the datacube is not defined")

        # Already neutral density
        if self.is_neutral_density: return self.copy()

        # Inform the user
        log.info("Creating a datacube in the corresponding neutral density unit ...")

        # Show current unit
        log.debug("Current unit: " + tostr(self.unit, add_physical_type=True))

        # Create
        return self.converted_to(self.corresponding_neutral_density_unit)

    # -----------------------------------------------------------------

    @property
    def wavelength_array(self):
        return self.wavelengths(unit=self.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return self.wavelength_grid.unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit_string(self):
        return tostr(self.wavelength_unit)

    # -----------------------------------------------------------------

    @property
    def nans(self):

        """
        This function ...
        :return: 
        """

        # Initialize image
        masks = Image("nans")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame nans
            nans = self.frames[frame_name].nans
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(nans, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def largest_nans(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("largest_nans")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame nans
            nans = self.frames[frame_name].nans.largest()
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(nans, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def infs(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("infs")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame infs
            infs = self.frames[frame_name].infs
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(infs, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def largest_infs(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("largest_infs")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame infs
            infs = self.frames[frame_name].infs.largest()
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(infs, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def invalid(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("invalid")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame invalid mask
            invalid = self.frames[frame_name].invalid
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(invalid, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def valid(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("valid")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame valid mask
            valid = self.frames[frame_name].valid
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(valid, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def zeroes(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("zeroes")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame zeroes
            zeroes = self.frames[frame_name].zeroes
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(zeroes, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    @property
    def nonzeroes(self):

        """
        This function ...
        :return:
        """

        # Initialize image
        masks = Image("nonzeroes")

        # Loop over the frames
        for frame_name in self.frame_names:

            # Get frame nonzeroes
            nonzeroes = self.frames[frame_name].nonzeroes
            mask_name = frame_name.replace("frame", "mask")

            # Add mask
            masks.add_mask(nonzeroes, mask_name)

        # Return the masks
        return masks

    # -----------------------------------------------------------------

    def saveto(self, path, add_metadata=True, origin=None, add_masks=True, add_segments=True, add_regions=False):

        """
        This function ...
        :param path:
        :param add_metadata:
        :param origin:
        :param add_masks:
        :param add_segments:
        :param add_regions:
        :return:
        """

        # Set extra info
        extra_info = dict()
        extra_info["wavunit"] = self.wavelength_unit_string
        extra_info["wvlngths"] = self.wavelengths_string(add_unit=False)

        # Call the base class implementation
        super(DataCube, self).saveto(path, add_metadata=add_metadata, origin=origin, add_masks=add_masks,
                                     add_segments=add_segments, add_regions=add_regions, extra_info=extra_info,
                                     add_plane_names=False, add_filter=False)

# -----------------------------------------------------------------

def _do_one_filter_convolution_from_file(datacube_path, wavelengthgrid_path, result_path, unit, fltrname):

    """
    This function ...
    :param datacube_path:
    :param wavelengthgrid_path:
    :param result_path:
    :param unit:
    :param fltrname:
    :return:
    """

    message_prefix = "[convolution with " + fltrname + " filter] "

    # Inform the user
    log.info(message_prefix + "Loading filter ...")

    # Resurrect the filter
    fltr = BroadBandFilter(fltrname)

    # Inform the user
    log.info(message_prefix + "Loading wavelength grid ...")

    # Resurrect the wavelength grid
    wavelength_grid = WavelengthGrid.from_file(wavelengthgrid_path)

    # Inform the user
    log.info(message_prefix + "Loading datacube ...")

    # Resurrect the datacube
    datacube = DataCube.from_file(datacube_path, wavelength_grid)

    # Inform the user
    log.info(message_prefix + "Getting wavelength array ...")

    # Get the array of wavelengths
    wavelengths = datacube.wavelengths(asarray=True, unit="micron")

    # Inform the user
    log.info(message_prefix + "Converting datacube to 3D array ...")

    # Convert the datacube to a numpy array where wavelength is the third dimension (index=2)
    array = datacube.asarray(axis=2)

    # Check shape of the datacube
    if array.shape[2] != len(wavelengths): raise ValueError("The wavelength axis of the 3D array must be the last one")

    # Debugging
    log.debug("Image shape: " + str(array.shape))

    # Slice the cube array to only the required range of wavelengths for this filter
    use = (fltr.min.to("micron").value <= wavelengths) * (wavelengths <= fltr.max.to("micron").value)
    #print(use.shape)
    array = array[:, :, use]
    #array = array[use, :, :]
    wavelengths = wavelengths[use]

    # Inform the user
    log.info(message_prefix + "Starting convolution ...")

    # Do the convolution, time it
    with time.elapsed_timer() as elapsed:

        # Convolve
        data = fltr.convolve(wavelengths, array)

        # Show time
        log.success("Convolved the datacube with the " + str(fltr) + " filter in " + str(elapsed()) + " seconds")

    # Inform the user
    log.info(message_prefix + "Convolution completed")

    # Create frame, set properties
    frame = Frame(data, unit=unit, filter=fltr, wcs=datacube.wcs)

    # Inform
    log.info(message_prefix + "Saving result to " + result_path + " ...")

    # Save the frame with the index as name
    frame.saveto(result_path)

    # Success
    if fs.is_file(result_path): log.success(message_prefix + "Succesfully saved the convolved frame for the '" + fltrname + "' filter to '" + result_path + "'")
    else: raise RuntimeError("Something went wrong saving the resulting frame")

# -----------------------------------------------------------------

def _do_one_filter_convolution(fltr, wavelengths, array, unit, wcs):

    """
    This function ...
    :param fltr:
    :param wavelengths:
    :param array:
    :return:
    """

    # Debugging
    log.debug("Convolving the datacube with the " + str(fltr) + " filter ...")

    # Check shape
    if array.shape[2] != len(wavelengths): raise ValueError("The wavelength axis of the 3D array must be the last one")

    # Debugging
    log.debug("Image shape: " + str(array.shape[1]) + " x " + str(array.shape[0]))

    # Slice the cube array to just the wavaelength range required for the filter
    use = (fltr.min.to("micron").value <= wavelengths) * (wavelengths <= fltr.max.to("micron").value)
    #print(use.shape)
    array = array[:, :, use]
    #array = array[use, :, :]
    wavelengths = wavelengths[use]

    # Debugging
    log.debug("The number of wavelengths for this filter is " + str(len(wavelengths)))

    # Calculate the observed image frame, time it
    with time.elapsed_timer() as elapsed:

        # Do the convolution
        data, wavelength_grid = fltr.convolve(wavelengths, array, return_grid=True)

        # Show time
        log.success("Convolved the datacube with the " + str(fltr) + " filter in " + str(elapsed()) + " seconds")

    # Create and return the frame and the wavelength grid
    return Frame(data, unit=unit, filter=fltr, wcs=wcs), wavelength_grid

# -----------------------------------------------------------------

def needs_spectral_convolution(fltr, spectral_convolution):

    """
    This function ...
    :param fltr:
    :param spectral_convolution: flag or list of Filters
    :return:
    """

    # Broad band filter
    if isinstance(fltr, BroadBandFilter):

        # Single boolean
        if types.is_boolean_type(spectral_convolution): return spectral_convolution

        # Sequence of filters: return whether the filter is in it
        elif types.is_sequence(spectral_convolution): return fltr in spectral_convolution

        # Invalid
        else: raise ValueError("Invalid option for 'spectral_convolution'")

    # Narrow band filters: no spectral convolution
    else: return False

# -----------------------------------------------------------------
