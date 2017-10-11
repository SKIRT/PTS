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
from collections import defaultdict
import numpy as np

# Import the relevant PTS classes and modules
from .image import Image
from .frame import Frame
from ...core.data.sed import SED, ObservedSED
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

# -----------------------------------------------------------------

parallel_filter_convolution_dirname = "datacube-parallel-filter-convolution"

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

def load_skirt_datacubes(output_path=None, contribution="total"):

    """
    This function ...
    :param output_path:
    :param contribution:
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
        datacube = DataCube.from_file_and_sed_file(datacube_path, sed_path)

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

    def __init__(self, name="untitled"):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(DataCube, self).__init__(name)

        # The wavelength grid
        self.wavelength_grid = None

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_output(cls, instrument_name, output_path=None, contribution="total"):

        """
        This function ...
        :param instrument_name:
        :param output_path:
        :param contribution:
        :return:
        """

        # Load paths
        datacube_paths_instruments = load_all_skirt_datacube_paths(output_path=output_path)
        sed_paths_instruments = load_skirt_sed_paths(output_path=output_path)

        # Get the desired datacube path
        datacube_path = datacube_paths_instruments[instrument_name][contribution]

        # Get the SED path for the instrument
        sed_path = sed_paths_instruments[instrument_name]

        # Return
        return cls.from_file_and_sed_file(datacube_path, sed_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_sed_file(cls, image_path, sed_path):

        """
        This function ...
        :param image_path:
        :param sed_path:
        :return:
        """

        from ...core.data.sed import load_sed

        # Load the SED (can be ObservedSED, and can be from SKIRT output)
        sed = load_sed(sed_path)

        # Create
        return cls.from_file_and_sed(image_path, sed)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_sed(cls, image_path, sed):

        """
        This function ...
        :param image_path:
        :param sed:
        :return:
        """

        # Get the wavelength grid
        wavelength_grid = WavelengthGrid.from_sed(sed)

        # Create the datacube
        return cls.from_file(image_path, wavelength_grid)

    # -----------------------------------------------------------------

    @classmethod
    def from_file_and_wavelength_grid_file(cls, image_path, wavelengths_path):

        """
        This function ...
        :param image_path:
        :param wavelengths_path:
        :return:
        """

        # Get the wavelength grid
        wavelength_grid = WavelengthGrid.from_file(wavelengths_path)

        # Create the datacube
        return cls.from_file(image_path, wavelength_grid)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, image_path, wavelength_grid):

        """
        This function ...
        :param image_path:
        :param wavelength_grid:
        :return:
        """

        # Call the corresponding base class function
        datacube = super(DataCube, cls).from_file(image_path, always_call_first_primary=False, no_filter=True,
                                                  density=True, density_strict=True) # IMPORTANT: ASSUME THAT DATACUBES ARE ALWAYS DEFINED IN SPECTRAL DENSITY UNITS!

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
    def from_files(cls, paths):

        """
        This function ...
        :param paths: paths of frames
        :return:
        """

        frames = [Frame.from_file(path) for path in paths]
        return cls.from_frames(frames)

    # -----------------------------------------------------------------

    @classmethod
    def from_frames(cls, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Create a datacube instance
        datacube = cls()

        # The indices of the frames, sorted on wavelength
        sorted_indices = sorted(range(len(frames)), key=lambda i: frames[i].filter.pivotwavelength())

        # The list of wavelengths
        wavelengths = []

        # Add the frames
        nframes = 0
        for index in sorted_indices:

            # Add the frame
            frame_name = "frame" + str(nframes)
            datacube.add_frame(frames[index], frame_name)

            # Add the wavelength
            wavelengths.append(frames[index].filter.pivotwavelength())

            # Increment the number of frames
            nframes += 1

        # Create the wavelength grid
        datacube.wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit="micron")

        # Return the datacube
        return datacube

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.wavelengths(unit, asarray, add_unit)

    # -----------------------------------------------------------------

    def wavelength_indices(self, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        return self.wavelength_grid.wavelength_indices(min_wavelength, max_wavelength)

    # -----------------------------------------------------------------

    def asarray(self, axis=3):

        """
        This function ...
        :return:
        """

        # Get a list that contains the frames
        frame_list = self.frames.as_list()

        # Stack the frames into a 3D numpy array
        if axis == 3: return np.dstack(frame_list)
        elif axis == 2: return np.hstack(frame_list)
        elif axis == 1: return np.vstack(frame_list)
        elif axis == 0: return np.stack(frame_list)
        else: raise ValueError("'axis' parameter should be integer 0-3")

    # -----------------------------------------------------------------

    def get_frame_index_for_wavelength(self, wavelength):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.closest_wavelength_index(wavelength)

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

    def get_frame_for_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        index = self.get_frame_index_for_wavelength(wavelength)
        return self.frames[index]

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

        # Initialize the SED
        sed = ObservedSED(photometry_unit=self.unit)

        # Create a mask from the region (or shape)
        mask = region.to_mask(self.xsize, self.ysize)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixels that belong to the region
            flux = np.sum(self.frames[frame_name][mask]) * self.unit

            # Get error
            if errorcube is not None:

                # Get the error in the pixel
                error = errorcube.frames[frame_name][mask] * self.unit
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

        # Initialize the SED
        sed = ObservedSED(photometry_unit=self.unit)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixel
            flux = self.frames[frame_name][y, x] * self.unit

            # Get error
            if errorcube is not None:

                # Get the error in the pixel
                error = errorcube.frames[frame_name][y, x] * self.unit
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
        if isinstance(mask, basestring): inverse_mask = self.masks[mask].inverse()
        elif isinstance(mask, Mask): inverse_mask = mask.inverse()
        elif mask is None: inverse_mask = None
        else: raise ValueError("Mask must be string or Mask (or None) instead of " + str(type(mask)))

        # Initialize the SED
        sed = SED(photometry_unit=self.unit)

        # Loop over the wavelengths
        for index in self.wavelength_indices(min_wavelength, max_wavelength):

            # Get the wavelength
            wavelength = self.wavelength_grid[index]

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Calculate the total flux
            if mask is not None: total_flux = self.frames[frame_name].sum() * self.unit
            else: total_flux = np.sum(self.frames[frame_name][inverse_mask]) * self.unit

            # Add an entry to the SED
            sed.add_point(wavelength, total_flux)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def convolve_with_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with the " + str(fltr) + " filter ...")

        # Convert the datacube to a numpy array where wavelength is the third dimension
        array = self.asarray()

        # Calculate the observed image frame
        data = fltr.convolve(self.wavelengths(asarray=True), array)
        frame = Frame(data)

        # Set the filter
        frame.filter = fltr

        # Set the unit of the frame
        frame.unit = self.unit

        # Set the wcs of the frame
        frame.wcs = self.wcs

        # Return the resulting frame
        return frame

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
            index = self.get_frame_index_for_wavelength(fltr.pivot)

            # Get a copy of the frame
            frame = self.frames[index].copy()

        # Set filter to be sure
        frame.filter = fltr

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def frames_for_filters(self, filters, convolve=False, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param convolve:
        :param nprocesses:
        :param check_previous_sessions:
        :return:
        """

        # Inform the user
        log.info("Getting frames for " + str(len(filters)) + " different filters ...")

        # Limit the number of processes to the number of filters
        nprocesses = min(nprocesses, len(filters))

        frames = []
        for_convolution = []

        # Loop over the filters
        for fltr in filters:

            # Needs spectral convolution?
            if needs_spectral_convolution(fltr, convolve):

                # Debugging
                log.debug("The frame for the " + str(fltr) + " filter will be calculated by convolving spectrally")

                # Add to list
                for_convolution.append(fltr)

                # Add placeholder
                frames.append(None)

            # No spectral convolution
            else:

                # Debugging
                log.debug("Getting the frame for the " + str(fltr) + " filter ...")

                # Get the index of the wavelength closest to that of the filter
                index = self.get_frame_index_for_wavelength(fltr.pivot)

                # Make a copy of the frame
                frame = self.frames[index].copy()

                # Set the filter
                frame.filter = fltr

                # Get the frame
                frames.append(frame)

        # Calculate convolved frames
        if len(for_convolution) > 0:
            # Debugging
            log.debug(str(len(for_convolution)) + " filters require spectral convolution")
            convolved_frames = self.convolve_with_filters(for_convolution, nprocesses=nprocesses, check_previous_sessions=check_previous_sessions)
        else:
            # Debugging
            log.debug("Spectral convolution will be used for none of the filters")
            convolved_frames = []

        # Add the convolved frames
        for fltr, frame in zip(for_convolution, convolved_frames):

            # Set filter to be sure
            frame.filter = fltr

            # Set the frame
            index = filters.index(fltr)
            frames[index] = frame

        # Return the list of frames
        return frames

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with " + str(len(filters)) + " different filters ...")

        # Limit the number of processes to the number of filters
        nprocesses = min(nprocesses, len(filters))

        # PARALLEL EXECUTION
        if nprocesses > 1: return self.convolve_with_filters_parallel(filters, nprocesses=nprocesses, check_previous_sessions=check_previous_sessions)

        # SERIAL EXECUTION
        else: return self.convolve_with_filters_serial(filters)

    # -----------------------------------------------------------------

    def convolve_with_filters_serial(self, filters):

        """
        Thisj function ...
        :param filters:
        :return:
        """

        # Debugging
        log.debug("Convolving the datacube with " + str(len(filters)) + " different filters on one process ...")

        # Initialize list to contain the output frames per filter
        nfilters = len(filters)
        frames = [None] * nfilters

        # Debugging
        log.debug("Converting the datacube into a single 3D array ...")

        # Convert the datacube to a numpy array where wavelength is the third dimension
        array = self.asarray()

        # Get the array of wavelengths
        wavelengths = self.wavelengths(asarray=True, unit="micron")

        # Loop over the filters
        for index in range(nfilters):

            # Get the current filter
            fltr = filters[index]

            # Do the filter convolution, put frame in the frames list
            _do_one_filter_convolution(fltr, wavelengths, array, frames, index, self.unit, self.wcs)

        # Return the list of resulting frames
        return frames

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

    def convolve_with_filters_parallel(self, filters, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
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
        return frames

    # -----------------------------------------------------------------

    def to_wavelength_density(self, new_unit, wavelength_unit):

        """
        This function ...
        :param new_unit:
        :param wavelength_unit:
        :return:
        """

        # Inform the user
        log.info("Converting the datacube from neutral flux density to flux density per unit of wavelength (in " + wavelength_unit + ")")

        # Get list of wavelengths in desired unit
        wavelengths = self.wavelength_grid.wavelengths(unit=wavelength_unit, add_unit=False)

        # Convert the frames from neutral surface brightness to wavelength surface brightness
        for l in range(self.nframes):

            # Get the wavelength
            wavelength = wavelengths[l]

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(l)

            # Divide this frame by the wavelength in micron
            self.frames[frame_name] /= wavelength

        # Set the new unit of the datacube
        self.unit = new_unit

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

    # Convert the datacube to a numpy array where wavelength is the third dimension
    array = datacube.asarray()

    # Inform the user
    log.info(message_prefix + "Starting convolution ...")

    # Do the convolution
    data = fltr.convolve(wavelengths, array)

    # Inform the user
    log.info(message_prefix + "Convolution completed")

    # Create frame, set properties
    frame = Frame(data)
    frame.unit = unit
    frame.filter = fltr
    frame.wcs = datacube.wcs

    # Inform
    log.info(message_prefix + "Saving result to " + result_path + " ...")

    # Save the frame with the index as name
    frame.saveto(result_path)

    # Success
    if fs.is_file(result_path): log.success(message_prefix + "Succesfully saved the convolved frame for the '" + fltrname + "' filter to '" + result_path + "'")
    else: raise RuntimeError("Something went wrong saving the resulting frame")

# -----------------------------------------------------------------

def _do_one_filter_convolution(fltr, wavelengths, array, frames, index, unit, wcs):

    """
    This function ...
    :param fltr:
    :param wavelengths:
    :param array:
    :param frames:
    :param index:
    :return:
    """

    # Debugging
    log.debug("Convolving the datacube with the " + str(fltr) + " filter ...")

    # Calculate the observed image frame
    data = fltr.convolve(wavelengths, array)
    frame = Frame(data)

    # Debugging
    log.success("Convolved the datacube with the " + str(fltr) + " filter ...")

    # Set the unit of the frame
    frame.unit = unit

    # Set the filter
    frame.filter = fltr

    # Set the wcs
    frame.wcs = wcs

    # Add the frame to the list
    frames[index] = frame

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
