#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.output Contains the SimulationOutput, ExtractionOutput, PlottingOutput and MiscOutput classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..basics.map import Map
from ..tools.utils import lazyproperty, memoize_method
from ..units.unit import parse_unit as u
from ..basics.log import log

# -----------------------------------------------------------------

cached_filename = "cached.dat"
cached_type_name = "cached"

# -----------------------------------------------------------------

def get_cache_path(directory_path):

    """
    This function ...
    :param directory_path:
    :return:
    """

    # Determine file path
    filepath = fs.join(directory_path, cached_filename)
    if not fs.is_file(filepath): return None

    # Get the cache path
    return fs.get_first_line(filepath)

# -----------------------------------------------------------------

def write_cache_path(directory_path, cache_path):

    """
    This function ...
    :param directory_path:
    :param cache_path:
    :return:
    """

    # Determine file path
    filepath = fs.join(directory_path, cached_filename)

    # Write
    fs.write_line(filepath, cache_path)

# -----------------------------------------------------------------

def get_all_output_cwd(**kwargs):

    """
    This function ...
    :return:
    """

    return get_all_output(fs.cwd(), **kwargs)

# -----------------------------------------------------------------

def get_all_output(path, out_name="out", extr_name="extr", plot_name="plot", misc_name="misc"):

    """
    This function ...
    :param out_name:
    :param extr_name:
    :param plot_name:
    :param misc_name:
    :return:
    """

    # Get directory names
    dirnames = fs.directories_in_path(path, returns="name")

    # There is an output directory
    if out_name in dirnames: output = get_output(fs.join(path, out_name))
    else: output = get_output(path)

    # Get extraction output
    if extr_name in dirnames: extraction = get_extraction(fs.join(path, extr_name))
    else: extraction = get_extraction(path)

    # Get plotting output
    if plot_name in dirnames: plotting = get_plotting(fs.join(path, plot_name))
    else: plotting = get_plotting(path)

    # Get misc output
    if misc_name in dirnames: misc = get_misc(fs.join(path, misc_name))
    else: misc = get_misc(path)

    # Return
    return output, extraction, plotting, misc

# -----------------------------------------------------------------

def get_output_cwd():

    """
    This function ...
    :return:
    """

    return SimulationOutput.from_cwd()

# -----------------------------------------------------------------

def get_output(path):

    """
    This function ...
    :param path:
    :return:
    """

    return SimulationOutput.from_directory(path)

# -----------------------------------------------------------------

def get_extraction_output_cwd():

    """
    This function ...
    :return:
    """

    return ExtractionOutput.from_cwd()

# -----------------------------------------------------------------

def get_extraction(path):

    """
    This function ...
    :param path:
    :return:
    """

    return ExtractionOutput.from_directory(path)

# -----------------------------------------------------------------

def get_plotting_cwd():

    """
    This function ...
    :return:
    """

    return PlottingOutput.from_cwd()

# -----------------------------------------------------------------

def get_plotting(path):

    """
    This function ...
    :param path:
    :return:
    """

    return PlottingOutput.from_directory(path)

# -----------------------------------------------------------------

def get_misc_cwd():

    """
    This function ...
    :return:
    """

    return MiscOutput.from_cwd()

# -----------------------------------------------------------------

def get_misc(path):

    """
    This function ...
    :param path:
    :return:
    """

    return MiscOutput.from_directory(path)

# -----------------------------------------------------------------

# The various output types
output_types = Map()
output_types.isrf = "isrf"
output_types.absorption = "abs"
output_types.temperature = "temp"
output_types.seds = "sed"
output_types.images = "image"
output_types.total_images = "image-total"
output_types.count_images = "image-counts"
output_types.direct_images = "image-direct"
output_types.transparent_images = "image-transparent"
output_types.scattered_images = "image-scattered"
output_types.dust_images = "image-dust"
output_types.dust_scattered_images = "image-dustscattered"
output_types.cell_temperature = "celltemp"
output_types.logfiles = "log"
output_types.wavelengths = "wavelengths"
output_types.grid = "grid"
output_types.gdensity = "grho"
output_types.tdensity = "trho"
output_types.cell_properties = "cellprops"
output_types.stellar_density = "stellar_density"
output_types.tree = "tree"
output_types.convergence = "convergence"
output_types.dust_mass = "mass"
output_types.dust_mix_properties = "mean"
output_types.dust_optical_properties = "optical"
output_types.dust_grain_sizes = "size"
output_types.parameters = "parameters"

# -----------------------------------------------------------------

# Description of the output types
output_type_choices = dict()
output_type_choices[output_types.isrf] = "interstellar radiation field strength"
output_type_choices[output_types.absorption] = "absorption luminosities"
output_type_choices[output_types.temperature] = "temperature"
output_type_choices[output_types.seds] = "all SEDs"
output_type_choices[output_types.images] = "all datacubes"
output_type_choices[output_types.total_images] = "datacubes of total emission"
output_type_choices[output_types.count_images] = "datacubes of photon counts for total emission"
output_type_choices[output_types.direct_images] = "datacubes of direct emission"
output_type_choices[output_types.transparent_images] = "datacubes of transparent emission"
output_type_choices[output_types.scattered_images] = "datacubes of scattered emission"
output_type_choices[output_types.dust_images] = "datacubes of dust emission"
output_type_choices[output_types.dust_scattered_images] = "datacubes of scattered dust emission"
output_type_choices[output_types.cell_temperature] = "temperature per dust cell"
output_type_choices[output_types.logfiles] = "log files"
output_type_choices[output_types.wavelengths] = "wavelength files"
output_type_choices[output_types.grid] = "grid files"
output_type_choices[output_types.gdensity] = "grid dust density"
output_type_choices[output_types.tdensity] = "theoretical dust density"
output_type_choices[output_types.cell_properties] = "dust cell properties"
output_type_choices[output_types.stellar_density] = "stellar component density in the dust grid"
output_type_choices[output_types.tree] = "dust grid tree data file"
output_type_choices[output_types.convergence] = "convergence file"
output_type_choices[output_types.dust_mass] = "dust population masses"
output_type_choices[output_types.dust_mix_properties] = "combined dust mix properties"
output_type_choices[output_types.dust_optical_properties] = "optical dust population properties"
output_type_choices[output_types.dust_grain_sizes] = "dust grain size information"
output_type_choices[output_types.parameters] = "parameter files"

# -----------------------------------------------------------------

def get_output_type(filename):

    """
    This function ...
    :param filename:
    :return:
    """

    # Cached file
    if filename == cached_filename: return cached_type_name

    ## ISRF
    elif filename.endswith("_ds_isrf.dat"): return output_types.isrf

    ## Absorption
    elif filename.endswith("_ds_abs.dat"): return output_types.absorption

    ## Temperature
    elif "_ds_temp" in filename and filename.endswith(".fits"): return output_types.temperature

    ## SED
    elif filename.endswith("_sed.dat"): return output_types.seds

    ## Total datacubes
    elif filename.endswith("_total.fits"): return output_types.total_images

    ## Counts datacubes
    elif filename.endswith("_counts.fits"): return output_types.count_images

    ## Direct datacubes
    elif filename.endswith("_direct.fits"): return output_types.direct_images

    ## Transparent datacubes
    elif filename.endswith("_transparent.fits"): return output_types.transparent_images

    ## Scattered datacubes
    elif filename.endswith("_scattered.fits"): return output_types.scattered_images

    ## Dust datacubes
    elif filename.endswith("_dust.fits"): return output_types.dust_images

    ## Dust scattered datacubes
    elif filename.endswith("_dustscattered.fits"): return output_types.dust_scattered_images

    ## Cell temperature data
    elif filename.endswith("_ds_celltemps.dat"): return output_types.cell_temperature

    ## Log files
    elif "_log" in filename and filename.endswith(".txt"): return output_types.logfiles

    ## Wavelength files
    elif filename.endswith("_wavelengths.dat"): return output_types.wavelengths

    ## Grid structure data
    elif "_ds_grid" in filename and filename.endswith(".dat"): return output_types.grid

    ## Grid dust density
    elif "_ds_grho" in filename and filename.endswith(".fits"): return output_types.gdensity

    ## Theoretical dust density
    elif "_ds_trho" in filename and filename.endswith(".fits"): return output_types.tdensity

    ## Dust cell properties
    elif "_ds_cellprops" in filename and filename.endswith(".dat"): return output_types.cell_properties

    ## Stellar density
    elif "_ds_stellar" in filename and filename.endswith(".dat"): return output_types.stellar_density

    ## Dust grid onvergence
    elif filename.endswith("_ds_convergence.dat"): return output_types.convergence

    ## Dust grid tree data
    elif "_ds_tree" in filename and filename.endswith(".dat"): return output_types.tree

    ## Dust mass
    elif "ds_mix_" in filename and "_mass" in filename: return output_types.dust_mass

    ## Dust mix properties
    elif "ds_mix_" in filename and "_mean" in filename: return output_types.dust_mix_properties

    ## Dust optical properties
    elif "ds_mix_" in filename and "_opti" in filename: return output_types.dust_optical_properties

    ## Dust grain sizes
    elif "ds_mix_" in filename and "_size" in filename: return output_types.dust_grain_sizes

    ## Parameter files
    elif filename.endswith("_parameters.xml") or filename.endswith("_parameters.tex"): return output_types.parameters

    # No match
    else: return None

# -----------------------------------------------------------------

parent_output_types = dict()
parent_output_types[output_types.total_images] = output_types.images
parent_output_types[output_types.direct_images] = output_types.images
parent_output_types[output_types.transparent_images] = output_types.images
parent_output_types[output_types.scattered_images] = output_types.images
parent_output_types[output_types.dust_images] = output_types.images
parent_output_types[output_types.dust_scattered_images] = output_types.images

# -----------------------------------------------------------------

def get_parent_type(output_type):

    """
    This function ...
    :param output_type:
    :return:
    """

    if output_type in parent_output_types: return parent_output_types[output_type]
    else: return output_type

# -----------------------------------------------------------------

other_name = "other"

# -----------------------------------------------------------------

class Output(object):

    """
    This function ...
    """

    __metaclass__ = ABCMeta
    _output_types = None
    _output_type_choices = None
    _with_directories = False
    _with_extensions = None
    _without_extensions = None
    _recursive = False

    # -----------------------------------------------------------------

    @staticmethod
    def get_type(filename, directory=None):

        """
        This function ...
        :param filename:
        :param directory:
        :return:
        """

        raise RuntimeError("This method should be implemented in the derived class")

    # -----------------------------------------------------------------

    @staticmethod
    def get_directory_type(dirname, directory=None):

        """
        This function ...
        :param dirname:
        :param directory:
        :return:
        """

        raise RuntimeError("This method should be implemented in the derived class")

    # -----------------------------------------------------------------

    def load_files(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Loading file paths ...")

        # Get flag
        get_prefix = kwargs.pop("get_prefix", False)
        prefix = None

        # Path of cached directory
        cached_path = None

        # Loop over the filepaths, categorize
        for filepath in args:

            # Get filename and directory path
            filename = fs.name(filepath)
            directory = fs.directory_of(filepath)

            # Get prefix
            if get_prefix:
                file_prefix = filename.split("_")[0]

                # Set and check prefix
                if prefix is None: prefix = file_prefix
                elif prefix != file_prefix: raise ValueError("Cannot add files with different simulation prefixes")

            # Get the output type
            output_type = self.get_type(filename, directory=directory)

            # Cached?
            if output_type == cached_type_name:
                if cached_path is not None: raise RuntimeError("Something went wrong")
                cached_path = get_cache_path(directory)
                continue

            # Add the type
            if output_type is None: self.files[other_name].append(filepath)
            else: self.files[output_type].append(filepath)

        # Load cached files and directories
        if cached_path is not None:

            # Files
            filepaths = fs.files_in_path(cached_path, startswith=prefix, extension=self._with_extensions,
                                         not_extension=self._without_extensions, recursive=self._recursive)
            self.load_files(*filepaths, get_prefix=get_prefix)

            # Directories
            if self._with_directories:
                dirpaths = fs.directories_in_path(cached_path, recursive=self._recursive)
                self.load_directories(*dirpaths)

        # Return the prefix
        if get_prefix: return prefix

    # -----------------------------------------------------------------

    def load_directories(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        # Debugging
        log.debug("Loading directory paths ...")

        # Loop over the directory paths, categorize
        for dirpath in args:

            # Get directory name and directory path
            dirname = fs.name(dirpath)
            directory = fs.directory_of(dirpath)

            # Get the output type
            output_type = self.get_directory_type(dirname, directory=directory)

            # Add the type
            if output_type is None: self.directories[other_name].append(dirpath)
            else: self.directories[output_type].append(dirpath)

    # -----------------------------------------------------------------

    def __init__(self, *args):

        """
        The constructor ...
        :param args:
        """

        # Dictionary of filepaths
        self.files = defaultdict(list)

        # Dictionary of directory paths
        self.directories = defaultdict(list)

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, *args, **kwargs):

        """
        This function ...
        :return:
        """

        raise RuntimeError("This method should be implemented in the derived class")

    # -----------------------------------------------------------------

    @classmethod
    def from_cwd(cls, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        return cls.from_directory(fs.cwd(), **kwargs)

    # -----------------------------------------------------------------

    @memoize_method
    def get_nfiles(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return len(self.files[output_type])

    # -----------------------------------------------------------------

    @memoize_method
    def get_ndirectories(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return len(self.directories[output_type])

    # -----------------------------------------------------------------

    @memoize_method
    def has_files(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return output_type in self.files and self.get_nfiles(output_type) > 0

    # -----------------------------------------------------------------

    @memoize_method
    def has_directories(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return output_type in self.directories and self.get_ndirectories(output_type) > 0

    # -----------------------------------------------------------------

    @memoize_method
    def has_single_file(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return output_type in self.files and self.get_nfiles(output_type) == 1

    # -----------------------------------------------------------------

    @memoize_method
    def has_single_directory(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        return output_type in self.directories and self.get_ndirectories(output_type) == 1

    # -----------------------------------------------------------------

    @memoize_method
    def get_files(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        if not self.has_files(output_type): return []
        return self.files[output_type]

    # -----------------------------------------------------------------

    @memoize_method
    def get_directories(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        if not self.has_directories(output_type): return []
        return self.directories[output_type]

    # -----------------------------------------------------------------

    @memoize_method
    def get_single_file(self, output_type):

        """
        This function ...
        :param output_type:
        :return:
        """

        if not self.has_single_file(output_type): raise IOError("Doesn't have a single '" + output_type + "' file")
        else: return self.get_files(output_type)[0]

    # -----------------------------------------------------------------

    @memoize_method
    def get_single_directory(self, output_type):

        """
        This function ...
        :return:
        """

        if not self.has_single_directory(output_type): raise IOError("Doesn't have a single '" + output_type + "' directory")
        else: return self.get_directories(output_type)[0]

    # -----------------------------------------------------------------

    @property
    def disk_size(self):

        """
        This function ...
        :return:
        """

        # Initialize the total size
        size = 0. * u("byte")

        # Loop over all the directories not inside another directory
        for dirpath in self.all_directory_paths_not_in_directory: size += fs.directory_size(dirpath)

        # Loop over all the files not inside a directory (loose files)
        for filepath in self.loose_file_paths: size += fs.file_size(filepath)

        # Return the total size
        return size

    # -----------------------------------------------------------------

    @lazyproperty
    def directory(self):

        """
        This property returns the common directory that all files and directories are in, if applicable
        :return:
        """

        path = None

        # Check files
        for output_type in self.files:
            for filepath in self.files[output_type]:
                if path is None: path = fs.directory_of(filepath)
                elif path != fs.directory_of(filepath): return None #raise ValueError("Output files are not in the same directory")

        # Check directories
        for output_type in self.directories:
            for dirpath in self.directories[output_type]:
                if path is None: path = fs.directory_of(dirpath)
                elif path != fs.directory_of(dirpath): return None

        # Return the directory path
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def root_directory(self):

        """
        This property returns the most specific directory path that contains all of the files and directories
        :return:
        """

        return fs.common_directory(self.all_paths)

    # -----------------------------------------------------------------

    def relative_path(self, path):

        """
        Thisf unction ...
        :param path:
        :return:
        """

        return fs.relative_to(path, self.root_directory)

    # -----------------------------------------------------------------

    @lazyproperty
    def nother_files(self):

        """
        This function ...
        :return:
        """

        return len(self.files[other_name])

    # -----------------------------------------------------------------

    @lazyproperty
    def nother_directories(self):

        """
        This function ...
        :return:
        """

        return len(self.directories[other_name])

    # -----------------------------------------------------------------

    @lazyproperty
    def has_other_files(self):

        """
        This function ...
        :return:
        """

        return other_name in self.files and self.nother_files > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def has_other_directories(self):

        """
        This function ...
        :return:
        """

        return other_name in self.directories and self.nother_directories > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def other_files(self):

        """
        This function ...
        :return:
        """

        if not self.has_other_files: return []
        else: return self.files[other_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def other_directories(self):

        """
        This function ...
        :return:
        """

        if not self.has_other_directories: return []
        else: return self.directories[other_name]

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        # Files
        for output_type in self.files:
            for path in self.files[output_type]: yield path

        # Directories
        for output_type in self.directories:
            for path in self.directories[output_type]: yield path

    # -----------------------------------------------------------------

    @lazyproperty
    def all_paths(self):

        """
        This function ...
        :return:
        """

        return list(self.__iter__())

    # -----------------------------------------------------------------

    def iter_files(self):

        """
        This function ...
        :return:
        """

        for output_type in self.files:
            for path in self.files[output_type]: yield path

    # -----------------------------------------------------------------

    def iter_files_not_in_directory(self):

        """
        This function ...
        :return:
        """

        # Loop over all file paths
        for path in self.iter_files():

            # Skip files that are contained by any of the directories
            if fs.any_contains_path(self.all_directory_paths, path): continue

            # Give the file path
            yield path

    # -----------------------------------------------------------------

    def iter_directories(self):

        """
        This function ...
        :return:
        """

        for output_type in self.directories:
            for path in self.directories[output_type]: yield path

    # -----------------------------------------------------------------

    def iter_directories_not_with_file(self):

        """
        This function ...
        :return:
        """

        # Loop over all directory paths
        for path in self.iter_directories():

            # Skip directories that contain any of the files
            if fs.contains_any_path(path, self.all_file_paths): continue

            # Give the directory path
            yield path

    # -----------------------------------------------------------------

    def iter_directories_not_in_directory(self):

        """
        This function ...
        :return:
        """

        # Loop over all directory paths
        for path in self.iter_directories():

            # Skip directories that are contained by any of the (other) directories
            if fs.any_contains_path(self.all_directory_paths, path): continue

            # Give the directory path
            yield path

    # -----------------------------------------------------------------

    @lazyproperty
    def all_file_paths(self):

        """
        This function ...
        :return:
        """

        return list(self.iter_files())

    # -----------------------------------------------------------------

    @lazyproperty
    def all_file_paths_not_in_directory(self):

        """
        This function ...
        :return:
        """

        return list(self.iter_files_not_in_directory())

    # -----------------------------------------------------------------

    @property
    def loose_file_paths(self):

        """
        This function ...
        :return:
        """

        return self.all_file_paths_not_in_directory

    # -----------------------------------------------------------------

    @lazyproperty
    def all_directory_paths(self):

        """
        This function ...
        :return:
        """

        return list(self.iter_directories())

    # -----------------------------------------------------------------

    @lazyproperty
    def all_directory_paths_not_with_file(self):

        """
        This function ...
        :return:
        """

        return list(self.iter_directories_not_with_file())

    # -----------------------------------------------------------------

    @lazyproperty
    def all_directory_paths_not_in_directory(self):

        """
        This function ...
        :return:
        """

        return list(self.iter_directories_not_in_directory())

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        #return self.nfiles + self.ndirectories
        return self.nfiles

    # -----------------------------------------------------------------

    @property
    def nfiles(self):

        """
        This function ...
        :return:
        """

        total = 0

        # Files
        for output_type in self._output_types:
            if self.has_files(output_type): total += self.get_nfiles(output_type)
        if self.has_other_files: total += self.nother_files

        return total

    # -----------------------------------------------------------------

    @property
    def ndirectories(self):

        """
        This function ...
        :return:
        """

        total = 0

        # Directories
        for output_type in self._output_types:
            if self.has_directories(output_type): total += self.get_ndirectories(output_type)
        if self.has_other_directories: total += self.nother_directories

        return total

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.to_string()

    # -----------------------------------------------------------------

    def to_string(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        from ..tools import formatting as fmt

        lines = []

        # Loop over the output types
        for output_type in self._output_types:
            if not self.has_files(output_type): continue

            # Empty line
            lines.append(line_prefix)

            # Get number of files
            nfiles = self.get_nfiles(output_type)

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_type].capitalize() + fmt.reset + " (" + str(nfiles) + "):"
            lines.append(line_prefix + title)
            lines.append(line_prefix)

            # Add paths
            for path in self.files[output_type]: lines.append(line_prefix + " - " + self.relative_path(path))

        # Other
        if self.has_other_files:

            # Empty line
            lines.append(line_prefix)

            # Get number of files
            nfiles = self.nother_files

            # Add title
            title = fmt.green + fmt.underlined + "Other output" + fmt.reset + " (" + str(nfiles) + "):"
            lines.append(line_prefix + title)
            lines.append(line_prefix)

            # Add paths
            for path in self.other_files: lines.append(line_prefix + " - " + self.relative_path(path))

        # Add new line
        lines.append(line_prefix)

        # Return
        return "\n".join(lines)

    # -----------------------------------------------------------------

    def show(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        print(self.to_string(line_prefix=line_prefix))

# -----------------------------------------------------------------

class SimulationOutput(Output):

    """
    This class ...
    """

    _output_types = output_types
    _output_type_choices = output_type_choices
    _without_extensions = "ski"

    # -----------------------------------------------------------------

    def __init__(self, *args):

        """
        This function ...
        :param args:
        """

        # Call the constructor of the base class
        super(SimulationOutput, self).__init__()

        # Load the file paths, set the simulation prefix
        self.prefix = self.load_files(*args, get_prefix=True)

    # -----------------------------------------------------------------

    @staticmethod
    def get_type(filename, directory=None):

        """
        This function ...
        :param filename:
        :param directory:
        :return:
        """

        return get_output_type(filename)

    # -----------------------------------------------------------------

    @staticmethod
    def get_directory_type(dirname, directory=None):

        """
        This function ...
        :param dirname:
        :param directory:
        :return:
        """

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, prefix=None):

        """
        This function ...
        :param path:
        :param prefix:
        :return:
        """

        filepaths = fs.files_in_path(path, startswith=prefix, not_extension=cls._without_extensions)
        return cls(*filepaths)

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_directory(cls, path, remote, prefix=None):

        """
        This function ...
        :param path:
        :param remote:
        :param prefix:
        :return:
        """

        filepaths = remote.files_in_path(path, startswith=prefix, not_extension=cls._without_extensions)
        return cls(*filepaths)

    # -----------------------------------------------------------------

    @property
    def nisrf(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.isrf)

    # -----------------------------------------------------------------

    @property
    def has_isrf(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.isrf)

    # -----------------------------------------------------------------

    @property
    def has_single_isrf(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.isrf)

    # -----------------------------------------------------------------

    @property
    def isrf(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.isrf)

    # -----------------------------------------------------------------

    @property
    def single_isrf(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.isrf)

    # -----------------------------------------------------------------

    @property
    def nabsorption(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.absorption)

    # -----------------------------------------------------------------

    @property
    def has_absorption(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.absorption)

    # -----------------------------------------------------------------

    @property
    def has_single_absorption(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.absorption)

    # -----------------------------------------------------------------

    @property
    def absorption(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.absorption)

    # -----------------------------------------------------------------

    @property
    def single_absorption(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.absorption)

    # -----------------------------------------------------------------

    @property
    def ntemperature(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.temperature)

    # -----------------------------------------------------------------

    @property
    def has_temperature(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.temperature)

    # -----------------------------------------------------------------

    @property
    def temperature(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.temperature)

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def has_single_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def seds(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def single_sed(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.seds)

    # -----------------------------------------------------------------

    @lazyproperty
    def nimages(self):

        """
        This function ...
        :return:
        """

        total = 0
        if self.has_total_images: total += self.ntotal_images
        if self.has_count_images: total += self.ncount_images
        if self.has_direct_images: total += self.ndirect_images
        if self.has_transparent_images: total += self.ntransparent_images
        if self.has_scattered_images: total += self.nscattered_images
        if self.has_dust_images: total += self.ndust_images
        if self.has_dust_scattered_images: total += self.ndust_scattered_images
        return total

    # -----------------------------------------------------------------

    @lazyproperty
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.nimages > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def images(self):

        """
        This function ...
        :return:
        """

        paths = []
        if self.has_total_images: paths.extend(self.total_images)
        if self.has_count_images: paths.extend(self.count_images)
        if self.has_direct_images: paths.extend(self.direct_images)
        if self.has_transparent_images: paths.extend(self.transparent_images)
        if self.has_scattered_images: paths.extend(self.scattered_images)
        if self.has_dust_images: paths.extend(self.dust_images)
        if self.has_dust_scattered_images: paths.extend(self.dust_scattered_images)
        return paths

    # -----------------------------------------------------------------

    @property
    def ntotal_images(self):

        """
        Thisf unctin ...
        :return:
        """

        return self.get_nfiles(self._output_types.total_images)

    # -----------------------------------------------------------------

    @property
    def has_total_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.total_images)

    # -----------------------------------------------------------------

    @property
    def total_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.total_images)

    # -----------------------------------------------------------------

    @property
    def ncount_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.count_images)

    # -----------------------------------------------------------------

    @property
    def has_count_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.count_images)

    # -----------------------------------------------------------------

    @property
    def count_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.count_images)

    # -----------------------------------------------------------------

    @property
    def ndirect_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.direct_images)

    # -----------------------------------------------------------------

    @property
    def has_direct_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.direct_images)

    # -----------------------------------------------------------------

    @property
    def direct_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.direct_images)

    # -----------------------------------------------------------------

    @property
    def ntransparent_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.transparent_images)

    # -----------------------------------------------------------------

    @property
    def has_transparent_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.transparent_images)

    # -----------------------------------------------------------------

    @property
    def transparent_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.transparent_images)

    # -----------------------------------------------------------------

    @property
    def nscattered_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.scattered_images)

    # -----------------------------------------------------------------

    @property
    def has_scattered_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.scattered_images)

    # -----------------------------------------------------------------

    @property
    def scattered_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.scattered_images)

    # -----------------------------------------------------------------

    @property
    def ndust_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_images)

    # -----------------------------------------------------------------

    @property
    def has_dust_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_images)

    # -----------------------------------------------------------------

    @property
    def dust_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_images)

    # -----------------------------------------------------------------

    @property
    def ndust_scattered_images(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_scattered_images)

    # -----------------------------------------------------------------

    @property
    def has_dust_scattered_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_scattered_images)

    # -----------------------------------------------------------------

    @property
    def dust_scattered_images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_scattered_images)

    # -----------------------------------------------------------------

    @property
    def ncell_temperature(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.cell_temperature)

    # -----------------------------------------------------------------

    @property
    def has_cell_temperature(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.cell_temperature)

    # -----------------------------------------------------------------

    @property
    def cell_temperature(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.cell_temperature)

    # -----------------------------------------------------------------

    @property
    def nlogfiles(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.logfiles)

    # -----------------------------------------------------------------

    @property
    def has_logfiles(self):

        """
        This funtion ...
        :return:
        """

        return self.has_files(self._output_types.logfiles)

    # -----------------------------------------------------------------

    @property
    def logfiles(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.logfiles)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_single_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.wavelengths)

    # -----------------------------------------------------------------

    @property
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.wavelengths)

    # -----------------------------------------------------------------

    @property
    def single_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.wavelengths)

    # -----------------------------------------------------------------

    @property
    def ngrid(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.grid)

    # -----------------------------------------------------------------

    @property
    def has_grid(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_files(self._output_types.grid)

    # -----------------------------------------------------------------

    @property
    def grid(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.grid)

    # -----------------------------------------------------------------

    @property
    def ngdensity(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.gdensity)

    # -----------------------------------------------------------------

    @property
    def has_gdensity(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.gdensity)

    # -----------------------------------------------------------------

    @property
    def gdensity(self):

        """
        This fucntion ...
        :return:
        """

        return self.get_files(self._output_types.gdensity)

    # -----------------------------------------------------------------

    @property
    def ntdensity(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.tdensity)

    # -----------------------------------------------------------------

    @property
    def has_tdensity(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.tdensity)

    # -----------------------------------------------------------------

    @property
    def tdensity(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.tdensity)

    # -----------------------------------------------------------------

    @property
    def ncell_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.cell_properties)

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.cell_properties)

    # -----------------------------------------------------------------

    @property
    def has_single_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.cell_properties)

    # -----------------------------------------------------------------

    @property
    def cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.cell_properties)

    # -----------------------------------------------------------------

    @property
    def single_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.cell_properties)

    # -----------------------------------------------------------------

    @property
    def nstellar_density(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.stellar_density)

    # -----------------------------------------------------------------

    @property
    def has_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.stellar_density)

    # -----------------------------------------------------------------

    @property
    def has_single_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.stellar_density)

    # -----------------------------------------------------------------

    @property
    def stellar_density(self):

        """
        Thisn function ...
        :return:
        """

        return self.get_files(self._output_types.stellar_density)

    # -----------------------------------------------------------------

    @property
    def single_stellar_density(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_single_file(self._output_types.stellar_density)

    # -----------------------------------------------------------------

    @property
    def ntree(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.tree)

    # -----------------------------------------------------------------

    @property
    def has_tree(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_files(self._output_types.tree)

    # -----------------------------------------------------------------

    @property
    def tree(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.tree)

    # -----------------------------------------------------------------

    @property
    def nconvergence(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.convergence)

    # -----------------------------------------------------------------

    @property
    def has_convergence(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.convergence)

    # -----------------------------------------------------------------

    @property
    def convergence(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.convergence)

    # -----------------------------------------------------------------

    @property
    def ndust_mass(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_mass)

    # -----------------------------------------------------------------

    @property
    def has_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_mass)

    # -----------------------------------------------------------------

    @property
    def dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_mass)

    # -----------------------------------------------------------------

    @property
    def ndust_mix_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_mix_properties)

    # -----------------------------------------------------------------

    @property
    def has_dust_mix_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_mix_properties)

    # -----------------------------------------------------------------

    @property
    def dust_mix_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_mix_properties)

    # -----------------------------------------------------------------

    @property
    def ndust_optical_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_optical_properties)

    # -----------------------------------------------------------------

    @property
    def has_dust_optical_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_optical_properties)

    # -----------------------------------------------------------------

    @property
    def dust_optical_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_optical_properties)

    # -----------------------------------------------------------------

    @property
    def ndust_grain_sizes(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.dust_grain_sizes)

    # -----------------------------------------------------------------

    @property
    def has_dust_grain_sizes(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.dust_grain_sizes)

    # -----------------------------------------------------------------

    @property
    def dust_grain_sizes(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.dust_grain_sizes)

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.parameters)

    # -----------------------------------------------------------------

    @property
    def has_parameters(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.parameters)

    # -----------------------------------------------------------------

    @property
    def has_single_parameters(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.parameters)

    # -----------------------------------------------------------------

    @property
    def parameters(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.parameters)

    # -----------------------------------------------------------------

    @property
    def single_parameters(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.parameters)

# -----------------------------------------------------------------

# The various extraction output types
extraction_output_types = Map()
extraction_output_types.timeline = "timeline"
extraction_output_types.memory = "memory"

# -----------------------------------------------------------------

# Description of the extraction output types
extraction_output_type_choices = dict()
extraction_output_type_choices[extraction_output_types.timeline] = "simulation phase timeline table"
extraction_output_type_choices[extraction_output_types.memory] = "simulation phase memory usage table"

# -----------------------------------------------------------------

def get_extraction_output_type(filename):

    """
    This function ...
    :param filename:
    :return:
    """

    # Cached file
    if filename == cached_filename: return cached_type_name

    ## Timeline
    elif filename == "timeline.dat": return extraction_output_types.timeline

    ## Memory
    elif filename == "memory.dat": return extraction_output_types.memory

    # No match
    else: return None

# -----------------------------------------------------------------

class ExtractionOutput(Output):

    """
    This class ...
    """

    _output_types = extraction_output_types
    _output_type_choices = extraction_output_type_choices

    # -----------------------------------------------------------------

    def __init__(self, *args):

        """
        This function ...
        :param args:
        """

        # Call the constructor of the base class
        super(ExtractionOutput, self).__init__()

        # Load the file paths
        self.load_files(*args)

    # -----------------------------------------------------------------

    @staticmethod
    def get_type(filename, directory=None):

        """
        This function ...
        :param filename:
        :param directory:
        :return:
        """

        return get_extraction_output_type(filename)

    # -----------------------------------------------------------------

    @staticmethod
    def get_directory_type(dirname, directory=None):

        """
        This function ...
        :param dirname:
        :param directory:
        :return:
        """

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the filepaths
        filepaths = fs.files_in_path(path)

        # Load
        return cls(*filepaths)

    # -----------------------------------------------------------------

    @property
    def ntimeline(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.timeline)

    # -----------------------------------------------------------------

    @property
    def has_timeline(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.timeline)

    # -----------------------------------------------------------------

    @property
    def has_single_timeline(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.timeline)

    # -----------------------------------------------------------------

    @property
    def timeline(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.timeline)

    # -----------------------------------------------------------------

    @property
    def single_timeline(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.timeline)

    # -----------------------------------------------------------------

    @property
    def nmemory(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.memory)

    # -----------------------------------------------------------------

    @property
    def has_memory(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.memory)

    # -----------------------------------------------------------------

    @property
    def has_single_memory(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.memory)

    # -----------------------------------------------------------------

    @property
    def memory(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.memory)

    # -----------------------------------------------------------------

    @property
    def single_memory(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.memory)

# -----------------------------------------------------------------

# The various plotting output types
plotting_output_types = Map()
plotting_output_types.seds = "seds"

# Description of the plotting output types
plotting_output_type_choices = dict()
plotting_output_type_choices[plotting_output_types.seds] = "plots of the simulated SEDs"

# -----------------------------------------------------------------

def get_plotting_output_type(filename):

    """
    This function ...
    :param filename:
    :return:
    """

    # Cached file
    if filename == cached_filename: return cached_type_name

    # SEDs
    elif "sed" in filename and (filename.endswith("png") or filename.endswith("pdf")): return plotting_output_types.seds

    # No match
    else: return None

# -----------------------------------------------------------------

class PlottingOutput(Output):

    """
    This function ...
    """

    _output_types = plotting_output_types
    _output_type_choices = plotting_output_type_choices
    _with_extensions = ["pdf", "png"]

    # -----------------------------------------------------------------

    def __init__(self, *args):

        """
        This function ...
        :param args:
        """

        # Call the constructor of the base class
        super(PlottingOutput, self).__init__()

        # Load the file paths
        self.load_files(*args)

    # -----------------------------------------------------------------

    @staticmethod
    def get_type(filename, directory=None):

        """
        This function ...
        :param filename:
        :param directory:
        :return:
        """

        return get_plotting_output_type(filename)

    # -----------------------------------------------------------------

    @staticmethod
    def get_directory_type(dirname, directory=None):

        """
        This function ...
        :param dirname:
        :param directory:
        :return:
        """

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get files
        filepaths = fs.files_in_path(path, extension=cls._with_extensions)
        return cls(*filepaths)

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.seds)

    # -----------------------------------------------------------------

    @property
    def seds(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.seds)

# -----------------------------------------------------------------

# The various misc output types
misc_output_types = Map()
misc_output_types.fluxes = "fluxes"
misc_output_types.image_fluxes = "image fluxes"
misc_output_types.images = "images"
misc_output_types.images_for_fluxes = "images for fluxes"
misc_output_types.fluxes_plot = "fluxes plot"
misc_output_types.images_fluxes_plot = "image fluxes plot"
misc_output_types.images_plot = "images plot"
misc_output_types.images_for_fluxes_plot = "images for fluxes plot"

# Description of the misc output types
misc_output_type_choices = dict()
misc_output_type_choices[misc_output_types.fluxes] = "mock fluxes"
misc_output_type_choices[misc_output_types.image_fluxes] = "mock fluxes from images"
misc_output_type_choices[misc_output_types.images] = "mock images"
misc_output_type_choices[misc_output_types.images_for_fluxes] = "images created for mock fluxes"
misc_output_type_choices[misc_output_types.fluxes_plot] = "plots of the mock fluxes"
misc_output_type_choices[misc_output_types.images_fluxes_plot] = "plots of the mock fluxes from images"
misc_output_type_choices[misc_output_types.images_plot] = "plot of the mock images"
misc_output_type_choices[misc_output_types.images_for_fluxes_plot] = "plot of the images created for mock fluxes"

# -----------------------------------------------------------------

def get_misc_output_type(filename, directory):

    """
    This function ...
    :param filename:
    :param directory:
    :return:
    """

    # Cached file
    if filename == cached_filename: return cached_type_name

    # Check whether directory is defined
    if directory is None: raise ValueError("Cannot determine misc output type if directory path is not specified")

    # Get directory name
    dirname = fs.name(directory)
    dirdirname = fs.name(fs.directory_of(directory))
    dirdirdirname = fs.name(fs.directory_of(fs.directory_of(directory)))

    # SED with mock fluxes
    if "fluxes" in filename and filename.endswith(".dat"):

        # From simulated fluxes or images?
        if dirname == "fluxes": return misc_output_types.fluxes
        elif dirname == "image fluxes": return misc_output_types.image_fluxes
        else: return None

    # Images in 'images' directory
    elif filename.endswith(".fits") and dirname == "images": return misc_output_types.images

    # Images in directory in 'images' directory
    elif filename.endswith(".fits") and dirdirname == "images":

        if dirdirdirname == "image fluxes": return misc_output_types.images_for_fluxes
        else: return misc_output_types.images

    # SED plot
    elif "fluxes" in filename and (filename.endswith(".png") or filename.endswith(".pdf")):

        # From simulated fluxes or images?
        if dirname == "fluxes": return misc_output_types.fluxes_plot
        elif dirname == "image fluxes": return misc_output_types.images_fluxes_plot
        else: return None

    # No match
    else: return None

# -----------------------------------------------------------------

def get_misc_directory_output_type(dirname, directory):

    """
    This function ...
    :param dirname:
    :param directory:
    :return:
    """

    # Get directory name
    dirdirname = fs.name(directory)
    dirdirdirname = fs.name(fs.directory_of(directory))

    # Fluxes directory
    if dirname == "fluxes": return misc_output_types.fluxes

    # Fluxes from images directory
    elif dirname == "image fluxes": return misc_output_types.image_fluxes

    # Images directory
    elif dirname == "images":

        if dirdirname == "image fluxes": return misc_output_types.images_for_fluxes
        else: return misc_output_types.images

    # No match
    else: return None

# -----------------------------------------------------------------

class MiscOutput(Output):

    """
    This function ...
    """

    _output_types = misc_output_types
    _output_type_choices = misc_output_type_choices
    _with_directories = True
    _recursive = True

    # -----------------------------------------------------------------

    def __init__(self, files, directories):

        """
        This function ...
        :param args:
        :param directories:
        """

        # Call the constructor of the base class
        super(MiscOutput, self).__init__()

        # Load the file paths
        self.load_files(*files)

        # Load the directory paths
        self.load_directories(*directories)

    # -----------------------------------------------------------------

    @staticmethod
    def get_type(filename, directory=None):

        """
        This function ...
        :param filename:
        :param directory:
        :return:
        """

        return get_misc_output_type(filename, directory=directory)

    # -----------------------------------------------------------------

    @classmethod
    def get_directory_type(cls, dirname, directory=None):

        """
        This function ...
        :param dirname:
        :param directory:
        :return:
        """

        return get_misc_directory_output_type(dirname, directory=directory)

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get paths
        filepaths = fs.files_in_path(path, recursive=cls._recursive)
        dirpaths = fs.directories_in_path(path, recursive=cls._recursive)

        # Create the object
        return cls(filepaths, dirpaths)

    # -----------------------------------------------------------------

    @property
    def nfluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.fluxes)

    # -----------------------------------------------------------------

    @property
    def has_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.fluxes)

    # -----------------------------------------------------------------

    @property
    def has_single_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.fluxes)

    # -----------------------------------------------------------------

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.fluxes)

    # -----------------------------------------------------------------

    @property
    def single_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.fluxes)

    # -----------------------------------------------------------------

    @property
    def nimage_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.image_fluxes)

    # -----------------------------------------------------------------

    @property
    def has_image_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.image_fluxes)

    # -----------------------------------------------------------------

    @property
    def has_single_image_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.image_fluxes)

    # -----------------------------------------------------------------

    @property
    def image_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.image_fluxes)

    # -----------------------------------------------------------------

    @property
    def single_image_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.image_fluxes)

    # -----------------------------------------------------------------

    @property
    def nimages(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.images)

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.images)

    # -----------------------------------------------------------------

    @property
    def images(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.images)

    # -----------------------------------------------------------------

    @property
    def nimages_for_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.images_for_fluxes)

    # -----------------------------------------------------------------

    @property
    def has_images_for_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.images_for_fluxes)

    # -----------------------------------------------------------------

    @property
    def images_for_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.images_for_fluxes)

    # -----------------------------------------------------------------

    @property
    def nfluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_single_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def single_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def nimages_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.images_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_images_fluxes_plot(self):

        """
        THis function ...
        :return:
        """

        return self.has_files(self._output_types.images_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_single_images_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.images_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def images_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.images_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def single_images_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.images_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def nimages_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.images_plot)

    # -----------------------------------------------------------------

    @property
    def has_images_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.images_plot)

    # -----------------------------------------------------------------

    @property
    def has_single_images_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.images_plot)

    # -----------------------------------------------------------------

    @property
    def images_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.images_plot)

    # -----------------------------------------------------------------

    @property
    def single_images_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.images_plot)

    # -----------------------------------------------------------------

    @property
    def nimages_for_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_nfiles(self._output_types.images_for_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_images_for_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_files(self._output_types.images_for_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def has_single_images_for_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.has_single_file(self._output_types.images_for_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def images_for_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_files(self._output_types.images_for_fluxes_plot)

    # -----------------------------------------------------------------

    @property
    def single_images_for_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return self.get_single_file(self._output_types.images_for_fluxes_plot)

# -----------------------------------------------------------------
