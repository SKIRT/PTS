#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.output Contains the SimulationOutput class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..basics.map import Map

# -----------------------------------------------------------------

# The various output types
output_types = Map()
output_types.isrf = "isrf"
output_types.absorption = "abs"
output_types.temperature = "temp"
output_types.seds = "sed"
output_types.images = "image"
output_types.total_images = "image-total"
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
output_types.tree = "tree"
output_types.convergence = "convergence"

# -----------------------------------------------------------------

# Description of the output types
output_type_choices = dict()
output_type_choices[output_types.isrf] = "interstellar radiation field strength"
output_type_choices[output_types.absorption] = "absorption luminosities"
output_type_choices[output_types.temperature] = "temperature"
output_type_choices[output_types.seds] = "all SEDs"
output_type_choices[output_types.images] = "all datacubes"
output_type_choices[output_types.total_images] = "datacubes of total emission"
output_type_choices[output_types.direct_images] = "datacubes of direct emission"
output_type_choices[output_types.transparent_images] = "datacubes of transparent emission"
output_type_choices[output_types.scattered_images] = "datacubes of scattered emission"
output_type_choices[output_types.dust_images] = "datacubes of dust emission"
output_type_choices[output_types.dust_scattered_images] = "datacubes of scattered dust emission"
output_type_choices[output_types.temperature] = "temperature per dust cell"
output_type_choices[output_types.logfiles] = "log files"
output_type_choices[output_types.wavelengths] = "wavelength files"
output_type_choices[output_types.grid] = "grid files"
output_type_choices[output_types.gdensity] = "grid dust density"
output_type_choices[output_types.tdensity] = "theoretical dust density"
output_type_choices[output_types.tree] = "dust grid tree data file"
output_type_choices[output_types.convergence] = "convergence file"

# -----------------------------------------------------------------

def get_output_type(filename):

    """
    This function ...
    :param filename:
    :return:
    """

    ## ISRF
    if filename.endswith("_ds_isrf.dat"): return output_types.isrf

    ## Absorption
    elif filename.endswith("_ds_abs.dat"): return output_types.absorption

    ## Temperature
    elif "_ds_temp" in filename and filename.endswith(".fits"): return output_types.temperature

    ## SED
    elif filename.endswith("_sed.dat"): return output_types.seds

    ## Total datacubes
    elif filename.endswith("_total.fits"): return output_types.total_images

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

    ## Dust grid onvergence
    elif filename.endswith("_ds_convergence.dat"): return output_types.convergence

    ## Dust grid tree data
    elif "_ds_tree" in filename and filename.endswith(".dat"): return output_types.tree

    # No match
    return None

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

class SimulationOutput(object):

    """
    This class ...
    """

    def __init__(self, *args):

        """
        The constructor ...
        :param args:
        """

        # The simulation prefix
        self.prefix = None

        # Dictionary of filepaths
        self.paths = defaultdict(list)

        # Loop over the filepaths, categorize
        for filepath in args:
            filename = fs.name(filepath)
            prefix = filename.split("_")[0]
            if self.prefix is None: self.prefix = prefix
            elif self.prefix != prefix: raise ValueError("Cannot add files with different simulation prefixes")
            output_type = get_output_type(filename)
            if output_type is None: self.paths[other_name].append(filepath)
            else: self.paths[output_type].append(filepath)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, paths):

        """
        This function ...
        :param paths:
        :return:
        """

        return cls(*paths)

    # -----------------------------------------------------------------

    @classmethod
    def from_cwd(cls, prefix=None):

        """
        This function ...
        :param prefix:
        :return:
        """

        return cls.from_directory(fs.cwd(), prefix=prefix)

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, prefix=None):

        """
        This function ...
        :param path:
        :param prefix:
        :return:
        """

        return cls.from_paths(fs.files_in_path(path, startswith=prefix))

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

        return cls.from_paths(remote.files_in_path(path, startswith=prefix))

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for path in self.isrf: yield path
        for path in self.absorption: yield path
        for path in self.temperature: yield path
        for path in self.seds: yield path
        for path in self.total_images: yield path
        for path in self.direct_images: yield path
        for path in self.transparent_images: yield path
        for path in self.scattered_images: yield path
        for path in self.dust_images: yield path
        for path in self.dust_scattered_images: yield path
        for path in self.cell_temperature: yield path
        for path in self.logfiles: yield path
        for path in self.wavelengths: yield path
        for path in self.grid: yield path
        for path in self.gdensity: yield path
        for path in self.tdensity: yield path
        for path in self.tree: yield path
        for path in self.convergence: yield path
        for path in self.other: yield path

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        total = 0
        if self.has_isrf: total += self.nisrf
        if self.has_absorption: total += self.nabsorption
        if self.has_temperature: total += self.ntemperature
        if self.has_seds: total += self.nseds
        if self.has_total_images: total += self.ntotal_images
        if self.has_direct_images: total += self.ndirect_images
        if self.has_transparent_images: total += self.ntransparent_images
        if self.has_scattered_images: total += self.nscattered_images
        if self.has_dust_images: total += self.ndust_images
        if self.has_dust_scattered_images: total += self.ndust_scattered_images
        if self.has_cell_temperature: total += self.ncell_temperature
        if self.has_logfiles: total += self.nlogfiles
        if self.has_wavelengths: total += self.nwavelengths
        if self.has_grid: total += self.ngrid
        if self.has_gdensity: total += self.ngdensity
        if self.has_tdensity: total += self.ntdensity
        if self.has_tree: total += self.ntree
        if self.has_convergence: total += self.nconvergence
        if self.has_other: total += self.nother
        return total

    # -----------------------------------------------------------------

    @property
    def nisrf(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.isrf])

    # -----------------------------------------------------------------

    @property
    def has_isrf(self):

        """
        This function ...
        :return:
        """

        return output_types.isrf in self.paths and self.nisrf > 0

    # -----------------------------------------------------------------

    @property
    def isrf(self):

        """
        This function ...
        :return:
        """

        if not self.has_isrf: return []
        return self.paths[output_types.isrf]

    # -----------------------------------------------------------------

    @property
    def nabsorption(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.absorption])

    # -----------------------------------------------------------------

    @property
    def has_absorption(self):

        """
        This function ...
        :return:
        """

        return output_types.absorption in self.paths and self.nabsorption > 0

    # -----------------------------------------------------------------

    @property
    def absorption(self):

        """
        This function ...
        :return:
        """

        if not self.has_absorption: return []
        return self.paths[output_types.absorption]

    # -----------------------------------------------------------------

    @property
    def ntemperature(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.temperature])

    # -----------------------------------------------------------------

    @property
    def has_temperature(self):

        """
        This function ...
        :return:
        """

        return output_types.temperature in self.paths and self.ntemperature > 0

    # -----------------------------------------------------------------

    @property
    def temperature(self):

        """
        This function ...
        :return:
        """

        if not self.has_temperature: return []
        return self.paths[output_types.temperature]

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.seds])

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        This function ...
        :return:
        """

        return output_types.seds in self.paths and self.nseds > 0

    # -----------------------------------------------------------------

    @property
    def seds(self):

        """
        This function ...
        :return:
        """

        if not self.has_seds: return []
        return self.paths[output_types.seds]

    # -----------------------------------------------------------------

    @property
    def nimages(self):

        """
        This function ...
        :return:
        """

        total = 0
        if self.has_total_images: total += self.ntotal_images
        if self.has_direct_images: total += self.ndirect_images
        if self.has_transparent_images: total += self.ntransparent_images
        if self.has_scattered_images: total += self.nscattered_images
        if self.has_dust_images: total += self.ndust_images
        if self.has_dust_scattered_images: total += self.ndust_scattered_images
        return total

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.nimages > 0

    # -----------------------------------------------------------------

    @property
    def images(self):

        """
        This function ...
        :return:
        """

        paths = []
        if self.has_total_images: paths.extend(self.total_images)
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

        return len(self.paths[output_types.total_images])

    # -----------------------------------------------------------------

    @property
    def has_total_images(self):

        """
        This function ...
        :return:
        """

        return output_types.total_images in self.paths and self.ntotal_images > 0

    # -----------------------------------------------------------------

    @property
    def total_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_total_images: return []
        else: return self.paths[output_types.total_images]

    # -----------------------------------------------------------------

    @property
    def ndirect_images(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.direct_images])

    # -----------------------------------------------------------------

    @property
    def has_direct_images(self):

        """
        This function ...
        :return:
        """

        return output_types.direct_images in self.paths and self.ndirect_images > 0

    # -----------------------------------------------------------------

    @property
    def direct_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_direct_images: return []
        else: return self.paths[output_types.direct_images]

    # -----------------------------------------------------------------

    @property
    def ntransparent_images(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.transparent_images])

    # -----------------------------------------------------------------

    @property
    def has_transparent_images(self):

        """
        This function ...
        :return:
        """

        return output_types.transparent_images in self.paths and self.ntransparent_images > 0

    # -----------------------------------------------------------------

    @property
    def transparent_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_transparent_images: return []
        else: return self.paths[output_types.transparent_images]

    # -----------------------------------------------------------------

    @property
    def nscattered_images(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.scattered_images])

    # -----------------------------------------------------------------

    @property
    def has_scattered_images(self):

        """
        This function ...
        :return:
        """

        return output_types.scattered_images in self.paths and self.nscattered_images > 0

    # -----------------------------------------------------------------

    @property
    def scattered_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_scattered_images: return []
        else: return self.paths[output_types.scattered_images]

    # -----------------------------------------------------------------

    @property
    def ndust_images(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.dust_images])

    # -----------------------------------------------------------------

    @property
    def has_dust_images(self):

        """
        This function ...
        :return:
        """

        return output_types.dust_images in self.paths and self.ndust_images > 0

    # -----------------------------------------------------------------

    @property
    def dust_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_dust_images: return []
        else: return self.paths[output_types.dust_images]

    # -----------------------------------------------------------------

    @property
    def ndust_scattered_images(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.dust_scattered_images])

    # -----------------------------------------------------------------

    @property
    def has_dust_scattered_images(self):

        """
        This function ...
        :return:
        """

        return output_types.dust_scattered_images in self.paths and self.ndust_scattered_images > 0

    # -----------------------------------------------------------------

    @property
    def dust_scattered_images(self):

        """
        This function ...
        :return:
        """

        if not self.has_dust_scattered_images: return []
        else: return self.paths[output_types.dust_scattered_images]

    # -----------------------------------------------------------------

    @property
    def ncell_temperature(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.cell_temperature])

    # -----------------------------------------------------------------

    @property
    def has_cell_temperature(self):

        """
        This function ...
        :return:
        """

        return output_types.cell_temperature in self.paths and self.ncell_temperature > 0

    # -----------------------------------------------------------------

    @property
    def cell_temperature(self):

        """
        This function ...
        :return:
        """

        if not self.has_cell_temperature: return []
        else: return self.paths[output_types.cell_temperature]

    # -----------------------------------------------------------------

    @property
    def nlogfiles(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.logfiles])

    # -----------------------------------------------------------------

    @property
    def has_logfiles(self):

        """
        This funtion ...
        :return:
        """

        return output_types.logfiles in self.paths and self.nlogfiles > 0

    #

    @property
    def logfiles(self):

        """
        This function ...
        :return:
        """

        if not self.has_logfiles: return []
        else: return self.paths[output_types.logfiles]

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.wavelengths])

    # -----------------------------------------------------------------

    @property
    def has_wavelengths(self):

        """
        This function ...
        :return:
        """

        return output_types.wavelengths in self.paths and self.nwavelengths > 0

    # -----------------------------------------------------------------

    @property
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        if not self.has_wavelengths: return []
        else: return self.paths[output_types.wavelengths]

    # -----------------------------------------------------------------

    @property
    def ngrid(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.grid])

    # -----------------------------------------------------------------

    @property
    def has_grid(self):

        """
        Thisf unction ...
        :return:
        """

        return output_types.grid in self.paths and self.ngrid > 0

    # -----------------------------------------------------------------

    @property
    def grid(self):

        """
        This function ...
        :return:
        """

        if not self.has_grid: return []
        else: return self.paths[output_types.grid]

    # -----------------------------------------------------------------

    @property
    def ngdensity(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.gdensity])

    # -----------------------------------------------------------------

    @property
    def has_gdensity(self):

        """
        This function ...
        :return:
        """

        return output_types.gdensity in self.paths and self.ngdensity > 0

    # -----------------------------------------------------------------

    @property
    def gdensity(self):

        """
        This fucntion ...
        :return:
        """

        if not self.has_gdensity: return []
        else: return self.paths[output_types.gdensity]

    # -----------------------------------------------------------------

    @property
    def ntdensity(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.tdensity])

    # -----------------------------------------------------------------

    @property
    def has_tdensity(self):

        """
        This function ...
        :return:
        """

        return output_types.tdensity in self.paths and self.ntdensity > 0

    # -----------------------------------------------------------------

    @property
    def tdensity(self):

        """
        This function ...
        :return:
        """

        if not self.has_tdensity: return []
        else: return self.paths[output_types.tdensity]

    # -----------------------------------------------------------------

    @property
    def ntree(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.tree])

    @property
    def has_tree(self):

        """
        Thisf unction ...
        :return:
        """

        return output_types.tree in self.paths and self.ntree > 0

    @property
    def tree(self):

        """
        This function ...
        :return:
        """

        if not self.has_tree: return []
        else: return self.paths[output_types.tree]

    # -----------------------------------------------------------------

    @property
    def nconvergence(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[output_types.convergence])

    @property
    def has_convergence(self):

        """
        This function ...
        :return:
        """

        return output_types.convergence in self.paths and self.nconvergence > 0

    @property
    def convergence(self):

        """
        This function ...
        :return:
        """

        if not self.has_convergence: return []
        else: return self.paths[output_types.convergence]

    # -----------------------------------------------------------------

    @property
    def nother(self):

        """
        This function ...
        :return:
        """

        return len(self.paths[other_name])

    # -----------------------------------------------------------------

    @property
    def has_other(self):

        """
        This function ...
        :return:
        """

        return other_name in self.paths and self.nother > 0

    # -----------------------------------------------------------------

    @property
    def other(self):

        """
        This function ...
        :return:
        """

        if not self.has_other: return []
        else: return self.paths[other_name]

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        from ..tools import formatting as fmt

        lines = []

        # ISRF
        if self.has_isrf:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.isrf].title() + fmt.reset + " (" + str(self.nisrf) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.isrf: lines.append(" - " + path)

        # Absorption
        if self.has_absorption:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.absorption].title() + fmt.reset + " (" + str(self.nabsorption) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.absorption: lines.append(" - " + path)

        # Temperature
        if self.has_temperature:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.temperature].title() + fmt.reset + " (" + str(self.ntemperature) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.temperature: lines.append(" - " + path)

        # SEDs
        if self.has_seds:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.seds].title() + fmt.reset + " (" + str(self.nseds) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.seds: lines.append(" - " + path)

        # Total images
        if self.has_total_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.total_images].title() + fmt.reset + " (" + str(self.ntotal_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.total_images: lines.append(" - " + path)

        # Direct images
        if self.has_direct_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.direct_images].title() + fmt.reset + " (" + str(self.ndirect_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.direct_images: lines.append(" - " + path)

        # Transparent images
        if self.has_transparent_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.transparent_images].title() + fmt.reset + " (" + str(self.ntransparent_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.transparent_images: lines.append(" - " + path)

        # Scattered images
        if self.has_scattered_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.scattered_images].title() + fmt.reset + " (" + str(self.nscattered_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.scattered_images: lines.append(" - " + path)

        # Dust images
        if self.has_dust_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.dust_images].title() + fmt.reset + " (" + str(self.ndust_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.dust_images: lines.append(" - " + path)

        # Dust scattered images
        if self.has_dust_scattered_images:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.dust_scattered_images].title() + fmt.reset + " (" + str(self.ndust_scattered_images) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.dust_scattered_images: lines.append(" - " + path)

        # Cell temperature
        if self.has_cell_temperature:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.cell_temperature].title() + fmt.reset + " (" + str(self.ncell_temperature) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.cell_temperature: lines.append(" - " + path)

        # Logfiles
        if self.has_logfiles:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.logfiles].title() + fmt.reset + " (" + str(self.nlogfiles) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.logfiles: lines.append(" - " + path)

        # Wavelengths
        if self.has_wavelengths:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.wavelengths].title() + fmt.reset + " (" + str(self.nwavelengths) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.wavelengths: lines.append(" - " + path)

        # Grid
        if self.has_grid:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.grid].title() + fmt.reset + " (" + str(self.ngrid) + "):"
            lines.append(title)
            lines.append("")

            # Add path
            for path in self.grid: lines.append(" - " + path)

        # Density
        if self.has_gdensity:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.gdensity].title() + fmt.reset + " (" + str(self.ngdensity ) + "):"
            lines.append(title)
            lines.append("")

            # Add path
            for path in self.gdensity: lines.append(" - " + path)

        # Density
        if self.has_tdensity:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.tdensity].title() + fmt.reset + " (" + str(self.ntdensity) + "):"
            lines.append(title)
            lines.append("")

            # Add path
            for path in self.tdensity: lines.append(" - " + path)

        # Tree
        if self.has_tree:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.tree].title() + fmt.reset + " (" + str(self.ntree) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.tree: lines.append(" - " + path)

        # Convergence
        if self.has_convergence:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + output_type_choices[output_types.convergence].title() + fmt.reset + " (" + str(self.nconvergence) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.convergence: lines.append(" - " + path)

        # Other
        if self.has_other:
            lines.append("")

            # Add title
            title = fmt.green + fmt.underlined + "Other output" + fmt.reset + " (" + str(self.nother) + "):"
            lines.append(title)
            lines.append("")

            # Add paths
            for path in self.other: lines.append(" - " + path)

        # Add new line
        lines.append("")

        # Return
        return "\n".join(lines)

# -----------------------------------------------------------------
