#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.data Contains the SimulationData class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from .output import SimulationOutput
from .table import SkirtTable, is_valid
from ..tools.utils import lazyproperty, LazyDictionary
from ..tools import filesystem as fs
from .skifile import SkiFile

# -----------------------------------------------------------------

class SimulationData(object):

    """
    This class ...
    """

    def __init__(self, output, coordinate_systems=None, distances=None):

        """
        The constructor ...
        :param output:
        :param coordinate_systems:
        :param distances:
        """

        # The simulation output
        self.output = output

        # Set the coordinate systems for the instruments
        self.coordinate_systems = coordinate_systems

        # Set the distances for the instruments
        self.distances = distances

    # -----------------------------------------------------------------

    @classmethod
    def from_output(cls, output, coordinate_systems=None, distances=None):

        """
        Thisnfunction ...
        :param output:
        :param coordinate_systems:
        :param distances:
        :return:
        """

        return cls(output, coordinate_systems=coordinate_systems, distances=distances)

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, prefix=None, coordinate_systems=None, distances=None):

        """
        Thinsfunction ...
        :param path:
        :param prefix:
        :param coordinate_systems:
        :param distances:
        :return:
        """

        # Get the output
        output = SimulationOutput.from_directory(path, prefix)

        # Create
        return cls(output, coordinate_systems=coordinate_systems, distances=distances)

    # -----------------------------------------------------------------

    @classmethod
    def from_cwd(cls, prefix=None, coordinate_systems=None, distances=None):

        """
        This fucntion ...
        :param prefix:
        :param coordinate_systems:
        :param distances:
        :return:
        """

        # Get the output
        output = SimulationOutput.from_cwd(prefix=prefix)

        # Create
        return cls(output, coordinate_systems=coordinate_systems, distances=distances)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, paths, coordinate_systems=None, distances=None):

        """
        This function ...
        :param paths:
        :param coordinate_systems:
        :param distances:
        :return:
        """

        # Get the output
        output = SimulationOutput(*paths)

        # Create
        return cls(output, coordinate_systems=coordinate_systems, distances=distances)

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):

        """
        This function ...
        :return:
        """

        return self.output.prefix

    # -----------------------------------------------------------------

    @property
    def ncoordinate_systems(self):

        """
        This function ...
        :return:
        """

        return len(self.coordinate_systems) if self.coordinate_systems is not None else 0

    # -----------------------------------------------------------------

    @property
    def has_coordinate_systems(self):

        """
        Thisfunction ...
        :return:
        """

        return self.ncoordinate_systems > 0

    # -----------------------------------------------------------------

    def has_coordinate_system_for_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        return self.has_coordinate_systems and instrument_name in self.coordinate_systems

    # -----------------------------------------------------------------

    def get_coordinate_system_for_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        if not self.has_coordinate_system_for_instrument(instrument_name): return None
        return self.coordinate_systems[instrument_name]

    # -----------------------------------------------------------------

    @property
    def has_skifile(self):

        """
        Thisf unction ...
        :return:
        """

        return self.output.has_single_parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def skifile(self):

        """
        This function ...
        :return:
        """

        if self.output.has_single_parameters: return SkiFile(self.output.single_parameters)
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_distances(self):

        """
        This function ...
        :return:
        """

        return self.distances is not None

    # -----------------------------------------------------------------

    def has_distance_for_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        return self.has_distances and instrument_name in self.distances

    # -----------------------------------------------------------------

    def get_distance_for_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        if self.has_distance_for_instrument(instrument_name): return self.distances[instrument_name]
        elif self.has_skifile: return self.skifile.get_instrument_distance(instrument_name)
        else: return None

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.to_string()

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.output.has_single_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_properties_path(self):

        """
        This function ...
        :return:
        """

        return self.output.single_cell_properties

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_properties(self):

        """
        This function ...
        :return:
        """

        return is_valid(self.cell_properties_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_properties(self):

        """
        This function ...
        :return:
        """

        return SkirtTable.from_file(self.cell_properties_path)

    # -----------------------------------------------------------------

    @property
    def has_isrf(self):

        """
        This function ...
        :return:
        """

        return self.output.has_isrf

    # -----------------------------------------------------------------

    @property
    def isrf_path(self):

        """
        This function ...
        :return:
        """

        return self.output.single_isrf

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_isrf(self):

        """
        This function ...
        :return:
        """

        return is_valid(self.isrf_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def isrf(self):

        """
        Thins function ...
        :return:
        """

        return SkirtTable.from_file(self.isrf_path)

    # -----------------------------------------------------------------

    @property
    def has_absorption(self):

        """
        Thisj function ...
        :return:
        """

        return self.output.has_absorption

    # -----------------------------------------------------------------

    @property
    def absorption_path(self):

        """
        This function ...
        :return:
        """

        return self.output.single_absorption

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_absorption(self):

        """
        This function ...
        :return:
        """

        return is_valid(self.absorption_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption(self):

        """
        This function ..
        :return:
        """

        return SkirtTable.from_file(self.absorption_path)

    # -----------------------------------------------------------------

    @property
    def has_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.output.has_wavelengths

    # -----------------------------------------------------------------

    @property
    def wavelengths_path(self):

        """
        This function ...
        :return:
        """

        return self.output.single_wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_wavelengths(self):

        """
        This function ...
        :return:
        """

        return is_valid(self.wavelengths_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        from .wavelengthgrid import WavelengthGrid
        return WavelengthGrid.from_skirt_output(self.wavelengths[0])

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        This function ...
        :return:
        """

        return self.output.has_seds

    # -----------------------------------------------------------------

    @property
    def sed_paths(self):

        """
        This function ...
        :return:
        """

        return self.output.seds

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_seds(self):

        """
        This function ...
        :return:
        """

        for path in self.sed_paths:
            if not is_valid(path): return False
        return True

    # -----------------------------------------------------------------

    def is_valid_sed(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        # Get the SED file path
        path = self.seds[instrument_name].get_raw("total")

        # Check validity
        return is_valid(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_paths_instruments(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        paths_instruments = dict()

        # Loop over the paths
        for path in self.sed_paths:

            # Get the instrument name
            instrument_name = get_sed_instrument_name(path, self.simulation_prefix)

            # Set the SED path for this instrument
            paths_instruments[instrument_name] = path

        # Return the dictionary
        return paths_instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def seds(self):

        """
        This function ...
        :return:
        """

        from ..data.sed import SED
        from ..plot.simulationseds import number_of_columns, contributions

        # Initialize dictionary
        seds_instruments = dict()

        # Loop over the paths
        for path in self.sed_paths:

            # Get the instrument name
            instrument_name = get_sed_instrument_name(path, self.simulation_prefix)

            # Check how many columns the SED file contains
            ncols = number_of_columns(path)

            # Lazy dictionary of SEDs for this instrument
            seds = LazyDictionary(SED.from_skirt)

            # Check the type of the Instrument / SED
            if ncols == 2: seds["total"] = path

            # More columns
            else:

                # Loop over the different contributions
                for contribution in contributions:

                    # Set
                    seds.set(contribution, path, contribution=contribution)

            # Set the SEDs for this instrument
            seds_instruments[instrument_name] = seds

        # Return the SEDs
        return seds_instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def has_instruments(self):

        """
        This function ...
        :return:
        """

        return self.has_seds or self.has_images

    # -----------------------------------------------------------------

    @lazyproperty
    def instrument_names(self):

        """
        This function ...
        :return:
        """

        # Have SED files
        if self.has_seds: return self.seds.keys()

        # Have datacubes
        elif self.has_images: return self.images.keys()

        # Neither SED files or datacubes
        else: raise ValueError("Cannot determine instrument names")

    # -----------------------------------------------------------------

    @lazyproperty
    def has_wavelength_grid(self):

        """
        This function ..
        :return:
        """

        return self.has_wavelengths or self.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        Thisf unction ...
        :return:
        """

        # Wavelength grid file
        if self.has_wavelengths: return self.wavelength_grid

        # There are SED files
        elif self.has_seds:

            from .wavelengthgrid import WavelengthGrid

            # Get one sed
            sed = self.seds[self.instrument_names[0]]["total"]

            # Create wavelength grid
            wavelength_grid = WavelengthGrid.from_sed(sed)

            # Return the wavelength grid
            return wavelength_grid

        # No wavelength grid file, no SEDs
        else: raise ValueError("Cannot get wavelength grid")

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.output.has_images

    # -----------------------------------------------------------------

    @property
    def image_paths(self):

        """
        This function ...
        :return:
        """

        return self.output.images

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_images(self):

        """
        This function ...
        :return:
        """

        from ...magic.core.fits import is_valid as is_valid_fits

        for path in self.image_paths:
            if not is_valid_fits(path): return False

        return True

    # -----------------------------------------------------------------

    def is_valid_image(self, instrument_name, contribution):

        """
        This function ...
        :param instrument_name:
        :param contribution:
        :return:
        """

        # Get the file path
        path = self.images[instrument_name].get_raw(contribution)

        from ...magic.core.fits import is_valid as is_valid_fits

        # Check validity
        return is_valid_fits(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def image_paths_instruments(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        paths_instruments = defaultdict(list)

        # Loop over the paths
        for path in self.image_paths:

            # Get the instrument name
            instrument_name = get_datacube_instrument_name(path, self.simulation_prefix)

            # Set the datacube path for this instrument
            paths_instruments[instrument_name].append(path)

        # Return the dictionary
        return paths_instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def images(self):

        """
        This function ...
        :return:
        """

        from ..misc.datacubes import get_datacube_instrument_name
        from ...magic.core.datacube import DataCube

        # Initialize dictionary
        images_instruments = dict()

        # Loop over the paths
        for path in self.image_paths:

            # Get the instrument name
            instrument_name = get_datacube_instrument_name(path, self.simulation_prefix)

            # Get the contribution
            filename = fs.strip_extension(fs.name(path))
            contribution = filename.split("_")[-1]

            # Initialize lazy dictionary for the datacubes of this instrument
            if instrument_name not in images_instruments:
                distance = self.get_distance_for_instrument(instrument_name)
                wcs = self.get_coordinate_system_for_instrument(instrument_name)
                images = LazyDictionary(DataCube.from_file, wavelength_grid=self.wavelength_grid, distance=distance, wcs=wcs)
                images_instruments[instrument_name] = images

            # Add to dictionary
            images_instruments[instrument_name][contribution] = path

        # Return the images
        return images_instruments

    # -----------------------------------------------------------------

    @property
    def has_stellar_density(self):

        """
        This funtion ...
        :return:
        """

        return self.output.has_stellar_density

    # -----------------------------------------------------------------

    @property
    def stellar_density_path(self):

        """
        This function ...
        :return:
        """

        return self.output.single_stellar_density

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_stellar_density(self):

        """
        This function ...
        :return:
        """

        return is_valid(self.stellar_density_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_density(self):

        """
        This function ...
        :return:
        """

        return SkirtTable.from_file(self.stellar_density_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_any(self):

        """
        This function ...
        :return:
        """

        # Check
        if self.has_isrf: return True
        elif self.has_absorption: return True
        elif self.has_stellar_density: return True
        elif self.has_seds: return True
        elif self.has_images: return True

        # Nothing
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def valid(self):

        """
        This function ...
        :return:
        """

        # ISRF
        if self.has_isrf and not self.valid_isrf: return False

        # Absorption
        if self.has_absorption and not self.valid_absorption: return False

        # Stellar density
        if self.has_stellar_density and not self.valid_stellar_density: return False

        # SEDs
        if self.has_seds and not self.valid_seds: return False

        # Images
        if self.has_images and not self.valid_images: return False

        # ALl checks passed
        return True

    # -----------------------------------------------------------------

    def to_string(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        from ..tools import formatting as fmt

        lines = []

        # ISRF
        if self.has_isrf:

            lines.append(line_prefix)

            # Make line
            if self.valid_isrf: line = fmt.green + fmt.underlined + "ISRF" + fmt.reset
            else: line = fmt.red + fmt.underlined + "ISRF: invalid" + fmt.reset

            lines.append(line_prefix + line)
            lines.append(line_prefix)

        # Absorption
        if self.has_absorption:

            lines.append(line_prefix)

            # Make line
            if self.valid_absorption: line = fmt.green + fmt.underlined + "Absorption" + fmt.reset
            else: line = fmt.red + fmt.underlined + "Absorption: invalid" + fmt.reset

            lines.append(line_prefix + line)
            lines.append(line_prefix)

        # Stellar density
        if self.has_stellar_density:

            lines.append(line_prefix)

            # Make line
            if self.valid_stellar_density: line = fmt.green + fmt.underlined + "Stellar density" + fmt.reset
            else: line = fmt.red + fmt.underlined + "Stellar density: invalid" + fmt.reset

            lines.append(line_prefix + line)
            lines.append(line_prefix)

        # SEDs
        if self.has_seds:

            lines.append(line_prefix)

            # Make line
            if self.valid_seds: line = fmt.green + fmt.underlined + "SEDs:" + fmt.reset
            else: line = fmt.red + fmt.underlined + "SEDs:" + fmt.reset

            lines.append(line_prefix + line)
            lines.append(line_prefix)

            # Loop over the instruments
            for instrument_name in self.seds:

                # Check whether file is valid
                is_valid = self.is_valid_sed(instrument_name)

                #
                if is_valid: line = " - " + fmt.green + instrument_name + fmt.reset
                else: line = " - " + fmt.red + instrument_name + ": invalid" + fmt.reset

                # Add the line
                lines.append(line_prefix + line)

                # Valid? -> show contributions
                if is_valid:

                    lines.append(line_prefix)

                    # Loop over the contributions
                    for contribution in self.seds[instrument_name]:

                        line = "    * " + contribution

                        # Add the line
                        lines.append(line_prefix + line)

                    lines.append(line_prefix)

        # Images
        if self.has_images:

            lines.append(line_prefix)

            # Make line
            if self.valid_images: line = fmt.green + fmt.underlined + "Images:" + fmt.reset
            else: line = fmt.red + fmt.underlined + "Images:" + fmt.reset

            lines.append(line_prefix + line)
            lines.append(line_prefix)

            # Loop over the instruments
            for instrument_name in self.images:

                lines.append(line_prefix + " - " + instrument_name + ":")
                lines.append(line_prefix)

                # Loop over the contributions
                for contribution in self.images[instrument_name]:

                    if self.is_valid_image(instrument_name, contribution): line = "    * " + fmt.green + contribution + fmt.reset
                    else: line = "    * " + fmt.red + contribution + ": invalid" + fmt.reset

                    # Add the line
                    lines.append(line_prefix + line)

                lines.append(line_prefix)

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

def get_sed_instrument_name(sed_path, prefix):

    """
    This function ...
    :param sed_path:
    :param prefix:
    :return:
    """

    return fs.name(sed_path).split("_sed.dat")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------

def get_datacube_instrument_name(datacube_path, prefix):

    """
    This function ...
    :param datacube_path:
    :param prefix:
    :return:
    """

    # For all
    return fs.name(datacube_path).split(prefix + "_")[1].rsplit("_", 1)[0]

# -----------------------------------------------------------------
