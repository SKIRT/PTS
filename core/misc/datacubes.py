#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.datacubes Contains the abstract DatacubesMiscMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..basics.log import log
from ..simulation.simulation import createsimulations
from ..tools import filesystem as fs
from ..tools.utils import lazyproperty
from ..simulation.wavelengthgrid import WavelengthGrid
from ..simulation.output import SimulationOutput
from ..simulation.skifile import SkiFile
from ..basics.containers import DefaultOrderedDict

# -----------------------------------------------------------------

all_contributions = ["total", "direct", "scattered", "dust", "dustscattered", "transparent", "counts"]

# -----------------------------------------------------------------

def find_datacubes_in_cwd(return_prefix=False, contributions=("total",)):

    """
    This function ...
    :param return_prefix:
    :param contributions:
    :return:
    """

    return find_datacubes(fs.cwd(), return_prefix=return_prefix, contributions=contributions)

# -----------------------------------------------------------------

def find_datacubes(output_path, return_prefix=False, contributions=("total",)):

    """
    This function ...
    :param output_path:
    :param return_prefix:
    :param contributions:
    :return:
    """

    # Create simulation output
    output = SimulationOutput.from_directory(output_path)

    # Return the paths
    paths = []
    if "total" in contributions: paths.extend(output.total_images)
    if "direct" in contributions: paths.extend(output.direct_images)
    if "scattered" in contributions: paths.extend(output.scattered_images)
    if "dust" in contributions: paths.extend(output.dust_images)
    if "dustscattered" in contributions: paths.extend(output.dust_scattered_images)
    if "transparent" in contributions: paths.extend(output.transparent_images)
    if "counts" in contributions: paths.extend(output.count_images)

    # Show the datacubes
    log.debug("Found datacubes:")
    for path in paths:
        filename = fs.name(path)
        log.debug(" - " + filename)

    # Return the paths
    if return_prefix: return paths, output.prefix
    else: return paths

# -----------------------------------------------------------------

def find_seds_in_cwd(return_prefix=False):

    """
    This function ...
    :param return_prefix:
    :return:
    """

    return find_seds(fs.cwd(), return_prefix=return_prefix)

# -----------------------------------------------------------------

def find_seds(output_path, return_prefix=False):

    """
    This function ...
    :param output_path:
    :param return_prefix:
    :return:
    """

    # Create simulation output
    output = SimulationOutput.from_directory(output_path)
    paths = output.seds

    # Show the SEDs
    log.debug("Found SEDS:")
    for path in paths:
        filename = fs.name(path)
        log.debug(" - " + filename)

    # Return the paths
    if return_prefix: return paths, output.prefix
    else: return paths

# -----------------------------------------------------------------

def load_parameters_in_cwd(prefix=None):

    """
    This function ...
    :param prefix:
    :return:
    """

    return load_parameters(fs.cwd(), prefix=prefix)

# -----------------------------------------------------------------

def load_parameters(output_path, prefix=None):

    """
    This function ...
    :param output_path:
    :param prefix:
    :return:
    """

    parameters_path = find_parameters_path(output_path, prefix=prefix)
    return SkiFile(parameters_path)

# -----------------------------------------------------------------

def load_instruments(output_path, prefix=None):

    """
    This function ...
    :param output_path:
    :param prefix:
    :return:
    """

    # Load the parameters
    ski = load_parameters(output_path, prefix=prefix)

    # Load the instruments
    instruments = OrderedDict()
    for name in ski.instrumentnames(): instruments[name] = ski.get_instrument_object(name)

    # Return
    return instruments

# -----------------------------------------------------------------

def load_instrument_distances(output_path, prefix=None):

    """
    This function ...
    :param output_path:
    :param prefix:
    :return:
    """

    # Load the parameters
    ski = load_parameters(output_path, prefix=prefix)

    # Load the distances
    distances = OrderedDict()
    for name in ski.instrumentnames(): distances[name] = ski.get_instrument_object(name).distance

    # Return
    return distances

# -----------------------------------------------------------------

def find_parameters_path(output_path, prefix=None):

    """
    This function ...
    :param output_path:
    :param prefix:
    :return:
    """

    # Look for single ski file
    filepaths = fs.files_in_path(output_path, extension="ski")
    nfiles = len(filepaths)
    if nfiles == 1: return filepaths[0]
    elif prefix is not None:
        # Find ski file path with prefix
        for filepath in filepaths:
            name = fs.strip_extension(fs.name(filepath))
            if name == prefix: return filepath

    # Create simulation output
    output = SimulationOutput.from_directory(output_path, prefix=prefix)

    # Look for single ski file that corresponds to the output prefix
    if prefix is None: prefix = output.prefix
    for filepath in filepaths:
        name = fs.strip_extension(fs.name(filepath))
        if name == prefix: return filepath

    # Has parameters file
    if output.has_single_parameters: return output.single_parameters

    # No parameter file, look for ski file in directory up
    else:

        simulation_path = fs.directory_of(output_path)
        filepaths = fs.files_in_path(simulation_path, extension="ski")
        #nfiles = len(filepaths)
        #if nfiles == 1 and fs.strip_extension(fs.name(filepaths[0])) == prefix: return filepaths[0]
        for filepath in filepaths:
            name = fs.strip_extension(fs.name(filepath))
            if name == prefix: return filepath

    # Nothing found
    return None

# -----------------------------------------------------------------

def load_wavelength_grid_in_cwd(input_path=None):

    """
    This function ...
    :param input_path:
    :return:
    """

    return load_wavelength_grid(fs.cwd(), input_path=input_path)

# -----------------------------------------------------------------

def load_wavelength_grid(output_path, input_path=None):

    """
    This function finds the wavelengths of a simulation from the output directory
    :param output_path:
    :param input_path:
    :return:
    """

    from ..data.sed import load_sed

    # Create simulation output
    output = SimulationOutput.from_directory(output_path)

    # Has wavelengths file in the output
    if output.has_single_wavelengths: return WavelengthGrid.from_skirt_output(output.single_wavelengths)

    # Has SED
    elif output.has_seds:

        sed = load_sed(output.seds[0])
        return WavelengthGrid.from_sed(sed)

    # Has parameters
    elif output.has_single_parameters:

        # Get ski
        parameters_path = output.single_parameters
        ski = SkiFile(parameters_path)
        return ski.get_wavelengths(input_path, as_grid=True)

    # Nothing
    else: raise IOError("No wavelength grid information can be found")

# -----------------------------------------------------------------

def load_wavelengths_in_cwd(unit="micron", input_path=None):

    """
    This function ...
    :param unit:
    :param input_path:
    :return:
    """

    return load_wavelengths(fs.cwd(), unit=unit, input_path=input_path)

# -----------------------------------------------------------------

def load_wavelengths(output_path, unit="micron", input_path=None):

    """
    This function ...
    :param output_path:
    :param unit:
    :param input_path:
    :return:
    """

    # Get the grid
    grid = load_wavelength_grid(output_path, input_path=input_path)

    # Return the wavelengths as an array
    return grid.wavelengths(asarray=True, unit=unit)

# -----------------------------------------------------------------

def load_seds_in_cwd(contributions=("total",), wavelength_unit=None, photometry_unit=None, instrument_names=None):

    """
    This function ...
    :param contributions:
    :param wavelength_unit:
    :param photometry_unit:
    :param instrument_names:
    :return:
    """

    return load_seds(fs.cwd(), contributions=contributions, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit, instrument_names=instrument_names)

# -----------------------------------------------------------------

def load_seds(output_path, contributions=("total",), wavelength_unit=None, photometry_unit=None, instrument_names=None):

    """
    This function ...
    :param output_path:
    :param contributions:
    :param wavelength_unit:
    :param photometry_unit:
    :param instrument_names:
    :return:
    """

    from ..data.sed import load_multiple_seds
    from .fluxes import get_sed_instrument_name

    # Find paths
    paths, prefix = find_seds(output_path, return_prefix=True)

    # Find instrument distances
    distances = load_instrument_distances(output_path, prefix=prefix)

    # Create dictionary for the SEDs
    seds = DefaultOrderedDict(OrderedDict)

    # Loop over the paths
    for path in paths:

        # Get instrument name
        instr_name = get_sed_instrument_name(path, prefix=prefix)
        if instrument_names is not None and instr_name not in instrument_names: continue

        # Load SEDs in this file
        seds_instrument = load_multiple_seds(path, as_dict=True, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)

        # Loop over the SEDs
        for contribution in seds_instrument:
            if contribution not in contributions: continue
            sed = seds_instrument[contribution]
            sed.distance = distances[instr_name]
            seds[instr_name][contribution] = sed

    # Return the SEDs
    return seds

# -----------------------------------------------------------------

def load_datacube_in_cwd(instrument_name, contribution="total", input_path=None):

    """
    This function ...
    :param instrument_name:
    :param contribution:
    :param input_path:
    :return:
    """

    return load_datacube(fs.cwd(), instrument_name, contribution=contribution, input_path=input_path)

# -----------------------------------------------------------------

def load_datacube(output_path, instrument_name, contribution="total", input_path=None):

    """
    This function ...
    :param output_path:
    :param instrument_name:
    :param contribution:
    :param input_path:
    :return:
    """

    datacubes = load_datacubes(output_path, contributions=(contribution,), instrument_names=(instrument_name,), input_path=input_path)
    return datacubes[instrument_name][contribution]

# -----------------------------------------------------------------

def load_datacubes_in_cwd(contributions=("total",), input_path=None, instrument_names=None):

    """
    Thisn function ...
    :param contributions:
    :param input_path:
    :param instrument_names:
    :return:
    """

    return load_datacubes(fs.cwd(), contributions=contributions, input_path=input_path, instrument_names=instrument_names)

# -----------------------------------------------------------------

def load_datacubes(output_path, contributions=("total",), input_path=None, instrument_names=None):

    """
    This function ...
    :param output_path:
    :param contributions:
    :param input_path:
    :param instrument_names:
    :return:
    """

    # Load datacube class
    from ...magic.core.datacube import DataCube

    # Find paths
    paths, prefix = find_datacubes(output_path, return_prefix=True, contributions=contributions)

    # Find wavelength grid
    wavelength_grid = load_wavelength_grid(output_path, input_path=input_path)

    # Find instrument distances
    distances = load_instrument_distances(output_path, prefix=prefix)
    #print(distances)

    # Create dictionary to contain the datacubes
    datacubes = DefaultOrderedDict(OrderedDict)

    # Loop over the paths
    for path in paths:

        # Get instrument name
        instr_name = get_datacube_instrument_name(path, prefix)
        if instrument_names is not None and instr_name not in instrument_names: continue

        # Get the contribution
        contribution = get_datacube_contribution(path)

        # Load the datacube
        datacube = DataCube.from_file(path, wavelength_grid)

        # Set the distance
        datacube.distance = distances[instr_name]

        # Add the datacube
        datacubes[instr_name][contribution] = datacube

    # Return the dictionary of datacubes
    return datacubes

# -----------------------------------------------------------------

class DatacubesMiscMaker(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(DatacubesMiscMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # The paths to the datacubes produced by SKIRT
        self.datacube_paths = None

        # The instrument names for which to make stuff
        self.instrument_names = None

        # The wavelength grid of the simulation
        self.wavelength_grid = None

        # The distances for the different instruments
        self.distances = None

    # -----------------------------------------------------------------

    @property
    def ndatacubes(self):

        """
        This function ...
        :return:
        """

        return len(self.datacube_paths)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DatacubesMiscMaker, self).setup(**kwargs)

        # Initialize (get the datacube paths and wavelengths)
        self.initialize(**kwargs)

        # Get output directory
        output_path = kwargs.pop("output_path", None)
        self.config.output = output_path

        # Get instrument (datacube names)
        self.get_instrument_names(**kwargs)

        # Get the distances
        self.get_distances(**kwargs)

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Initialize datacube paths, simulation prefix and wavelengths
        if "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.initialize_from_simulation(simulation)

        # Simulation directory is given
        elif "simulation_path" in kwargs:
            simulation = createsimulations(kwargs.pop("simulation_output_path"), single=True)
            self.initialize_from_simulation(simulation)

        # Simulation output directory is given
        elif "simulation_output_path" in kwargs:
            self.initialize_from_output_path(kwargs.pop("simulation_output_path"))

        # Multiple datacube paths are given
        elif "datacube_paths" in kwargs:

            # Check prefix
            if kwargs.get("prefix", None) is None: raise ValueError("Simulation prefix must be specified if datacube paths are passed explicitly")

            # Get datacube paths and prefix
            self.datacube_paths = kwargs.pop("datacube_paths")
            self.simulation_prefix = kwargs.pop("prefix")

            # Get wavelengths
            if kwargs.get("wavelengths", None) is not None: self.wavelengths = kwargs.pop("wavelengths")
            elif kwargs.get("wavelength_grid", None) is not None: self.wavelengths = kwargs.pop("wavelength_grid").wavelengths(asarray=True, unit="micron")
            elif kwargs.get("sed", None) is not None: self.wavelengths = kwargs.pop("sed").wavelengths(asarray=True, unit="micron")
            else:
                for datacube_path in self.datacube_paths:
                    if "total.fits" not in datacube_path: continue
                    sed_filepath = datacube_path.replace('total.fits', 'sed.dat')
                    if not fs.is_file(sed_filepath): continue
                    from ..data.sed import load_sed
                    sed = load_sed(sed_filepath)
                    wavelength_grid = WavelengthGrid.from_sed(sed)
                    self.wavelengths = wavelength_grid.wavelengths(asarray=True, unit="micron")
                    break
                else: raise IOError("Cannot find any SED file corresponding to the datacubes (for the wavelengths)") # no break encountered

        # One datacube path is given
        elif "datacube_path" in kwargs:

            # Check prefix
            if kwargs.get("prefix", None) is None: raise ValueError("Simulation prefix must be specified if datacube paths are passed explicitly")

            # Get datacube path and prefix
            filepath = kwargs.pop("datacube_path")
            self.simulation_prefix = kwargs.pop("prefix")
            self.datacube_paths = [filepath]

            # Get wavelengths
            if kwargs.get("wavelengths", None) is not None: self.wavelengths = kwargs.pop("wavelengths")
            elif kwargs.get("wavelength_grid", None) is not None: self.wavelengths = kwargs.pop("wavelength_grid").wavelengths(asarray=True, unit="micron")
            elif kwargs.get("sed", None) is not None: self.wavelengths = kwargs.pop("sed").wavelengths(asarray=True, unit="micron")
            else:
                sed_filepath = filepath.replace('total.fits', 'sed.dat')
                if not fs.is_file(sed_filepath): raise IOError("Cannot find corresponding SED file (for the wavelengths)")
                from ..data.sed import load_sed
                sed = load_sed(sed_filepath)
                wavelength_grid = WavelengthGrid.from_sed(sed)
                self.wavelengths = wavelength_grid.wavelengths(asarray=True, unit="micron")

        # Look in working directory
        else: self.initialize_from_cwd()

        # Check whether we have datacubes
        if self.datacube_paths is None or self.ndatacubes == 0: raise ValueError("No datacubes")

    # -----------------------------------------------------------------

    def initialize_from_simulation(self, simulation):

        """
        Thisf unction ...
        :param simulation:
        :return:
        """

        # Debugging
        log.debug("Initializing from simulation object ...")

        # Obtain the paths to the 'total' FITS files created by the simulation
        self.datacube_paths = simulation.totalfitspaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # If ski file is defined
        if simulation.has_ski:

            # Initialize dictionary
            self.distances = dict()

            # Load the ski file
            ski = simulation.parameters()

            # Loop over the different simulated TOTAL datacubes
            for path in self.total_datacube_paths:

                # Get the name of the instrument
                instr_name = get_datacube_instrument_name(path, self.simulation_prefix)

                # Get the distance
                distance = ski.get_instrument_distance(instr_name)

                # Set the distance
                self.distances[instr_name] = distance

    # -----------------------------------------------------------------

    def initialize_from_cwd(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Simulation or simulation output path not specified: searching for simulation output in the current working directory ...")

        # Initialize from current working directory
        self.initialize_from_output_path(self.config.path)

    # -----------------------------------------------------------------

    def initialize_from_output_path(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Debugging
        log.debug("Looking for total datacube files in directory '" + output_path + "' ...")

        # Set the paths to the total FITS files created by the simulation
        self.datacube_paths, self.simulation_prefix = find_datacubes(output_path, return_prefix=True)

        # Set the list of wavelengths for the simulation
        self.wavelengths = load_wavelengths(output_path)

    # -----------------------------------------------------------------

    def get_instrument_names(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the instrument names ...")

        # Instrument names
        if kwargs.get("instrument_names", None) is not None:

            # Check
            if "instruments" in kwargs: raise ValueError("Cannot specify 'instruments' and 'instrument_names' simultaneously")

            # Get names of the instruments (datacubes) of which to create observed images
            self.instrument_names = kwargs.pop("instrument_names")

        # Instruments
        elif kwargs.get("instruments", None) is not None: self.instrument_names = kwargs.pop("instruments")

        # From config
        elif self.config.instruments is not None: self.instrument_names = self.config.instruments

    # -----------------------------------------------------------------

    def get_distances(self, **kwargs):

        """
        This funciton ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the instrument distances ...")

        # TODO: check with distance found from ski file

        # Instrument distances
        if kwargs.get("distances", None) is not None:

            self.distances = kwargs.pop("distances")

        # Single distance for all instruments
        elif kwargs.get("distance", None) is not None:

            self.distances = dict()
            distance = kwargs.pop("distance")
            for instr_name in self.instrument_names:
                self.distances[instr_name] = distance

    # -----------------------------------------------------------------

    @lazyproperty
    def total_datacube_paths(self):

        """
        Thisn function ...
        :return:
        """

        paths = []
        for path in self.datacube_paths:
            if not is_total_datacube(path): continue
            paths.append(path)
        return paths

    # -----------------------------------------------------------------

    def make_for_instrument(self, instr_name):

        """
        This function ...
        :param instr_name:
        :return:
        """

        # If a list of instruments is defined an this instrument is not in this list, skip it
        if self.instrument_names is None: return True
        elif instr_name in self.instrument_names: return True
        else: return False

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")
        
        # Construct the wavelength grid from the array of wavelengths
        self.wavelength_grid = WavelengthGrid.from_wavelengths(self.wavelengths, "micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelengths)

# -----------------------------------------------------------------

def is_total_datacube(datacube_path):

    """
    This function ...
    :param datacube_path:
    :return:
    """

    return get_datacube_contribution(datacube_path) == "total"

# -----------------------------------------------------------------

def get_datacube_contribution(datacube_path):

    """
    This function ...
    :param datacube_path:
    :return:
    """

    name = fs.strip_extension(fs.name(datacube_path))
    return name.split("_")[-1]

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
