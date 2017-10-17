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

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..basics.log import log
from ..simulation.simulation import createsimulations
from ..tools import filesystem as fs
from ..tools.utils import lazyproperty
from ..simulation.wavelengthgrid import WavelengthGrid

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

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DatacubesMiscMaker, self).setup(**kwargs)

        # Initialize datacube paths, simulation prefix and wavelengths
        if "simulation" in kwargs:
            simulation = kwargs.pop("simulation")
            self.initialize_from_simulation(simulation)
        elif "simulation_path" in kwargs:
            simulation = createsimulations(kwargs.pop("simulation_output_path"), single=True)
            self.initialize_from_simulation(simulation)
        elif "simulation_output_path" in kwargs:
            self.initialize_from_output_path(kwargs.pop("simulation_output_path"))
        else: self.initialize_from_cwd()

        # Get output directory
        output_path = kwargs.pop("output_path", None)
        self.config.output = output_path

        # Get instrument (datacube names)
        self.get_instrument_names(**kwargs)

        # Get the distances
        self.get_distances(**kwargs)

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
                instr_name = get_instrument_name(path, self.simulation_prefix)

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

        # Get datacube paths
        #datacube_paths = fs.files_in_path(self.config.path, extension="fits")
        total_datacube_paths = fs.files_in_path(output_path, extension="fits", endswith="_total")

        # Get SED paths
        sed_paths = fs.files_in_path(output_path, extension="dat", endswith="_sed")

        # Determine prefix
        prefix = None
        for path in total_datacube_paths:
            filename = fs.strip_extension(fs.name(path))
            if prefix is None:
                prefix = filename.split("_")[0]
            elif prefix != filename.split("_")[0]:
                raise IOError("Not all datacubes have the same simulation prefix")
        if prefix is None: raise IOError("No datacubes were found")

        # Set the paths to the toatl FITS files created by the simulation
        self.datacube_paths = total_datacube_paths

        from ..data.sed import load_sed

        # Load one of the SEDs
        if len(sed_paths) == 0: raise IOError("No SED files") # TODO: look for wavelength grid
        sed = load_sed(sed_paths[0])

        # Set the list of wavelengths for the simulation
        self.wavelengths = sed.wavelengths(asarray=True, unit="micron")

        # Set the simulation prefix
        self.simulation_prefix = prefix

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

    name = fs.strip_extension(fs.name(datacube_path))
    if name.endswith("_total"): return True
    else: return False

# -----------------------------------------------------------------

def get_instrument_name(datacube_path, prefix):

    """
    This function ...
    :param datacube_path:
    :param prefix:
    :return:
    """

    # ONLY FOR TOTAL
    #return fs.name(datacube_path).split("_total.fits")[0].split(prefix + "_")[1]

    # For all
    return fs.name(datacube_path).split(prefix + "_")[1].rsplit("_", 1)[0]

# -----------------------------------------------------------------
