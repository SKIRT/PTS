#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.simulation.sed Contains the ComponentSED class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import StaticSkirtSimulation
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import filesystem as fs
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.basics.log import log
from ...core.simulation.execute import run_simulation
from ...core.tools import introspection
from ...core.data.sed import SED
from ...core.simulation.skifile import SkiFile
from ..basics.instruments import SEDInstrument
from ..basics.models import DeprojectionModel3D
from ...core.prep.smile import get_panchromatic_template
from ..build.construct import add_new_stellar_component
from ...core.basics import containers
from ...core.simulation.wavelengthgrid import load_wavelength_grid

# -----------------------------------------------------------------

default_npackages = 1e5

# Instruments/orientations
earth_name = "earth"

# Wavelength file name
default_wavelengths_filename = "wavelengths.txt"

# -----------------------------------------------------------------

class ComponentSED(object):

    """
    This class ...
    """

    def __init__(self, name, component, wavelength_grid=None, description=None, path=None, npackages=default_npackages,
                 input_filepaths=None, distance=None, inclination=None, position_angle=None, wavelengths_filename=None):

        """
        The constructor ...
        :param name:
        :param component:
        :param wavelength_grid:
        :param description:
        :param path:
        :param npackages:
        :param input_filepaths:
        :param distance:
        :param inclination:
        :param position_angle:
        :param wavelengths_filename:
        """

        # Set the name and description
        self.name = name
        self.description = description

        # Set the path
        if path is None: path = introspection.create_temp_dir("sed__" + name)
        self.path = path

        # Set the wavelength grid
        if wavelength_grid is not None:

            if wavelengths_filename is not None: raise ValueError("Cannot specify wavelength grid and wavelength grid filename")
            self.wavelength_grid = wavelength_grid

            # Save the wavelengths file so SKIRT can use it
            wavelengths_filepath = fs.join(self.path, default_wavelengths_filename)
            self.wavelength_grid.saveto(wavelengths_filepath)

            # Add to input filepaths
            input_filepaths[default_wavelengths_filename] = wavelengths_filepath

        # Filename is given
        elif wavelengths_filename is not None:

            # Is actually a file path
            #if fs.is_absolute(wavelengths_filename) and fs.is_file(wavelengths_filename):
            if fs.is_file(wavelengths_filename):

                # Check if is input filepaths
                if input_filepaths is None:
                    input_filepaths = {default_wavelengths_filename: wavelengths_filename}
                    self.wavelengths_filename = default_wavelengths_filename
                elif wavelengths_filename not in input_filepaths.values():
                    input_filepaths[default_wavelengths_filename] = fs.absolute_path(wavelengths_filename)
                    self.wavelengths_filename = default_wavelengths_filename
                else: self.wavelengths_filename = containers.get_key_for_value(input_filepaths, wavelengths_filename)

                # Load
                self.wavelength_grid = load_wavelength_grid(wavelengths_filename)

            # Is a filename
            else:

                if input_filepaths is None: raise ValueError("Input filepaths are not defined, cannot locate wavelength file")
                if wavelengths_filename not in input_filepaths: raise ValueError("Path of the wavelength file is not defined")
                self.wavelengths_filename = wavelengths_filename

                # Load
                self.wavelength_grid = load_wavelength_grid(input_filepaths[self.wavelengths_filename])

        # No wavelength grid (file) specified
        else: raise ValueError("Must specify either wavelength grid or wavelengths filename")

        # Set properties
        if distance is not None: self.distance = distance
        if inclination is not None: self.inclination = inclination
        if position_angle is not None: self.position_angle = position_angle

        # Set the component
        self.component = component

        # Set the number of photon packages
        self.npackages = npackages

        # Set the input filepaths
        self.input_paths = input_filepaths

        # Run the simulation?
        self.simulation = self.get_simulation()

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):
        return self.name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):
        return self.name

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grid(self):
        return self.wavelength_grid is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def out_path(self):
        return fs.create_directory_in(self.path, "out")

    # -----------------------------------------------------------------

    @property
    def has_sed(self):
        return fs.has_files_in_path(self.out_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @property
    def has_simulation(self):
        return self.has_sed

    # -----------------------------------------------------------------

    def get_simulation(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        if self.has_simulation: return StaticSkirtSimulation(prefix=self.simulation_prefix, inpath=self.input_paths, outpath=self.out_path, ski_path=self.ski_path, name=self.simulation_name)

        # Run the simulation
        return self.run_simulation()

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_skifile: self.create_skifile()

        # Create the definition and return
        return SingleSimulationDefinition(self.ski_path, self.out_path, input_path=self.input_paths, name=self.name)

    # -----------------------------------------------------------------

    @property
    def has_skifile(self):
        return fs.is_file(self.ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_path(self):
        return fs.join(self.path, self.name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def skifile(self):

        """
        This function ...
        :return:
        """

        # Load the ski file if it already exists
        if self.has_skifile: return SkiFile(self.ski_path)

        # Create
        else: return self.create_skifile()

    # -----------------------------------------------------------------

    @property
    def has_component_model(self):
        return "model" in self.component and self.component.model is not None

    # -----------------------------------------------------------------

    @property
    def has_component_deprojection(self):
        return "deprojection" in self.component and self.component.deprojection is not None

    # -----------------------------------------------------------------

    @property
    def has_model(self):
        return self.has_component_model or self.has_component_deprojection

    # -----------------------------------------------------------------

    @property
    def model(self):
        if self.has_component_model: return self.component.model
        elif self.has_component_deprojection: return self.component.deprojection
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_deprojection(self):
        return self.has_model and isinstance(self.model, DeprojectionModel3D)

    # -----------------------------------------------------------------

    @property
    def deprojection(self):
        if self.has_deprojection: return self.model
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def distance(self):
        if not self.has_deprojection: return None
        else: return self.deprojection.distance

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def inclination(self):
        if not self.has_deprojection: return None
        else: return self.deprojection.inclination

    # -----------------------------------------------------------------

    @property
    def has_inclination(self):
        return self.inclination is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def position_angle(self):
        if not self.has_deprojection: return None
        else: return self.deprojection.position_angle

    # -----------------------------------------------------------------

    @property
    def has_position_angle(self):
        return self.position_angle is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):

        """
        This function ...
        :return:
        """

        # Check properties
        if not self.has_distance: raise ValueError("Distance is not defined")
        if not self.has_inclination:
            log.warning("Inclination is not defined, assuming isotropic model ...")
            inclination = Angle(0.0, "deg")
        else: inclination = self.inclination
        if not self.has_position_angle:
            log.warning("Position angle is not defined, assuming isotropic model ...")
            position_angle = Angle(0.0, "deg")
        else: position_angle = self.position_angle

        # Create and return the instrument
        return SEDInstrument.from_properties(self.distance, inclination, position_angle)

    # -----------------------------------------------------------------

    def create_skifile(self, npackages=default_npackages):

        """
        This function ...
        :param npackages:
        :return:
        """

        # Check whether the wavelength grid is defined
        if not self.has_wavelength_grid: raise ValueError("Wavelength grid path must be set")

        # Create a ski template
        ski = get_panchromatic_template()

        # Add the old stellar bulge component
        add_new_stellar_component(ski, self.name, self.component)

        # Add the instrument
        ski.add_instrument(earth_name, self.sed_instrument)

        # Set the wavelength grid
        ski.set_file_wavelength_grid(self.wavelengths_filename)

        # Set the number of photon packages
        ski.setpackages(npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.ski_path, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty
    def description_string(self):

        """
        This function ...
        :return:
        """

        if self.description is None: return ""
        else: return self.description + " "

    # -----------------------------------------------------------------

    def run_simulation(self):

        """
        This function ...
        :return:
        """

        # Show message
        log.info("Running SKIRT for the " + self.description_string + "intrinsic SED ...")

        # Run simulation
        return run_simulation(self.definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def output(self):
        return self.simulation.output

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_filepath(self):
        return self.output.single_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sed(self):
        return SED.from_skirt(self.sed_filepath)

# -----------------------------------------------------------------
