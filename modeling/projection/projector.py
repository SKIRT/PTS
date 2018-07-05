#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.project Contains the Projector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...modeling.basics.instruments import FrameInstrument
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.launch.launcher import SKIRTLauncher
from ...magic.core.frame import Frame
from ...modeling.basics.models import DeprojectionModel3D
from ...core.tools import introspection

# -----------------------------------------------------------------

map_filename = "map.fits"

# -----------------------------------------------------------------

# TODO: class should be connected to the ComponentProjections class in the other module

# -----------------------------------------------------------------

def project(model, projections, output_path=None, simulation_path=None, npackages=1e7, parallelization=None, coordinate_systems=None):

    """
    This function ...
    :param model:
    :param projections:
    :param output_path:
    :param simulation_path:
    :param npackages:
    :param parallelization
    :param coordinate_systems:
    :return:
    """

    # Create the projector
    projector = Projector()

    # Set output path
    if output_path is not None:
        projector.config.output = output_path
        projector.config.write = True
    else: projector.config.write = False

    # Set the simulation path
    projector.config.simulation_path = simulation_path

    # Set simulation options
    projector.config.npackages = npackages
    if parallelization is not None: projector.config.parallelization = parallelization

    # Run the projector
    projector.run(model=model, projections=projections, coordinate_systems=coordinate_systems)

    # Return the projected images
    return projector.projected

# -----------------------------------------------------------------

class Projector(Configurable):

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
        super(Projector, self).__init__(*args, **kwargs)

        # The model to project
        self.model = None

        # The projections
        self.projections = None

        # Coordinate systems
        self.coordinate_systems = None

        # The instruments
        self.instruments = dict()

        # The ski file
        self.ski = None

        # The simulation path and the simulation input path
        self.simulation_path = None
        self.in_path = None

        # The SKIRT launcher
        self.launcher = SKIRTLauncher()

        # The projected maps
        self.projected = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Create the instruments
        self.create_instruments()

        # 3. Set the input path
        if self.is_deprojection: self.set_input_path()

        # 4. Create the ski file
        self.create_ski()

        # 5. Write the ski file
        self.write_ski()

        # 6. Launch SKIRT
        self.launch()

        # 7. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Projector, self).setup(**kwargs)

        # Get the model
        self.model = kwargs.pop("model")

        # Get the projections
        self.projections = kwargs.pop("projections")

        # Coordinate systems?
        self.coordinate_systems = kwargs.pop("coordinate_systems", None)

        # Set the simulation path
        if self.config.simulation_path is not None: self.simulation_path = self.config.simulation_path
        else: self.simulation_path = introspection.create_unique_temp_dir("deprojection")

    # -----------------------------------------------------------------

    @property
    def has_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        return self.coordinate_systems is not None

    # -----------------------------------------------------------------

    def has_coordinate_system(self, instrument_name):

        """
        Thisn function ...
        :param instrument_name:
        :return:
        """

        return self.has_coordinate_systems and instrument_name in self.coordinate_systems

    # -----------------------------------------------------------------

    @lazyproperty
    def smile(self):

        """
        This function ...
        :return:
        """

        return SKIRTSmileSchema()

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the ski file template ...")

        # Create
        ski = self.smile.create_oligochromatic_template()

        # Remove the existing instruments
        ski.remove_all_instruments()

        # Remove the stellar system
        # ski.remove_stellar_system()

        # Remove the dust system
        if ski.has_dust_system: ski.remove_dust_system()

        # Set the number of photon packages
        ski.setpackages(self.config.npackages)

        # Enable writing options
        # ski.enable_all_writing_options()

        # Return the ski template
        return ski

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projections
        for name in self.projections:

            # Create frame instrument
            instrument = FrameInstrument.from_projection(self.projections[name])

            # Set the instrument
            self.instruments[name] = instrument

    # -----------------------------------------------------------------

    @property
    def is_deprojection(self):

        """
        This function ...
        :return:
        """

        return isinstance(self.model, DeprojectionModel3D)

    # -----------------------------------------------------------------

    def set_input_path(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the input path ...")

        # No map?
        if not self.model.has_map: raise ValueError("The map cannot be loaded")

        # .. because file is present?
        map_path = self.model.filepath
        if fs.is_file(map_path):

            self.in_path = fs.directory_of(map_path)  # directory of input map
            filename = fs.name(map_path)

        # .. or because map is loaded
        elif self.model.map_is_loaded:

            # Create input directory
            self.in_path = fs.create_directory_in(self.simulation_path, "in")

            # Save the map
            map_path = fs.join(self.in_path, map_filename)
            filename = map_filename
            self.model.map.saveto(map_path)  # save the map

        # We shouldn't get here
        else: raise RuntimeError("We shouldn't get here")

        # Set the model map filename
        self.model.filename = filename

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a ski file ...")

        # Make copy of ski file template
        ski = self.ski_template.copy()

        # Add the instruments
        for instrument_name in self.instruments: ski.add_instrument(instrument_name, self.instruments[instrument_name])

        # Remove the dust system
        if ski.has_dust_system: ski.remove_dust_system()

        # Add the stellar component
        ski.create_new_stellar_component(geometry=self.model, luminosities=[1])

        # Set the ski file
        self.ski = ski

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save the ski file
        self.ski.saveto(self.ski_path, fix=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_path(self):

        """
        Thisn function ...
        :return:
        """

        # Determine path
        return fs.join(self.simulation_path, "projection.ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.simulation_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This funtion ...
        :return:
        """

        # Create simulation definition
        return SingleSimulationDefinition(self.ski_path, self.out_path, self.in_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Set launcher settings
        self.launcher.config.show_progress = True

        # Run
        self.launcher.run(definition=self.definition, parallelization=self.config.parallelization)
        simulation = self.launcher.simulation

        # Loop over the instruments
        for instrument_name in self.instruments:

            # Determine the path to the simulated image file
            other_path = fs.join(self.out_path, simulation.prefix() + "_" + instrument_name + "_total.fits")
            other_map = Frame.from_file(other_path)

            # Get coordinate system for this map
            if self.has_coordinate_system(instrument_name): other_map.wcs = self.coordinate_systems[instrument_name]

            # Set the map
            self.projected[instrument_name] = other_map

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the projected maps
        self.write_projected()

    # -----------------------------------------------------------------

    def write_projected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projected maps ...")

        # Loop over the instruments
        for instrument_name in self.instruments:

            # Determine path
            path = self.output_path_file(instrument_name + ".fits")

            # Save
            self.projected[instrument_name].saveto(path)

# -----------------------------------------------------------------
