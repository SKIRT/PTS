#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.imagesrepresentation Contains the ImagesRepresentationBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.basics.configurable import Configurable
from ...core.prep.smile import SKIRTSmileSchema
from ...core.launch.launcher import SKIRTLauncher
from .construct import add_dust_component
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.prep.smile import expected_types

# -----------------------------------------------------------------

simulation_prefix = "dustgrid"
skifilename = simulation_prefix + ".ski"

# -----------------------------------------------------------------

class DustGridBuilder(Configurable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(DustGridBuilder, self).__init__(*args, **kwargs)

        # Smile
        self.smile = SKIRTSmileSchema()

        # Create the SKIRT launcher
        self.launcher = SKIRTLauncher()

        # The ski file template
        self.ski = None

        # The model
        self.definition = None

        # THe model representation
        self.representation = None

        # The dictionary of input map paths
        self.input_map_paths = dict()

        # Paths
        self.ski_path = None
        self.out_path = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the ski file
        self.create_ski()

        # 7. Writing
        self.write()

        # Lauch the simulation
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustGridBuilder, self).setup(**kwargs)

        # Determine path
        self.ski_path = self.output_path_file(skifilename)

        # Determine output path
        self.out_path = self.output_path #self.output_path_directory("out")

        # Get model definition and representation
        self.definition = kwargs.pop("definition")
        self.representation = kwargs.pop("representation")

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Create template ski file
        #ski = self.smile.create_panchromatic_template()
        self.ski = self.smile.create_oligochromatic_template()

        # Set components
        self.set_components()

        #self.set_instruments()

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # 1. Set stellar components
        #self.set_stellar_components()

        # 2. Set dust components
        self.set_dust_components()

    # -----------------------------------------------------------------

    def set_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar components ...")

    # -----------------------------------------------------------------

    def set_dust_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust components ...")

        # Loop over the dust components
        for name in self.definition.dust_component_names:

            # Load the component
            component = self.definition.get_dust_component(name)

            #print(name, component)

            # Add the dust component
            map_filename = add_dust_component(self.ski, name, component)

            # If map filename is defined, set path in dictionary
            if map_filename is not None: self.input_map_paths[map_filename] = component.map_path

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instrument
        #self.ski.add_instrument("earth", self.representation.sed_instrument)

        # Set the number of photon packages
        self.ski.setpackages(0)

        # Set the name of the wavelength grid file
        #self.ski.set_file_wavelength_grid("wavelengths.txt")
        # NO -> OLIGO

        # Set the dust emissivityex
        #self.set_dust_emissivity()

        # Set the dust grid
        self.ski.set_dust_grid(self.representation.dust_grid)

        # Set all-cells dust library
        #self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        #self.set_selfabsorption()

        # Disable all writing options
        #self.ski.disable_all_writing_options()

        # Enable writing options
        self.ski.enable_all_writing_options()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save
        self.ski.saveto(self.ski_path)

        from ...core.tools import filesystem as fs
        #fs.open_file(self.ski_path)

        replacements = dict()

        for a, b in expected_types.items():
            from_string = a + ' type=""'
            to_string = a + ' type="' + b + '"'
            replacements[from_string] = to_string

        fs.replace_strings(self.ski_path, replacements)

        fs.open_file(self.ski_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Create definition
        definition = SingleSimulationDefinition(self.ski_path, self.out_path, self.input_map_paths)

        # Run
        self.launcher.run(definition=definition)

# -----------------------------------------------------------------
