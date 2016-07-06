#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.bestmodellauncher Contains the BestModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy import constants

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import inspection, tables
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.basics.filter import Filter
from ..basics.models import SersicModel, DeprojectionModel
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..decomposition.decomposition import load_parameters
from ...magic.basics.skyregion import SkyRegion
from ..basics.instruments import SEDInstrument, FrameInstrument
from ..core.sun import Sun
from ..core.mappings import Mappings
from ...magic.tools import wavelengths
from ...core.tools.logging import log
from ..basics.projection import GalaxyProjection
from ..core.sed import ObservedSED
from .wavelengthgrids import WavelengthGridGenerator
from .dustgrids import DustGridGenerator
from ...core.basics.range import IntegerRange, RealRange, QuantityRange

# -----------------------------------------------------------------

template_ski_path = fs.join(inspection.pts_dat_dir("modeling"), "ski", "template.ski")

# -----------------------------------------------------------------

class BestModelLauncher(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(BestModelLauncher, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski = None

        # The structural parameters
        self.parameters = None

        # The projection system
        self.projection = None

        # The truncation ellipse
        self.ellipse = None

        # The geometric bulge model
        self.bulge = None

        # The deprojection model
        self.deprojection = None
        self.deprojections = dict()

        # The instrument
        self.instrument = None

        # The table of weights for each band
        self.weights = None

        # The observed SED
        self.observed_sed = None

        # Filters
        self.i1 = None
        self.fuv = None

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

        # Coordinate system
        self.reference_wcs = None

        # The ski files for simulating the contributions of the various stellar components
        self.ski_contributions = dict()

        # The ski file for generating simulated images
        self.ski_images = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.load_input()

        # 3. Create the wavelength grid
        self.create_wavelength_grids()

        # 4. Create the bulge model
        self.create_bulge_model()

        # 5. Create the deprojection model
        self.create_deprojection_model()

        # 6. Create the instrument
        self.create_instrument()

        # 7. Create the dust grids
        self.create_dust_grids()

        # 8. Adjust the ski file
        self.adjust_ski()

        # 9. Adjust the ski files for simulating the contributions of the various stellar components
        self.adjust_ski_contributions()

        # 10. Adjust the ski file for generating simulated images
        self.adjust_ski_images()

        # 11. Calculate the weight factor to give to each band
        self.calculate_weights()

        # 12. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingInitializer, self).setup()

        # Create filters
        self.i1 = Filter.from_string("I1")
        self.fuv = Filter.from_string("FUV")

        # Solar properties
        sun = Sun()
        self.sun_fuv = sun.luminosity_for_filter_as_unit(self.fuv) # Get the luminosity of the Sun in the FUV band
        self.sun_i1 = sun.luminosity_for_filter_as_unit(self.i1)   # Get the luminosity of the Sun in the IRAC I1 band

        # Reference coordinate system
        reference_path = fs.join(self.truncation_path, self.reference_image + ".fits")
        self.reference_wcs = CoordinateSystem.from_file(reference_path)

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # 1. Load the template ski file
        self.load_template()

        # 2. Load the structural parameters of the galaxy
        self.load_parameters()

        # 3. Load the projection system
        self.load_projection()

        # 4. Load the truncation ellipse
        self.load_truncation_ellipse()

        # 5. Load the observed SED
        self.load_observed_sed()

    # -----------------------------------------------------------------

    
    # -----------------------------------------------------------------

    def adjust_ski_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for simulating the contribution of the various stellar components ...")

        # Loop over the different contributions, create seperate ski file instance
        contributions = ["old", "young", "ionizing"]
        component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                           "young": "Young stars",
                           "ionizing": "Ionizing stars"}
        for contribution in contributions:

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # Remove other stellar components
            ski.remove_stellar_components_except(component_names[contribution])

            # Add the ski file to the dictionary
            self.ski_contributions[contribution] = ski

    # -----------------------------------------------------------------

    def adjust_ski_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for generating simulated images ...")

        # Create a copy of the ski file instance
        self.ski_images = self.ski.copy()

        # Remove all instruments
        self.ski_images.remove_all_instruments()

        # Create frame instrument to generate datacube
        frame_instrument = FrameInstrument.from_projection(self.projection)

        # Add the frame instrument
        self.ski_images.add_instrument("earth", frame_instrument)

        # Add the SED instrument
        self.ski_images.add_instrument("earth", self.instrument)
        
    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski files for simulating the contributions of the various stellar components
        self.write_ski_files_contributions()

        # Write the ski file for generating simulated images
        self.write_ski_file_images()

    # -----------------------------------------------------------------

    def write_ski_files_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files for simulating the contribution of the various stellar components ...")

        # Loop over the ski files
        for contribution in self.ski_contributions:

            # Determine the path to the ski file
            ski_path = fs.join(self.fit_best_contribution_paths[contribution], self.galaxy_name + ".ski")

            # Write the ski file
            self.ski_contributions[contribution].saveto(ski_path)

    # -----------------------------------------------------------------

    def write_ski_file_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file for creating simulated images ...")

        # Determine the path to the ski file
        ski_path = fs.join(self.fit_best_images_path, self.galaxy_name + ".ski")

        # Write the ski file
        self.ski_images.saveto(ski_path)

# -----------------------------------------------------------------
