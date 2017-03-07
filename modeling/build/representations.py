#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations Contains the RepresentationsGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ...core.tools.logging import log
from ...core.prep.wavelengthgrids import WavelengthGridGenerator
from ...core.prep.dustgrids import DustGridGenerator
from ...core.basics.range import RealRange, QuantityRange
from ...core.basics.unit import parse_unit as u
from ..build.component import get_stellar_component_names, get_dust_component_names, load_stellar_component, load_dust_component
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.basics.configuration import prompt_string
from ...core.basics.quantity import represent_quantity

# -----------------------------------------------------------------

class RepresentationGenerator(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(RepresentationBuilder, self).__init__(config, interactive)
        
        # The wavelength grid and dust grid generators
        self.wg_generator = None
        self.dg_generator = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 3. Create the wavelength grids
        self.create_wavelength_grids()
        
        # 6. Create the dust grids
        self.create_dust_grids()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RepresentationBuilder, self).setup(**kwargs)

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

# -----------------------------------------------------------------
