#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .tables import GenerationsTable, ChiSquaredTable, ParametersTable, BestParametersTable
from ...core.simulation.skifile import LabeledSkiFile
from ...core.basics.distribution import Distribution
from ..basics.instruments import load_instrument
from ..core.model import Model
from ...core.simulation.grids import load_grid
from ...core.simulation.skifile import SkiFile
from ...core.simulation.simulation import SkirtSimulation
from .tables import ModelProbabilitiesTable

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        self.runs_table_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------



# -----------------------------------------------------------------

def get_run_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fit_path = fs.join(modeling_path, "fit")
    return fs.directories_in_path(fit_path, returns="name")

# -----------------------------------------------------------------
