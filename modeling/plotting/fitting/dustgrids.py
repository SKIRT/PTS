#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.chisquared Contains the ChiSquaredPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.plot.grids import plotgrids
from ....core.simulation.simulation import SkirtSimulation

# -----------------------------------------------------------------

class DustGridsPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DustGridsPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust grids ...")

        # Loop over the different dust grids
        index = 0
        for path in fs.directories_in_path(self.dust_grids_path, sort=int):

            ski_path = fs.join(path, self.object_name + ".ski")
            simulation = SkirtSimulation(ski_path=ski_path, outpath=path)

            # Plot the grid
            plotgrids(simulation, output_path=self.plot_fitting_dust_grids_path, silent=True, prefix=str(index))

            # Increment the index
            index += 1

# -----------------------------------------------------------------
