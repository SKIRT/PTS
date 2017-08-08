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
from ...misc.geometryplotter import GeometryPlotter
from ...basics.models import load_3d_model

# -----------------------------------------------------------------

class GeometriesPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GeometriesPlotter, self).__init__(*args, **kwargs)

        # The geometries
        self.geometries = dict()

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the geometries ...")

        # Load all geometries in the directory
        for path, name in fs.files_in_path(self.fit_geometries_path, extension="mod", returns=["path", "name"]):

            # Load the geometry
            model = load_3d_model(path)

            # Add the geometry
            self.geometries[name] = model

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the geometries of the model ...")

        # Create the geometry plotter
        plotter = GeometryPlotter()

        # Add the geometries
        for label in self.geometries: plotter.add_geometry(self.geometries[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "geometries.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
