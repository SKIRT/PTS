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
from ....core.simulation.wavelengthgrid import WavelengthGrid
from ....core.data.transmission import TransmissionCurve
from ....core.plot.transmission import TransmissionPlotter
from ....core.plot.wavelengthgrid import WavelengthGridPlotter
from ....core.basics.emissionlines import EmissionLines
from ....core.data.seds import load_example_mappings_sed, load_example_bruzualcharlot_sed, load_example_zubko_sed
from ....core.data.sed import SED, ObservedSED

# -----------------------------------------------------------------

class WavelengthGridsPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(WavelengthGridsPlotter, self).__init__(*args, **kwargs)

        # The wavelength grids
        self.wavelength_grids = []

        # The transmission curves
        self.transmission_curves = dict()

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Load the transmission curves
        self.load_transmission_curves()

        # load wavelengthg grids
        self.load_wavelength_grids()

    # -----------------------------------------------------------------

    def load_transmission_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the transmission curves for the observation filters ...")

        # Load the observed SED
        sed = ObservedSED.from_file(self.observed_sed_path)
        for fltr in sed.filters():

            # Create the transmission curve
            transmission = TransmissionCurve.from_filter(fltr)

            # Normalize the transmission curve
            transmission.normalize(value=1.0, method="max")

            # Add the transmission curve to the dictionary
            self.transmission_curves[str(fltr)] = transmission

    # -----------------------------------------------------------------

    def load_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grids used for the fitting ...")

        # Loop over the files found in the fit/wavelength grids directory
        for path in fs.files_in_path(self.fitting_run.fit_wavelength_grids_path, extension="txt", sort=int, not_contains="grids"):

            # Load the wavelength grid
            grid = WavelengthGrid.from_skirt_input(path)

            # Add the grid
            self.wavelength_grids.append(grid)

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        # Loop over the different wavelength grids
        index = 0
        for grid in self.wavelength_grids:
            # Plot the low-resolution and high-resolution wavelength grids with the observation filters
            path = fs.join(self.plot_fitting_wavelength_grids_path, "wavelengths_" + str(index) + "_filters.pdf")
            self.plot_wavelengths_filters(grid, path)

            # Plot the wavelengths with the SEDs of stars and dust
            path = fs.join(self.plot_fitting_wavelength_grids_path, "wavelengths_" + str(index) + "_seds.pdf")
            self.plot_wavelengths_seds(grid, path)

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def plot_wavelengths_filters(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with")

        # Create the transmission plotter
        plotter = TransmissionPlotter()

        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the transmission curves
        for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label],
                                                                              label)

        # Add the wavelength points
        for wavelength in grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength,
                    min_transmission=0.0, max_transmission=1.05)

    # -----------------------------------------------------------------

    def plot_wavelengths_seds(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with SEDs of stars and dust ...")

        # Create the wavelength grid plotter
        plotter = WavelengthGridPlotter()

        # Set title
        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the wavelength grid
        plotter.add_wavelength_grid(grid, "fitting simulations")

        # Add MAPPINGS SFR SED
        mappings_sed = load_example_mappings_sed()
        plotter.add_sed(mappings_sed, "MAPPINGS")

        # Add Bruzual-Charlot stellar SED
        bc_sed = load_example_bruzualcharlot_sed()
        plotter.add_sed(bc_sed, "Bruzual-Charlot")

        # Add Zubko dust emission SED
        zubko_sed = load_example_zubko_sed()
        plotter.add_sed(zubko_sed, "Zubko")

        # Add emission lines
        emission_lines = EmissionLines()
        for line in emission_lines: plotter.add_emission_line(line)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength)

# -----------------------------------------------------------------
