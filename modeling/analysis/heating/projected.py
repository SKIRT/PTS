#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.projected Contains the ProjectedDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....core.tools import tables, inspection
from ...core.sed import SED
from ....core.simulation.wavelengthgrid import WavelengthGrid
from ....magic.core.image import Image

# -----------------------------------------------------------------

class ProjectedDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(ProjectedDustHeatingAnalyser, self).__init__(config)

        # -- Attributes --

        # The wavelength grid used for the simulations
        self.wavelength_grid = None

        # The number of wavelengths
        self.number_of_wavelengths = None

        # The SEDs
        self.total_sed = None
        self.evolved_sed = None
        self.unevolved_sed = None

        # The datacubes
        self.total_datacube = None
        self.evolved_datacube = None
        self.unevolved_datacube = None

        # The flux fractions
        self.fractions = None
        self.fraction_maps = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ProjectedDustHeatingAnalyser instance
        analyser = cls()

        # Set the modeling path
        analyser.config.path = arguments.path

        # Return the new instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the wavelength grid
        self.load_wavelength_grid()

        # 3. Load the simulated SEDs
        self.load_seds()

        # 4. Load the simulated data cubes
        self.load_datacubes()

        # 5. Calculate the heating fractions
        self.calculate_heating_fractions()

        # 6. Writing
        self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ProjectedDustHeatingAnalyser, self).setup()

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid file produced by SKIRT ...")

        # Determine the path to the wavelength grid file in heating/total
        wavelengths_path = fs.join(self.output_paths["total"], self.galaxy_name + "_wavelengths.dat")

        # Load the wavelength grid as a table
        self.wavelength_grid = WavelengthGrid.from_skirt_output(wavelengths_path)

        # Determine the number of wavelengths
        self.number_of_wavelengths = len(self.wavelength_grid)

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated SEDs ...")

        # Load the SED of the simulation with the total stellar population
        self.load_total_sed()

        # Load the SED of the simulation with the evolved stellar population
        self.load_evolved_sed()

        # Load the SED of the simulation with the unevolved stellar populations
        self.load_unevolved_sed()

    # -----------------------------------------------------------------

    def load_total_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SED of the simulation with the total stellar population ...")

        # Determine the path to the SED file
        path = fs.join(self.output_paths["total"], self.galaxy_name + "_earth_sed.dat")

        # Load the SED
        self.total_sed = SED.from_skirt(path)

    # -----------------------------------------------------------------

    def load_evolved_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SED of the simulation with the evolved stellar population ...")

        # Determine the path to the SED file
        path = fs.join(self.output_paths["old"], self.galaxy_name + "_earth_sed.dat")

        # Load the SED
        self.evolved_sed = SED.from_skirt(path)

    # -----------------------------------------------------------------

    def load_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SED of the simulation with the unevolved stellar populations ...")

        # Determine the path to the SED file
        path = fs.join(self.output_paths["unevolved"], self.galaxy_name + "_earth_sed.dat")

        # Load the SED
        self.unevolved_sed = SED.from_skirt(path)

    # -----------------------------------------------------------------

    def load_datacubes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated datacubes ...")

        #
        self.load_total_datacube()

        self.load_evolved_datacube()

        self.load_unevolved_datacube()

    # -----------------------------------------------------------------

    def load_total_datacube(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the datacube of the simulation with the total stellar population ...")

        # Determine the path to the datacube
        path = fs.join(self.output_paths["total"], self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.total_datacube = Image.from_file(path)

    # -----------------------------------------------------------------

    def load_evolved_datacube(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the datacube of the simulation with the evolved stellar population ...")

        # Determine the path to the datacube
        path = fs.join(self.output_paths["evolved"], self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.evolved_datacube = Image.from_file(path)

    # -----------------------------------------------------------------

    def load_unevolved_datacube(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the datacube of the simulation with the unevolved stellar populations ...")

        # Determine the path to the datacube
        path = fs.join(self.output_paths["unevolved"], self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.unevolved_datacube = Image.from_file(path)

    # -----------------------------------------------------------------

    def calculate_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the heating fractions ...")

        #
        self.calculate_sed_heating_fractions()

        #
        self.calculate_pixel_heating_fractions()

    # -----------------------------------------------------------------

    def calculate_sed_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # ...
        all_wavelengths = self.wavelength_grid.wavelengths(asarray=True, unit="micron")
        dust_wavelengths = all_wavelengths > 10.

        total_fluxes = self.total_sed.fluxes(asarray=True)
        evolved_fluxes = self.evolved_sed.fluxes(asarray=True)
        unevolved_fluxes = self.unevolved_sed.fluxes(asarray=True)

        # Calculate the heating fractions
        unevolved_fractions = 0.5 * (unevolved_fluxes + total_fluxes - evolved_fluxes) / total_fluxes
        evolved_fractions = 0.5 * (evolved_fluxes + total_fluxes - unevolved_fluxes) / total_fluxes

        # Create the table of flux fractions
        names = ["Wavelength", "Flux fraction from unevolved stellar populations",
                 "Flux fraction from evolved stellar population"]
        data = [all_wavelengths[dust_wavelengths], unevolved_fractions[dust_wavelengths],
                evolved_fractions[dust_wavelengths]]
        self.fractions = tables.new(data, names)

    # -----------------------------------------------------------------

    def calculate_pixel_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # The interesting wavelengths
        wavelengths = [12., 24., 70., 100., 250., 350., 500.]

        # Loop over the interesting wavelengths
        for wavelength in wavelengths:

            # Get the index of the wavelength grid point closest to this wavelength
            index = self.wavelength_grid.closest_wavelength_index(wavelength)
            wavelength = self.wavelength_grid.table["Wavelength"][index]

            # Determine the name of the corresponding frame in the datacube image
            frame_name = "frame" + str(index)

            total_fluxes = self.total_datacube.frames[frame_name]
            evolved_fluxes = self.evolved_datacube.frames[frame_name]
            unevolved_fluxes = self.unevolved_datacube.frames[frame_name]

            # Calculate the heating fractions
            unevolved_fractions = 0.5 * (unevolved_fluxes + total_fluxes - evolved_fluxes) / total_fluxes
            #evolved_fractions = 0.5 * (evolved_fluxes + total_fluxes - unevolved_fluxes) / total_fluxes

            # Add the fraction map to the dictionary
            self.fraction_maps[wavelength] = unevolved_fractions

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the table with the flux fractions
        self.write_fractions()

    # -----------------------------------------------------------------

    def write_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing a table with the flux fractions at different wavelengths ...")

        # Determine the path to the table of the flux fractions
        path = fs.join(self.analysis_heating_path, "fractions.dat")

        # Write the table
        tables.write(self.fractions, path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the flux fraction as a function of wavelength
        self.plot_fractions()

        # Plot the maps of the flux fraction at specific wavelengths
        self.plot_fraction_maps()

    # -----------------------------------------------------------------

    def plot_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the flux fractions as a function of wavelength ...")

        # Determine the path to the plot file
        path = fs.join(self.analysis_heating_path, "fractions.pdf")

        # Create the figure
        plt.figure(figsize=(7, 5))
        plt.ylabel('$F^\prime_{\lambda,\mathrm{unev.}}$ [$\%$]', fontsize=20)
        plt.xlabel('$\lambda/\mu\mathrm{m}$', fontsize=20)

        plt.xlim(10., 1.e3)
        plt.ylim(0., 60.)

        plt.xscale('log')
        plt.tick_params(labelsize=20)

        # plt.subplots_adjust(bottom=0.1)
        # plt.plot(dustwls,Fold, 'r-', label="old stars")

        plt.plot(self.fractions["Wavelength"], self.fractions["Flux fraction from unevolved stellar populations"], 'k-', label="Unevolved stellar populations")
        plt.plot(self.fractions["Wavelength"], self.fractions["Flux fraction from evolved stellar population"], "g-", label="Evolved stellar population")

        # plt.plot(dustwls,Fyoung_alternative1, 'r-', label="alt 1")
        # plt.plot(dustwls,Fyoung_alternative2, 'g-', label="alt 2")
        # plt.plot(dustwls,Fyoung_alternative3, 'c-', label="alt 3")

        #plt.plot(dustwls, Fyoung_alternative4, 'k-', label="alt 4")

        #plt.fill_between(dustwls, Fyoung, Fyoung_alternative4, color='grey', alpha='0.5')

        plt.tight_layout()
        # plt.legend(loc='upper left',numpoints=1,markerscale=1.5,fontsize=14)

        # Save the figure
        plt.savefig(path)
        plt.close()

    # -----------------------------------------------------------------

    def plot_fraction_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the heating fraction by unevolved stars at specific wavelengths ...")

        fig_x_size =
        fig_y_sie =

        # Create a figure
        self._figure = plt.figure(figsize=(fig_x_size, fig_y_size))
        self._figure.subplots_adjust(left=0.05, right=0.95)

        # Creat grid
        self._grid = AxesGrid(self._figure, 111,
                              nrows_ncols=(len(self.rows), 3),
                              axes_pad=0.02,
                              label_mode="L",
                              share_all=True,
                              cbar_location="right",
                              cbar_mode="single",
                              cbar_size="0.5%",
                              cbar_pad="0.5%",
                              )  # cbar_mode="single"

        for cax in self._grid.cbar_axes:
            cax.toggle_label(False)

# -----------------------------------------------------------------
