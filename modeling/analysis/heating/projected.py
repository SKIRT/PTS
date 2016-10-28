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
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....core.tools import tables
from ....core.data.sed import SED
from ....core.simulation.wavelengthgrid import WavelengthGrid
from ....magic.core.datacube import DataCube
from ....magic.plot.imagegrid import StandardImageGridPlotter

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

        self.total_output_path = None

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

        # The TIR maps
        self.total_tir_map = None
        self.evolved_tir_map = None
        self.unevolved_tir_map = None

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

        # 6. Calculate maps of the TIR luminosity
        self.calculate_tir_maps()

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

        # Determine the output path for the 'total' simulation
        self.total_output_path = self.analysis_run.heating_output_path_for_contribution("total")

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid file produced by SKIRT ...")

        # Determine the path to the wavelength grid file in heating/total
        wavelengths_path = fs.join(self.total_output_path, self.galaxy_name + "_wavelengths.dat")

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
        path = fs.join(self.total_output_path, self.galaxy_name + "_earth_sed.dat")

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
        path = fs.join(self.analysis_run.heating_output_path_for_contribution("old"), self.galaxy_name + "_earth_sed.dat")

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
        path = fs.join(self.analysis_run.heating_output_path_for_contribution("unevolved"), self.galaxy_name + "_earth_sed.dat")

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

        # Load total datacube
        self.load_total_datacube()

        # Load evolved datacube
        self.load_evolved_datacube()

        # Load unevolved datacube
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
        path = fs.join(self.total_output_path, self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.total_datacube = DataCube.from_file(path, self.wavelength_grid)

    # -----------------------------------------------------------------

    def load_evolved_datacube(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the datacube of the simulation with the evolved stellar population ...")

        # Determine the path to the datacube
        path = fs.join(self.analysis_run.heating_output_path_for_contribution("evolved"), self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.evolved_datacube = DataCube.from_file(path, self.wavelength_grid)

    # -----------------------------------------------------------------

    def load_unevolved_datacube(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the datacube of the simulation with the unevolved stellar populations ...")

        # Determine the path to the datacube
        path = fs.join(self.analysis_run.heating_output_path_for_contribution("unevolved"), self.galaxy_name + "_earth_total.fits")

        # Load the datacube
        self.unevolved_datacube = DataCube.from_file(path, self.wavelength_grid)

    # -----------------------------------------------------------------

    def calculate_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the heating fractions ...")

        # Calculate the flux fractions of the evolved and unevolved stellar populations over the wavelength spectrum
        self.calculate_sed_heating_fractions()

        # Calculate the flux fractions of the evolved and unevolved stellar population for selected wavelenghts in each pixel
        self.calculate_pixel_heating_fractions()

    # -----------------------------------------------------------------

    def calculate_sed_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the flux fractions over the wavelength spectrum ...")

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
        log.info("Calculating the flux fractions for selected wavelengths in each pixel ...")

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

    def calculate_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the TIR luminosity in each pixel ...")

        # Calculate the TIR maps
        self.calculate_total_tir_map()
        self.calculate_unevolved_tir_map()
        self.calculate_evolved_tir_map()

    # -----------------------------------------------------------------

    def calculate_total_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the total TIR luminosity in each pixel ...")

        cube = self.total_datacube.asarray()
        wavelengths = self.wavelength_grid.wavelengths(asarray=True, unit="micron")
        deltas = self.wavelength_grid.deltas(asarray=True, unit="micron")

        # Calculate the map
        self.total_tir_map = integrate_pixel_seds(cube, wavelengths, deltas)

    # -----------------------------------------------------------------

    def calculate_unevolved_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the unevolved TIR luminosity in each pixel ...")

        cube = self.unevolved_datacube.asarray()
        wavelengths = self.wavelength_grid.wavelengths(asarray=True, unit="micron")
        deltas = self.wavelength_grid.deltas(asarray=True, unit="micron")

        # Calculate the map
        self.unevolved_tir_map = integrate_pixel_seds(cube, wavelengths, deltas)

    # -----------------------------------------------------------------

    def calculate_evolved_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the evolved TIR luminosity in each pixel ...")

        cube = self.evolved_datacube.asarray()
        wavelengths = self.wavelength_grid.wavelengths(asarray=True, unit="micron")
        deltas = self.wavelength_grid.deltas(asarray=True, unit="micron")

        # Calculate the map
        self.evolved_tir_map = integrate_pixel_seds(cube, wavelengths, deltas)

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

        # Write the TIR maps
        self.write_tir_maps()

    # -----------------------------------------------------------------

    def write_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing a table with the flux fractions at different wavelengths ...")

        # Determine the path to the table of the flux fractions
        path = fs.join(self.projected_heating_path, "fractions.dat")

        # Write the table
        tables.write(self.fractions, path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def write_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the TIR maps ...")

        # Determine the path to the total TIR map
        path = fs.join(self.projected_heating_path, "tir_total.fits")

        # Write the total TIR map
        self.total_tir_map.save(path)

        # Determine the path to the unevolved TIR map
        path = fs.join(self.projected_heating_path, "tir_unevolved.fits")

        # Write the unevolved TIR map
        self.unevolved_tir_map.save(path)

        # Determine the path to the evolved TIR map
        path = fs.join(self.projected_heating_path, "tir_evolved.fits")

        # Write the evolved TIR map
        self.evolved_tir_map.save(path)

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
        path = fs.join(self.projected_heating_path, "fractions.pdf")

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

        # Create the image grid plotter
        plotter = StandardImageGridPlotter()

# -----------------------------------------------------------------

def integrate_pixel_seds(cube, wls, dwls):

    """
    This function ...
    :param cube:
    :param wls:
    :param dwls:
    :return:
    """

    Lsun = 3.846e26 # Watts
    MjySr_to_LsunMicron = 1.e6 * (36./206264.806247)**2 * 1.e-26 * 4*np.pi*(0.785e6*3.086e+16)**2 * 3.e14/(wls**2) / Lsun

    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])

    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = np.sum(sed * MjySr_to_LsunMicron * dwls)

    return slice

# -----------------------------------------------------------------
