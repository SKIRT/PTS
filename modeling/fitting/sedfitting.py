#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.sedfitting Contains the SEDFitter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import imageio
from matplotlib import pyplot as plt

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import tables, filesystem
from ...core.tools.logging import log
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

descriptions = {"FUV young": "FUV luminosity of the young stars",
                "FUV ionizing": "FUV luminosity of the ionizing stars",
                "Dust mass": "dust mass"}

# -----------------------------------------------------------------

class SEDFitter(FittingComponent):
    
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
        super(SEDFitter, self).__init__(config)

        # -- Attributes --

        # The table of chi squared values
        self.chi_squared = None

        # The model parameter table
        self.parameters = None

        # The best parameter values
        self.best_fuv_young = None
        self.best_fuv_ionizing = None
        self.best_dust_mass = None

        # The tables with the probability distributions for the different fit parameters
        self.distributions = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SEDFitter instance
        fitter = cls(arguments.config)

        # Set the modeling path
        fitter.config.path = arguments.path

        # Return the new instance
        return fitter

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the chi squared table
        self.load_chi_squared()

        # 3. Load the parameter table
        self.load_parameters()

        # 4. Calculate the probability distributions
        self.calculate_distributions()

        # 5. Make an animation of the fitting procedure
        self.animate()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDFitter, self).setup()

    # -----------------------------------------------------------------

    def load_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the table with the chi squared value for each model ...")

        # Load the chi squared table
        self.chi_squared = tables.from_file(self.chi_squared_table_path, format="ascii.ecsv")

        # Check whether the table is non-empty
        if len(self.chi_squared) == 0: raise RuntimeError("Could not find any chi squared value, it appears no simulations have been run yet")

        # Sort the table for decreasing chi squared value
        self.chi_squared.sort("Chi squared")
        self.chi_squared.reverse()

    # -----------------------------------------------------------------

    def load_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the table with the model parameters ...")

        # Load the parameter table
        self.parameters = tables.from_file(self.parameter_table_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def calculate_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probability distributions for each parameter ...")

        # Get the chi squared values
        chi_squared_values = self.chi_squared["Chi squared"]

        # Calculate the probability for each model
        probabilities = np.exp(-0.5 * chi_squared_values)

        # Initialize lists to contain the parameter values for each simulation
        young_lum_of_simulations = []
        ionizing_lum_of_simulations = []
        dust_mass_of_simulations = []

        # From the parameters table, get the entries that correspond to simulations that
        # are finished (and for which a chi squared value has been calculated)
        for i in range(len(self.chi_squared)):

            # Get the name of the simulation
            simulation_name = self.chi_squared["Simulation name"][i]

            # Find the index of this simulation in the parameters table
            j = tables.find_index(self.parameters, simulation_name)

            # "FUV young" "FUV ionizing" "Dust mass"
            young_lum_of_simulations.append(self.parameters["FUV young"][j])
            ionizing_lum_of_simulations.append(self.parameters["FUV ionizing"][j])
            dust_mass_of_simulations.append(self.parameters["Dust mass"][j])

        # Convert into NumPy arrays
        young_lum_of_simulations = np.array(young_lum_of_simulations)
        ionizing_lum_of_simulations = np.array(ionizing_lum_of_simulations)
        dust_mass_of_simulations = np.array(dust_mass_of_simulations)

        # FUV luminosity of young stellar population
        young_lum_unique = sorted(list(set(young_lum_of_simulations)))

        # FUV luminosity of ionizing stellar population
        ionizing_lum_unique = sorted(list(set(ionizing_lum_of_simulations)))

        # Dust mass
        dust_mass_unique = sorted(list(set(dust_mass_of_simulations)))

        # Initialize lists to contain the probability for each unique parameter value
        young_lum_probabilities = []
        ionizing_lum_probabilities = []
        dust_mass_probabilities = []

        # Loop over all unique parameter values of the FUV luminosity of the young stars
        for young_lum in young_lum_unique:
            simulation_indices = young_lum_of_simulations == young_lum
            young_lum_probabilities.append(np.sum(probabilities[simulation_indices]))

        # Loop over all unique parameter values of the FUV luminosity of the ionizing stars
        for ionizing_lum in ionizing_lum_unique:
            simulation_indices = ionizing_lum_of_simulations == ionizing_lum
            ionizing_lum_probabilities.append(np.sum(probabilities[simulation_indices]))

        # Loop over all unique parameter values of the dust mass
        for dust_mass in dust_mass_unique:
            simulation_indices = dust_mass_of_simulations == dust_mass
            dust_mass_probabilities.append(np.sum(probabilities[simulation_indices]))

        # Convert the probability lists into NumPy arrays and normalize them
        young_lum_probabilities = np.array(young_lum_probabilities) / sum(young_lum_probabilities)
        ionizing_lum_probabilities = np.array(ionizing_lum_probabilities) / sum(ionizing_lum_probabilities)
        dust_mass_probabilities = np.array(dust_mass_probabilities) / sum(dust_mass_probabilities)

        # Create the probability distributions for the different parameters
        self.distributions["FUV young"] = Distribution.from_probabilities(young_lum_probabilities, young_lum_unique, "FUV young")
        self.distributions["FUV ionizing"] = Distribution.from_probabilities(ionizing_lum_probabilities, ionizing_lum_unique, "FUV ionizing")
        self.distributions["Dust mass"] = Distribution.from_probabilities(dust_mass_probabilities, dust_mass_unique, "Dust mass")

    # -----------------------------------------------------------------

    def animate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating an animation of the SED fitting procedure ...")

        # Initialize a list to contain the frames of the animation
        frames = []

        # Loop over the entries of the chi squared table (sorted by decreasing chi squared)
        for i in range(len(self.chi_squared)):

            # Get the name of the simulation
            simulation_name = self.chi_squared["Simulation name"][i]

            # Determine the path to the corresponding SED plot file
            path = filesystem.join(self.fit_plot_path, simulation_name, "sed.png")

            # Load the image (as a NumPy array)
            image = imageio.imread(path)

            # Add the image to the list of frames
            frames.append(image)

        # Determine the path to the animation file
        path = self.full_output_path("fitting.gif")

        # Create and write the GIF file
        imageio.mimwrite(path, frames)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file of the best simulation
        self.write_best()

        # Write the probability distributions in table format
        self.write_distributions()

        # Plot the probability distributions as histograms
        self.plot_distributions()

    # -----------------------------------------------------------------

    def write_best(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the best model parameters ...")

        # Get the simulation name of the last entry in the chi squared table (the lowest chi squared value)
        simulation_name = self.chi_squared["Simulation name"][len(self.chi_squared)-1]

        # Determine the path to the simulation's ski file
        ski_path = filesystem.join(self.fit_out_path, simulation_name, self.galaxy_name + ".ski")

        # Copy the ski file to the fit/best directory
        filesystem.copy_file(ski_path, self.fit_best_path)

        # Find the corresponding index in the parameter table
        index = tables.find_index(self.parameters, simulation_name, "Simulation name")

        # Get the best parameter values
        self.best_fuv_young = self.parameters["FUV young"][index]
        self.best_fuv_ionizing = self.parameters["FUV ionizing"][index]
        self.best_dust_mass = self.parameters["Dust mass"][index]

        # Write a file with the best parameter values
        path = filesystem.join(self.fit_best_path, "parameters.dat")
        with open(path, 'w') as best_parameters:
            best_parameters.write("FUV young: " + str(self.best_fuv_young) + "\n")
            best_parameters.write("FUV ionizing: " + str(self.best_fuv_ionizing) + "\n")
            best_parameters.write("Dust mass: " + str(self.best_dust_mass) + "\n")

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the probability distributions ...")

        # Loop over the entries in the 'probabilities' table
        for parameter_name in self.distributions:

            # Debugging
            log.debug("Writing the probability distribution of the " + descriptions[parameter_name] + " ...")

            # Determine the path to the resulting table file
            #path = filesystem.join(self.fit_prob_path, parameter_name.lower().replace(" ", "_") + ".dat")

            # Write the table of probabilities for this parameter
            self.distributions[parameter_name].save(path)

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Loop over the different fit parameters
        for parameter_name in self.distributions:

            # Debugging
            log.debug("Plotting the probability distribution of the " + descriptions[parameter_name] + " ...")

            # Get the probability distributinon for
            distribution = self.distributions[parameter_name]
            description = descriptions[parameter_name]

            # Create a plot file for the probability distribution
            path = filesystem.join(self.fit_prob_path, parameter_name + ".pdf")
            distribution.plot(title="Probability of the " + description, path=path, logscale=True)

    # -----------------------------------------------------------------

    def plot_probabilities_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        fig = plt.figure(figsize=(10, 4))

        dust_mass_scale = 1.0
        fuv_young_scale = 1.0
        fuv_ionizing_scale = 1.0

        locplot = [[0.07, 0.15, 0.305, 0.81], [0.375, 0.15, 0.305, 0.81], [0.68, 0.15, 0.305, 0.81]]

        best_dust_mass = self.best_dust_mass / dust_mass_scale
        best_fuv_young = self.best_fuv_young / fuv_young_scale
        best_fuv_ionizing = self.best_fuv_ionizing / fuv_ionizing_scale

        dust_mass_50 = self.percentiles["Dust mass"][1]
        fuv_young_50 = self.percentiles["FUV young"][1]
        fuv_ionizing_50 = self.percentiles["FUV ionizing"][1]

        # Get the values for the different parameters
        dust_mass_values = self.probabilities["Dust mass"]["Dust mass"]
        fuv_young_values = self.probabilities["FUV young"]["FUV young"]
        fuv_ionizing_values = self.probabilities["FUV ionizing"]["FUV ionizing"]

        # Get the associated probabilities
        dust_mass_probs = self.probabilities["Dust mass"]["Probability"]
        fuv_young_probs = self.probabilities["FUV young"]["Probability"]
        fuv_ionizing_probs = self.probabilities["FUV ionizing"]["Probability"]

        fig_a = plt.axes(locplot[0])
        fig_a.set_ylabel('Probability', fontsize=18)
        fig_a.set_xlabel('M$_\mathrm{dust}\, [10^7 M_\odot]$', fontsize=18)
        #abcissa = np.array(params[0]) / MdustScale
        abcissa = dust_mass_values / dust_mass_scale
        width = 0.3e7
        fig_a.bar(abcissa - 0.5 * width, dust_mass_probs, width=width, color='g', ec='k')
        fig_a.plot([best_dust_mass, best_dust_mass], [0, 1], 'k--') # most probable value
        fig_a.plot([dust_mass_50, dust_mass_50], [0, 1], 'r--') # median
        #fig_a.set_xlim(2.5, 6.4)
        #fig_a.set_ylim(0, 0.4)

        fig_b = plt.axes(locplot[1])
        fig_b.set_ylabel('Probability', fontsize=18)
        fig_b.set_xlabel('$\lambda L_\lambda^{\mathrm{young}} \, [10^8 L_\odot]$', fontsize=18)
        #abcissa = np.array(params[1]) / LyoungScale
        abcissa = fuv_young_values / fuv_young_scale
        width = 1.05e16
        fig_b.bar(abcissa - 0.5 * width, fuv_young_probs, width=width, color='g', ec='k')
        fig_b.plot([best_fuv_young, best_fuv_young], [0, 1], 'k--') # most probable value
        fig_b.plot([fuv_young_50, fuv_young_50], [0, 1], 'r--') # median
        #fig_b.set_xlim(6, 19.5)
        #fig_b.set_ylim(0, 0.4)
        fig_b.get_yaxis().set_visible(False)

        fig_c = plt.axes(locplot[2])
        fig_c.set_ylabel('Probability', fontsize=18)
        fig_c.set_xlabel('$\lambda L_\lambda^{\mathrm{ion.}} \, [10^8 L_\odot]$', fontsize=18)
        #abcissa = np.array(params[2]) / LionizingScale
        abcissa = fuv_ionizing_values / fuv_ionizing_scale
        width = 0.10e8
        fig_c.bar(abcissa - 0.5 * width, fuv_ionizing_probs, width=width, color='g', ec='k')
        fig_c.plot([best_fuv_ionizing, best_fuv_ionizing], [0, 1], 'k--')  # most probable value
        fig_c.plot([fuv_ionizing_50, fuv_ionizing_50], [0, 1], 'r--') # median
        #fig_c.set_xlim(0.001, 2.9)
        #fig_c.set_ylim(0, 0.4)
        fig_c.get_yaxis().set_visible(False)

        # Save the figure
        path = filesystem.join(self.fit_prob_path, "probabilities.pdf")
        fig.savefig(path)

        #if False:
        #    "THIS DOES NOT WORK!"
        #    "Need to interpolate somehow because the 3 parameter vectors have different lengths..."
        #    fig = plt.figure(figsize=(10, 4))
        #    plotContours(params1, MdustProb, FyoungProb, FionizedProb)
        #    fig.savefig(outpath + "plotChi2Contours.pdf", format='pdf')

# -----------------------------------------------------------------

# PLOT CONTOURS: http://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

# -----------------------------------------------------------------

def plotContours(params, MdustProb, FyoungProb, FionizedProb):

    MdustScale     = 1.e7 # in 1e7 Msun
    LyoungScale    = 3.846e26 / (1.425e21*0.153) * 1e8 # in 1e8 Lsun
    LionizingScale = 3.53e14 # in Msun/yr

    MdustBestFit    = 69079833.3333     / MdustScale
    FyoungBestFit   = 1.69488344353e+15 / LyoungScale
    FionizedBestFit = 3.53e+14          / LionizingScale

    p1 = MdustProb[None,:]*FionizedProb[:,None]
    p2 = FyoungProb[None,:]*FionizedProb[:,None]
    p3 = MdustProb[None,:]*FyoungProb[:,None]

    locplot = [[0.07,0.15,0.305,0.81],[0.375,0.15,0.305,0.81],[0.68,0.15,0.305,0.81]]

    fig_a = plt.axes(locplot[0])
    fig_a.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_a.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0]) / MdustScale
    y = np.array(params[2]) / LionizingScale
    fig_a.imshow(p1, cmap='gray', interpolation=None,
               origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_a.set_aspect('auto')

    fig_b = plt.axes(locplot[1])
    fig_b.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_b.set_xlabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    x = np.array(params[1]) / LyoungScale
    y = np.array(params[2]) / LionizingScale
    fig_b.imshow(p2, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_b.set_aspect('auto')

    fig_c = plt.axes(locplot[2])
    fig_c.set_ylabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    fig_c.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0])/MdustScale
    y = np.array(params[1])/LyoungScale
    fig_c.imshow(p3, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_c.set_aspect('auto')

# -----------------------------------------------------------------
