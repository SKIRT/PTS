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
from scipy import interpolate
from scipy import integrate

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import tables, filesystem
from ...core.tools.logging import log

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

        # Load the parameter table
        self.load_parameters()

        # 3. Calculate the probability distributions
        self.calculate_probabilities()

        # 4. Make an animation of the fitting procedure
        self.animate()

        # Writing
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
        self.chi_squared = tables.from_file(self.chi_squared_table_path)

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
        self.parameters = tables.from_file(self.parameter_table_path)

    # -----------------------------------------------------------------

    def calculate_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probability distributions for each parameter ...")

        return

        # Degrees of freedom = datapoints - free parameters - 1


        params, probabilities = getPDFs(chi2list, parList, dof)

        MdustProb = np.array(probabilities1[0]) + np.array(probabilities2[0])
        FyoungProb = np.array(probabilities1[1]) + np.array(probabilities2[1])
        FionizedProb = np.array(probabilities1[2] + probabilities2[2])

        MdustProb = MdustProb / np.sum(MdustProb)
        FyoungProb = FyoungProb / np.sum(FyoungProb)
        FionizedProb = FionizedProb / np.sum(FionizedProb)

        if 1:
            fig = plt.figure(figsize=(10, 4))
            plotPDFs(fig, params1, MdustProb, FyoungProb, FionizedProb, correction)
            fig.savefig(outpath + "plotPDFs.pdf", format='pdf')

        if 0:
            "THIS DOES NOT WORK!"
            "Need to interpolate somehow because the 3 parameter vectors have different lengths..."
            fig = plt.figure(figsize=(10, 4))
            plotContours(fig, params1, MdustProb, FyoungProb, FionizedProb, correction)
            fig.savefig(outpath + "plotChi2Contours.pdf", format='pdf')

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
        fuv_young = self.parameters["FUV young"][index]
        fuv_ionizing = self.parameters["FUV ionizing"][index]
        dust_mass = self.parameters["Dust mass"][index]

        # Write a file with the best parameter values
        path = filesystem.join(self.fit_best_path, "parameters.dat")
        with open(path, 'w') as best_parameters:
            best_parameters.write("FUV young: " + str(fuv_young) + "\n")
            best_parameters.write("FUV ionizing: " + str(fuv_ionizing) + "\n")
            best_parameters.write("Dust mass: " + str(dust_mass) + "\n")

# -----------------------------------------------------------------

def getPDFs(chi2file, parfile, dof):

    input = np.loadtxt(chi2file)
    chi2list = input[:,1]/dof

    probabilities = np.exp(-0.5*chi2list)

    params = []
    data = open(parfile)
    for line in data:
        if line[0] != '#':
            strings = line.split()
            params.append(map(float, strings[1:4]))

    params = np.array(params)

    parSets = [[]]
    parSetProbabilities = [[]]
    for i in range(0, len(params[0,:])):
        parameter = params[:,i]
        parset = list(set(parameter))
        parset.sort()

        setProbabilities = []
        for j in range(0,len(parset)):
            idx = (parameter == parset[j])
            setProbabilities.append(np.sum(probabilities[idx]))

        parSets.append(parset)
        parSetProbabilities.append(setProbabilities)

    return parSets[1:], parSetProbabilities[1:]

# -----------------------------------------------------------------

def plotPDFs(fig, params, MdustProb, FyoungProb, FionizedProb, correction):

    MdustScale     = 1.e7 # in 1e7 Msun
    LyoungScale    = 3.846e26 / (1.425e21*0.153) * 1e8 # in 1e8 Lsun
    #LionizingScale = 3.53e14 # in Msun/yr
    LionizingScale = 3.846e26 / (1.425e21*0.153) * 1e8 # in 1e8 Lsun


    MdustBestFit    = 69079833.3333     / MdustScale * correction
    FyoungBestFit   = 1.69488344353e+15 / LyoungScale
    FionizedBestFit = 3.53e+14          / LionizingScale

    Mdust_50 = findPercentile(params[0],MdustProb, 50.) / MdustScale * correction
    Mdust_16 = findPercentile(params[0],MdustProb, 16.) / MdustScale * correction
    Mdust_84 = findPercentile(params[0],MdustProb, 84.) / MdustScale * correction

    Fyoung_50 = findPercentile(params[1],FyoungProb, 50.) / LyoungScale
    Fyoung_16 = findPercentile(params[1],FyoungProb, 16.) / LyoungScale
    Fyoung_84 = findPercentile(params[1],FyoungProb, 84.) / LyoungScale

    Fionized_50 = findPercentile(params[2],FionizedProb, 50.) / LionizingScale
    Fionized_16 = findPercentile(params[2],FionizedProb, 16.) / LionizingScale
    Fionized_84 = findPercentile(params[2],FionizedProb, 84.) / LionizingScale

    locplot = [[0.07,0.15,0.305,0.81],[0.375,0.15,0.305,0.81],[0.68,0.15,0.305,0.81]]

    fig_a = plt.axes(locplot[0])
    fig_a.set_ylabel('Probability',fontsize=18)
    fig_a.set_xlabel('M$_\mathrm{dust}\, [10^7 M_\odot]$',fontsize=18)
    abcissa = np.array(params[0])/MdustScale * correction
    width = 0.3
    fig_a.bar(abcissa-0.5*width, MdustProb,width=width, color='g',ec='k')
    fig_a.plot([MdustBestFit,MdustBestFit],[0,1],'k--')
    fig_a.plot([Mdust_50,Mdust_50],[0,1],'r--')
    fig_a.set_xlim(2.5,6.4)
    fig_a.set_ylim(0,0.4)

    fig_b = plt.axes(locplot[1])
    fig_b.set_ylabel('Probability',fontsize=18)
    fig_b.set_xlabel('$\lambda L_\lambda^{\mathrm{young}} \, [10^8 L_\odot]$',fontsize=18)
    abcissa = np.array(params[1])/LyoungScale
    width = 1.05
    fig_b.bar(abcissa-0.5*width, FyoungProb,width=width, color='g',ec='k')
    fig_b.plot([FyoungBestFit,FyoungBestFit],[0,1],'k--')
    fig_b.plot([Fyoung_50,Fyoung_50],[0,1],'r--')
    fig_b.set_xlim(6,19.5)
    fig_b.set_ylim(0,0.4)
    fig_b.get_yaxis().set_visible(False)

    fig_c = plt.axes(locplot[2])
    fig_c.set_ylabel('Probability',fontsize=18)
    fig_c.set_xlabel('$\lambda L_\lambda^{\mathrm{ion.}} \, [10^8 L_\odot]$',fontsize=18)
    abcissa = np.array(params[2])/LionizingScale
    width = 0.10
    fig_c.bar(abcissa-0.5*width, FionizedProb,width=width, color='g',ec='k')
    fig_c.plot([FionizedBestFit,FionizedBestFit],[0,1],'k--')
    fig_c.plot([Fionized_50,Fionized_50],[0,1],'r--')
    fig_c.set_xlim(0.001,2.9)
    fig_c.set_ylim(0,0.4)
    fig_c.get_yaxis().set_visible(False)

    print('Parameter  50th        16th        84th percentile')
    print('Mdust     '+str(Mdust_50)+'  '+str(Mdust_16)+'  '+str(Mdust_84))
    print('Fyoung    '+str(Fyoung_50)+'  '+str(Fyoung_16)+'  '+str(Fyoung_84))
    print('Fionized  '+str(Fionized_50)+'  '+str(Fionized_16)+'  '+str(Fionized_84))
    print('\nThis corresponds to the following input parameters for Mdust, Fyoung and Fionized:\n')
    print(str(Mdust_50*MdustScale)+'    '+str(Fyoung_50*LyoungScale)+'    '+str(Fionized_50*LionizingScale))

# -----------------------------------------------------------------

def plotContours(fig, params, MdustProb, FyoungProb, FionizedProb, correction):

    MdustScale     = 1.e7 # in 1e7 Msun
    LyoungScale    = 3.846e26 / (1.425e21*0.153) * 1e8 # in 1e8 Lsun
    LionizingScale = 3.53e14 # in Msun/yr

    MdustBestFit    = 69079833.3333     / MdustScale * correction
    FyoungBestFit   = 1.69488344353e+15 / LyoungScale
    FionizedBestFit = 3.53e+14          / LionizingScale

    p1 = MdustProb[None,:]*FionizedProb[:,None]
    p2 = FyoungProb[None,:]*FionizedProb[:,None]
    p3 = MdustProb[None,:]*FyoungProb[:,None]

    locplot = [[0.07,0.15,0.305,0.81],[0.375,0.15,0.305,0.81],[0.68,0.15,0.305,0.81]]

    fig_a = plt.axes(locplot[0])
    fig_a.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_a.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0])/MdustScale * correction
    y = np.array(params[2])/LionizingScale
    fig_a.imshow(p1, cmap='gray', interpolation=None,
               origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_a.set_aspect('auto')

    fig_b = plt.axes(locplot[1])
    fig_b.set_ylabel('SFR $[M_\odot \mathrm{yr}^{-1} ]$',fontsize=18)
    fig_b.set_xlabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    x = np.array(params[1])/LyoungScale
    y = np.array(params[2])/LionizingScale
    fig_b.imshow(p2, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_b.set_aspect('auto')

    fig_c = plt.axes(locplot[2])
    fig_c.set_ylabel('F$^{FUV}_\mathrm{young} [10^8 L_\odot]$',fontsize=18)
    fig_c.set_xlabel('M$_\mathrm{dust} [10^7 M_\odot]$',fontsize=18)
    x = np.array(params[0])/MdustScale * correction
    y = np.array(params[1])/LyoungScale
    fig_c.imshow(p3, cmap='gray', interpolation=None,
                 origin='lower', extent=[x[0],x[-1],y[0],y[-1]] )
    fig_c.set_aspect('auto')

# -----------------------------------------------------------------

def findPercentile(parameter, probability, percentile):

    npoints = 10000
    interpfunc = interpolate.interp1d(parameter,probability, kind='linear')

    parRange = np.linspace(min(parameter),max(parameter),npoints)
    interProb = interpfunc(parRange)


    cumInteg = np.zeros(npoints-1)

    for i in range(1,npoints-1):
        cumInteg[i] = cumInteg[i-1] + (0.5*(interProb[i+1] + interProb[i]) * (parRange[i+1] - parRange[i]))

    cumInteg = cumInteg / cumInteg[-1]
    idx = (np.abs(cumInteg-percentile/100.)).argmin()

    return parRange[idx]

# -----------------------------------------------------------------
