#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.fitter Contains the GridBlackBodyFitter and GeneticBlackBodyFitter classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import random
import numpy as np
from scipy.optimize import minimize
from scipy.constants import h,k,c
from multiprocessing import Pool

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....magic.misc.spire import SPIRE
from ....core.basics.configurable import Configurable

# PTS evolution classes and modules
from ....evolve.engine import GAEngine, RawScoreCriteria
from ....evolve.genomes.list1d import G1DList
from ....evolve import mutators
from ....evolve import initializators
from ....evolve import constants

# -----------------------------------------------------------------

# define temperature grid
T1dist = np.arange(10, 30, 1)
T2dist = np.arange(30, 60, 5)

# -----------------------------------------------------------------

k850 = 0.077  # Maarten (and SKIRT) use different value you so you might want to change to be consistent
bootstrapn = 16  # number of bootstrapped SEDS

# De dust mass absorption coefficient is golflengte afhankelijk, dus je kan ook bv k350 gebruiken en herschalen volgens:
# kv = kv0*(vo/v)^beta
# Of gewoon k350 gebruiken en kv = k850*(850/wa)**2 in the func en func2 functies in het script aanpassen naar: kv=k350*(350/wa)**2.

# -----------------------------------------------------------------

def blackbody_lam(lam, T):

    """
    Blackbody as a function of wavelength (um) and temperature (K).
    """

    # Get the wavelength in metres
    lam = 1e-6 * lam

    # in SI units
    return 2.*h*c / (lam**3. * (np.exp(h*c / (lam*k*T)) - 1))

# -----------------------------------------------------------------

def func(wa, D, Md, T1, beta):

    """
    Plot the input model and synthetic data
    """

    kv = k850 * (850/wa)**2
    flux = kv*10**Md/D**2*blackbody_lam(wa, T1)*2.08*10**11  # 2.08*10**11=1.98*10**30/9.5/10**44*10**2:  1.98*10**30=conversion from Msun to kg,  Mpc**2=9.5*10**44 m**2, Jy=10**-26 W m**-2 Hz-1
    return flux

# -----------------------------------------------------------------

def func2(wa, distance, z, Md, T1, Md2, T2):

    """
    Fit two blackbodies to the synthetic data
    """

    wa = wa/(1+z) # shift for redshift (this effectively takes care of the Kcorrection)
    kv1 = k850*(850/wa)**2    #kv=kv0*(vo/v)**beta
    kv2 = k850*(850/wa)**1.5  #changed to beta = 1.5 for consistency with magphys
    flux = kv1*10**Md/distance**2/(1+z)*blackbody_lam(wa, T1)*2.08*10**11+kv2*10**Md2/distance**2/(1+z)*blackbody_lam(wa, T2)*2.08*10**11   #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def chi_squared(Dust, wa, yorig, yerrorig, D, z, T1, T2, corr):

    """
    This function calculates the chi squared value
    :param Dust:
    :param wa:
    :param yorig:
    :param yerrorig:
    :param D:
    :param z:
    :param T1:
    :param T2:
    :return:
    """

    #y = yorig[:]
    #yerr = yerrorig[:]

    #y[2] = yorig[2] * corr[2] # KEcorr factor dependant on model SED temperature
    #y[3] = yorig[3] * corr[3]
    #y[4] = yorig[4] * corr[4]

    #yerr[2] = yerrorig[2] * corr[2]
    #yerr[3] = yerrorig[3] * corr[3]
    #yerr[4] = yerrorig[4] * corr[4]

    # Apply corrections to compare the observed data with the model SED
    y = yorig * corr
    yerr = yerrorig * corr

    y2 = func2(wa, D, z, Dust[0], T1, Dust[1], T2)

    som=0
    for i in range(len(wa)):
        som+=((y[i]-y2[i])/yerr[i])**2

    return som

# -----------------------------------------------------------------

def newdata(y,err):

    """
    Generate new SED based on observed SED and its errors
    """

    yb=y[:]
    for i in range(len(y)):
        yb[i]=random.gauss(y[i],err[i])
    return yb

# -----------------------------------------------------------------

class GridBlackBodyFitter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(GridBlackBodyFitter, self).__init__(config)

        # The SPIRE instance
        self.spire = SPIRE()

        # The pixel or galaxy SEDs
        self.seds = None

        # The distances and redshifts
        self.distances = None
        self.redshifts = None

        # The process pool
        self.pool = None

        # The dust masses
        self.dust_masses = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Do the fitting
        self.fit()

        # 3. Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GridBlackBodyFitter, self).setup(**kwargs)

        # Get the SEDs
        self.seds = kwargs.pop("seds")

        # Get the distances
        self.distances = kwargs.pop("distances")
        if not isinstance(self.distances, list): self.distances = [self.distances] * len(self.seds)

        # Get the redshifts
        self.redshifts = kwargs.pop("redshifts")
        if not isinstance(self.redshifts, list): self.redshifts = [self.redshifts] * len(self.seds)

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        results = []

        # Loop over the pixels
        for index in range(len(self.seds)):

            # Get the sed
            sed = self.seds[index]

            # Get the distance and redshift
            distance = self.distances[index]
            redshift = self.redshifts[index]

            # Get fluxes and errors in Jansky
            ydata = sed.fluxes(unit="Jy", asarray=True)
            yerr = sed.errors(unit="Jy", asarray=True)

            # Execute
            result = self.pool.apply_async(_fit_one_pixel, args=(ydata, yerr, distance, redshift, self.spire,))  # All simple types (strings)

            # Add the result
            results.append(result)

        # Get and set the dust masses
        self.dust_masses = [result.get() for result in results]

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        log.info("Writing ...")

# -----------------------------------------------------------------

class GeneticBlackBodyFitter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(GeneticBlackBodyFitter, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Do the fitting
        self.fit()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GeneticBlackBodyFitter, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        log.info("Fitting ...")

    # -----------------------------------------------------------------

    def fit_(self, wavelengths, ydata, yerr, D):

        """
        This function ...
        :param wavelengths:
        :param ydata:
        :param yerr:
        :param D:
        :return:
        """

        minima = [t1_min, t2_min, md_min, ratio_min]
        maxima = [t1_max, t2_max, md_max, ratio_max]

        # Create the first genome
        genome = G1DList(4)

        # Set genome options
        genome.setParams(minima=minima, maxima=maxima, bestrawscore=0.00, rounddecimal=2)
        genome.initializator.set(initializators.HeterogeneousListInitializerReal)
        # genome.mutator.set(mutators.HeterogeneousListMutatorRealRange)
        genome.mutator.set(mutators.HeterogeneousListMutatorRealGaussian)

        # Create the genetic algorithm engine
        engine = GAEngine(genome)

        # Set options for the engine
        engine.terminationCriteria.set(RawScoreCriteria)
        engine.setMinimax(constants.minimaxType["minimize"])
        engine.setGenerations(5)
        engine.setCrossoverRate(0.5)
        engine.setPopulationSize(100)
        engine.setMutationRate(0.5)

        # Initialize the genetic algorithm
        engine.initialize()

        ###

        engine.evolve()

        # Get best individual parameters
        best = engine.bestIndividual()
        best_t1 = best.genomeList[0]
        best_t2 = best.genomeList[1]
        best_md = best.genomeList[2]
        best_ratio = best.genomeList[3]

        # Return the best fitting parameters
        return best_t1, best_t2, best_md, best_ratio

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        log.info("Writing ...")

# -----------------------------------------------------------------

def _fit_one_pixel(wavelengths, ydata, yerr, distance, redshift, spire):

    """
    This function ...
    :param wavelengths:
    :param ydata:
    :param yerr:
    :return:
    """

    # ydata = [] # observed fluxes in Jy
    # yerr = [] # observed errors in Jy

    # Set initial chi squared value
    chi2 = float("inf")

    # T1prob=np.zeros(len(T1dist))
    ptot = 0

    # Lists for the best parameter values found during bootstrapping
    dust_mass_list = np.zeros(bootstrapn)
    tc_list = np.zeros(bootstrapn)
    tw_list = np.zeros(bootstrapn)

    # Fit bootstrapped SEDs and store in lists of Mds,Tc and Tw
    for iii in range(bootstrapn):

        Mds = T1_sav = T2_sav = None

        # Generate bootstrap SED
        ydatab = newdata(ydata, yerr)

        for ii in range(len(T1dist)):
            for T2 in T2dist:

                T1 = T1dist[ii]

                initial_guess = np.array([7., 6.])

                # KEcorr factor dependant on (cold) model SED temperature
                kecorr_250 = spire.get_ebeam_temperature_psw(T1)
                kecorr_350 = spire.get_ebeam_temperature_pmw(T1)
                kecorr_500 = spire.get_ebeam_temperature_plw(T1)
                corr = [1., 1., kecorr_250, kecorr_350, kecorr_500] # pacs 70, pacs 160, SPIRE 250, SPIRE 350, SPIRE 500

                # Do the minimization
                bnds = ((1, None), (0, None))
                popt = minimize(chi_squared, initial_guess, args=(wavelengths, ydatab, yerr, distance, redshift, T1, T2, corr), bounds=bnds, method='TNC', options={'maxiter': 20})

                Md_new = popt.x
                chi2_new = chi_squared(Md_new, wavelengths, ydatab, yerr, distance, redshift, T1, T2, corr)

                # If chi squared is lower
                if chi2_new < chi2:  # store parameters if chi2<previous chi2

                    chi2 = chi2_new
                    Mds = Md_new
                    T1_sav = T1
                    T2_sav = T2

        # Save best value for parameters for this bootstrap
        dust_mass_list[iii] = np.log10(10 ** Mds[0] + 10 ** Mds[1])
        tc_list[iii] = T1_sav
        tw_list[iii] = T2_sav

    # Best parameter values
    Mds = T1_sav = T2_sav = None

    # Fit observed fluxes (return best fit SED=lowest chisquared)
    for ii in range(len(T1dist)):
        for T2 in T2dist:

            # T1
            T1 = T1dist[ii]

            # Initial guess
            initial_guess = np.array([7., 6.])

            # KEcorr factor dependant on (cold) model SED temperature
            kecorr_250 = spire.f250(T1)
            kecorr_350 = spire.f350(T1)
            kecorr_500 = spire.f500(T1)
            corr = [1., 1., kecorr_250, kecorr_350, kecorr_500]  # pacs 70, pacs 160, SPIRE 250, SPIRE 350, SPIRE 500

            # Do minimization
            bnds = ((1, None), (0, None))
            popt = minimize(chi_squared, initial_guess, args=(wavelengths, ydatab, yerr, distance, redshift, T1, T2, corr), bounds=bnds, method='TNC', options={'maxiter': 20})
            #  minimize(function, [initial guesses], args=(function arguments),bounds=bnds,method='TNC',options={'maxiter':20})  maxiter limited to reduce runtime

            # Get the chi squared value
            Md_new = popt.x
            chi2_new = chi_squared(Md_new, wavelengths, ydatab, yerr, distance, redshift, T1, T2, corr)

            # Check if chi squared value is better
            if chi2 > chi2_new:

                # Update the values
                chi2 = chi2_new
                Mds = Md_new
                T1_sav = T1dist[ii]
                T2_sav = T2

        # OUTPUT

        ## BOOTSTRAPPING

        # Dust mass statistics
        # log_Md
        # log_Mderr
        mean_mds = np.mean(dust_mass_list)
        std_mds = np.std(dust_mass_list)

        # T1 statistics
        # Tc
        # Tcerr
        mean_tc = np.mean(tc_list)
        std_tc = np.std(tc_list)

        # T2 statistics
        # Tw
        # Twerr
        mean_tw = np.mean(tw_list)
        std_tw = np.std(tw_list)

        ## BEST FIT

        # log_Md_bestfit
        logmd = np.log10(10 ** Mds[0] + 10 ** Mds[1])

        # log_Md_c_bestfit
        log_md_c = Mds[0]

        # log_Md_w_bestfit
        log_md_w = Mds[1]

        # Tc_bestfit
        tc = T1_sav

        # Tw_bestfit
        tw = T2_sav

        # Chi2_bestfit
        #chi2 =

        # Return the values
        return mean_mds, std_mds, mean_tc, std_tc, mean_tw, std_tw, logmd, log_md_c, log_md_w, tc, tw, chi2

# -----------------------------------------------------------------
