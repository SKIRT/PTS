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
from abc import ABCMeta, abstractmethod
import random
import numpy as np
from scipy.optimize import minimize
from scipy.constants import h,k,c
from multiprocessing import Pool, current_process

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....magic.services.spire import SPIRE
from ....core.basics.configurable import Configurable
from ....core.basics.range import RealRange, QuantityRange
from ....core.units.parsing import parse_unit as u

# PTS evolution classes and modules
from ....evolve.core.engine import GeneticEngine, RawScoreCriteria
from ....evolve.genomes.list1d import G1DList
from ....evolve.core import mutators
from ....evolve.core import initializators
from ....evolve.core import constants

# -----------------------------------------------------------------

t1_range = QuantityRange(10., 30., "K")
t2_range = QuantityRange(30., 60., "K")

# Define temperature grids
T1dist = np.arange(t1_range.min.to("K").value, t1_range.max.to("K").value, 1)
T2dist = np.arange(t2_range.min.to("K").value, t2_range.max.to("K").value, 5)

# -----------------------------------------------------------------

beta = 2.0

k850 = 0.077  # Maarten (and SKIRT) use different value you so you might want to change to be consistent

# -----------------------------------------------------------------

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

def chi_squared_genome(genome, wavelengths, yorig, yerrorig, distance, redshift):

    """
    This function ...
    :param genome:
    :param wavelengths:
    :param yorig:
    :param yerrorig:
    :param distance:
    :param redshift:
    :return:
    """

    T1 = genome.genomeList[0]
    T2 = genome.genomeList[1]
    d1 = genome.genomeList[2]
    d2 = genome.genomeList[3]

    # KEcorr factor dependant on (cold) model SED temperature
    kecorr_250 = spire.get_kcol_temperature_psw(T1 * u("K"), beta, extended=True)
    kecorr_350 = spire.get_kcol_temperature_pmw(T1 * u("K"), beta, extended=True)
    kecorr_500 = spire.get_kcol_temperature_plw(T1 * u("K"), beta, extended=True)
    corr = [1., 1., kecorr_250, kecorr_350, kecorr_500]  # pacs 70, pacs 160, SPIRE 250, SPIRE 350, SPIRE 500

    # chi_squared, initial_guess, args=(wavelengths, ydatab, yerr, distance, redshift, T1, T2, corr)

    Dust = [d1, d2]

    # Calculate and return the chi squared value
    return chi_squared(Dust, wavelengths, yorig, yerrorig, distance, redshift, T1, T2, corr)

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

# The SPIRE instance
spire = SPIRE()

# -----------------------------------------------------------------

class BlackBodyFitter(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BlackBodyFitter, self).__init__(*args, **kwargs)

        # The pixel or galaxy SEDs
        self.seds = None

        # The distances and redshifts
        self.distances = None
        self.redshifts = None

        # The dust masses of the cold dust component
        self.cold_masses = None

        # The dust masses of the warm dust component
        self.warm_masses = None

        # The cold dust temperatures
        self.cold_temperatures = None

        # The warm dust temperatures
        self.warm_temperatures = None

        # The chi squared values
        self.chi_squared = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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
        super(BlackBodyFitter, self).setup(**kwargs)

        # Get the SEDs
        self.seds = kwargs.pop("seds")

        # Get the distances
        self.distances = kwargs.pop("distances")
        if not isinstance(self.distances, list): self.distances = [self.distances] * self.nseds

        # Get the redshifts
        self.redshifts = kwargs.pop("redshifts")
        if not isinstance(self.redshifts, list): self.redshifts = [self.redshifts] * self.nseds

        # Resize the result arrays
        self.cold_masses = np.zeros(self.nseds)
        self.warm_masses = np.zeros(self.nseds)
        self.cold_temperatures = np.zeros(self.nseds)
        self.warm_temperatures = np.zeros(self.nseds)
        self.chi_squared = np.zeros(self.nseds)

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        This function ...
        :return:
        """

        return len(self.seds)

    # -----------------------------------------------------------------

    @property
    def dust_masses(self):

        """
        This function ...
        :return:
        """

        return self.cold_masses + self.warm_masses

    # -----------------------------------------------------------------

    @abstractmethod
    def fit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

class GridBlackBodyFitter(BlackBodyFitter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GridBlackBodyFitter, self).__init__(*args, **kwargs)

        # The process pool
        self.pool = None

        # Error map
        self.cold_mass_errors = None
        self.warm_mass_errors = None
        self.cold_temperature_errors = None
        self.warm_temperature_errors = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GridBlackBodyFitter, self).setup(**kwargs)

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting black body spectra to the pixel SEDs ...")

        results = []

        # Loop over the pixels
        for index in range(len(self.seds)):

            # Debugging
            log.debug("Fitting SED " + str(index+1) + " of " + str(self.nseds) + " ...")

            # Get the sed
            sed = self.seds[index]

            # Get the distance and redshift
            distance = self.distances[index]
            redshift = self.redshifts[index]

            # Get the wavelengths
            wavelengths = sed.wavelengths(unit="micron", asarray=True)

            # Get fluxes and errors in Jansky
            ydata = sed.photometry(unit="Jy", asarray=True)
            yerr_low = sed.errors_min(unit="Jy", asarray=True)
            yerr_high = sed.errors_max(unit="Jy", asarray=True)

            # Execute
            result = self.pool.apply_async(_fit_one_pixel, args=(wavelengths, ydata, yerr_low, yerr_high, distance, redshift,))  # All simple types (strings)

            # Add the result
            results.append(result)

        # Get the results
        for index, result in enumerate(results):

            # Get the output
            tc, tw, md_c, md_w, chi2, bootstrap = result.get()

            # Set properties of this pixel
            self.cold_masses[index] = md_c
            self.warm_masses[index] = md_w
            self.cold_temperatures[index] = tc
            self.warm_temperatures[index] = tw
            self.chi_squared[index] = chi2

            self.cold_mass_errors[index] = bootstrap.cold_mass_error
            self.warm_mass_errors[index] = bootstrap.warm_mass_error
            self.cold_temperature_errors[index] = bootstrap.cold_temperature_error
            self.warm_temperature_errors[index] = bootstrap.warm_temperature_error

        # Close and join the process pool
        #self.pool.close()
        #self.pool.join()

    # -----------------------------------------------------------------

    @property
    def dust_mass_errors(self):

        """
        This function ...
        :return:
        """

        return np.sqrt(self.cold_mass_errors**2 + self.warm_mass_errors**2)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        log.info("Writing ...")

# -----------------------------------------------------------------

class GeneticBlackBodyFitter(BlackBodyFitter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GeneticBlackBodyFitter, self).__init__(*args, **kwargs)

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

        t1_min = t1_range.min.to("K").value
        t1_max = t1_range.max.to("K").value
        t2_min = t1_range.min.to("K").value
        t2_max = t2_range.max.to("K").value

        d1_min = 0.5
        d1_max = 10.0
        d2_min = 0.5
        d2_max = 10.0

        minima = [t1_min, t2_min, d1_min, d2_min]
        maxima = [t1_max, t2_max, d1_max, d2_max]

        # Loop over the pixels
        for index in range(len(self.seds)):

            # Get the sed
            sed = self.seds[index]

            # Get the distance and redshift
            distance = self.distances[index]
            redshift = self.redshifts[index]

            # Get the wavelengths
            wavelengths = sed.wavelengths(unit="micron", asarray=True)

            # Get fluxes and errors in Jansky
            ydata = sed.photometry(unit="Jy", asarray=True)
            yerr_low = sed.errors_min(unit="Jy", asarray=True)
            yerr_high = sed.errors_max(unit="Jy", asarray=True)

            # Observed errors in Jy
            yerr = 0.5 * (-yerr_low + yerr_high)

            # -----------------------------------------------------------------

            # Create the first genome
            genome = G1DList(4)

            # Set genome options
            genome.setParams(minima=minima, maxima=maxima, bestrawscore=0.00, rounddecimal=2)
            genome.initializator.set(initializators.HeterogeneousListInitializerReal)
            genome.mutator.set(mutators.HeterogeneousListMutatorRealRange)
            #genome.mutator.set(mutators.HeterogeneousListMutatorRealGaussian)

            # Set the evaluator function
            chi_squared_lambda = lambda genome: chi_squared_genome(genome, wavelengths, ydata, yerr, distance, redshift)
            genome.evaluator.set(chi_squared_lambda)

            # Create the genetic algorithm engine
            engine = GeneticEngine(genome)

            # Set options for the engine
            engine.terminationCriteria.set(RawScoreCriteria)
            #engine.setMinimax(constants.minimaxType["minimize"])
            engine.setMinimax("minimize")
            engine.setGenerations(5)
            engine.setCrossoverRate(0.5)
            engine.setPopulationSize(100)
            engine.setMutationRate(0.5)

            # Initialize the genetic algorithm
            engine.initialize()

            # Evolve
            engine.evolve()

            # Get best individual parameters
            best = engine.bestIndividual()
            best_t1 = best.genomeList[0]
            best_t2 = best.genomeList[1]
            best_d1 = best.genomeList[2]
            best_d2 = best.genomeList[3]

            chi2 = best.score

            #
            self.cold_masses[index] = best_d1
            self.warm_masses[index] = best_d2
            self.cold_temperatures[index] = best_t1
            self.warm_temperatures[index] = best_t2
            self.chi_squared[index] = chi2

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        log.info("Writing ...")

# -----------------------------------------------------------------

class BootstrapResult(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        self.cold_mass = None
        self.cold_mass_error = None
        self.warm_mass = None
        self.warm_mass_error = None
        self.cold_temperature = None
        self.cold_temperature_error = None
        self.warm_temperature = None
        self.warm_temperature_error = None

# -----------------------------------------------------------------

def _fit_one_pixel(wavelengths, ydata, yerr_low, yerr_high, distance, redshift):

    """
    This function ...
    :param wavelengths:
    :param ydata:
    :param yerr_low:
    :param yerr_high:
    :return:
    """

    print("Process " + str(current_process()) + " working on an SED ...")

    # Observed errors in Jy
    yerr = 0.5 * (-yerr_low + yerr_high)

    # Set initial chi squared value
    chi2 = float("inf")

    # T1prob=np.zeros(len(T1dist))
    ptot = 0

    # Lists for the best parameter values found during bootstrapping
    cold_dust_mass_list = np.zeros(bootstrapn)
    warm_dust_mass_list = np.zeros(bootstrapn)
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
                kecorr_250 = spire.get_kcol_temperature_psw(T1 * u("K"), beta, extended=True)
                kecorr_350 = spire.get_kcol_temperature_pmw(T1 * u("K"), beta, extended=True)
                kecorr_500 = spire.get_kcol_temperature_plw(T1 * u("K"), beta, extended=True)
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
        #dust_mass_list[iii] = np.log10(10 ** Mds[0] + 10 ** Mds[1])
        cold_dust_mass_list[iii] = Mds[0]
        warm_dust_mass_list[iii] = Mds[1]
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
            kecorr_250 = spire.get_kcol_temperature_psw(T1 * u("K"), beta, extended=True)
            kecorr_350 = spire.get_kcol_temperature_pmw(T1 * u("K"), beta, extended=True)
            kecorr_500 = spire.get_kcol_temperature_plw(T1 * u("K"), beta, extended=True)
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
    #mean_mds = np.mean(dust_mass_list)
    #std_mds = np.std(dust_mass_list)

    mean_md_cold = np.mean(cold_dust_mass_list)
    mean_md_warm = np.mean(warm_dust_mass_list)
    std_md_cold = np.std(cold_dust_mass_list)
    std_md_warm = np.std(warm_dust_mass_list)

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
    log_md = np.log10(10 ** Mds[0] + 10 ** Mds[1])

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

    print("")
    print("Result:")
    print("Mean cold Md:", mean_md_cold, "+/-", std_md_cold)
    print("Mean warm Md:", mean_md_warm, "+/-", std_md_warm)
    print("Mean Tc:", mean_tc, "+/-", std_tc)
    print("Mean Tw:", mean_tw, "+/-", std_tw)
    print("log_md:", log_md)
    print("log_md_c:", log_md_c)
    print("log_md_w:", log_md_w)
    print("tc:", tc)
    print("tw:", tw)
    print("chi2:", chi2)
    print("")

    # Return the values
    #return mean_mds, std_mds, mean_tc, std_tc, mean_tw, std_tw, log_md, log_md_c, log_md_w, tc, tw, chi2

    bootstrap = BootstrapResult()
    bootstrap.cold_mass = mean_md_cold
    bootstrap.cold_mass_error = std_md_cold
    bootstrap.cold_temperature = mean_tc
    bootstrap.cold_temperature_error = std_tc
    bootstrap.warm_temperature = mean_tw
    bootstrap.warm_temperature_error = std_tw

    # Return the output
    return tc, tw, 10.**log_md_c, 10.**log_md_w, chi2, bootstrap

# -----------------------------------------------------------------
