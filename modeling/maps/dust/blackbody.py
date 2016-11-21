#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.blackbody Contains the BlackBodyDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.constants import h,k,c
from multiprocessing import Pool, Process, Manager

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..component import MapsComponent
from ....core.tools.logging import log
from ....magic.core.frame import Frame
from ....magic.misc.spire import SPIRE

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

t1_min = 10.
t1_max = 30.
t1_step = 3.

t2_min = 30.
t2_max = 60.
t2_step = 10.

md_min = 4.
md_max = 6.
#md_step = 0.02
md_step = 0.1

ratio_min = 0.
ratio_max = 1.
ratio_guess = 0.5

# ----------------------------------------------------------------- NEW VERSION:

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

    kv = k850*(850/wa)**2
    flux = kv*10**Md/D**2*blackbody_lam(wa, T1)*2.08*10**11  # 2.08*10**11=1.98*10**30/9.5/10**44*10**2:  1.98*10**30=conversion from Msun to kg,  Mpc**2=9.5*10**44 m**2, Jy=10**-26 W m**-2 Hz-1
    return flux

# -----------------------------------------------------------------

def func2(wa, D, z, Md, T1, Md2, T2):

    """
    Fit two blackbodies to the synthetic data
    """

    wa = wa/(1+z) # shift for redshift (this effectively takes care of the Kcorrection)
    kv1 = k850*(850/wa)**2    #kv=kv0*(vo/v)**beta
    kv2 = k850*(850/wa)**1.5  #changed to beta = 1.5 for consistency with magphys
    flux = kv1*10**Md/D**2/(1+z)*blackbody_lam(wa, T1)*2.08*10**11+kv2*10**Md2/D**2/(1+z)*blackbody_lam(wa, T2)*2.08*10**11   #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def chi_squared(Dust, wa, yorig, yerrorig, D, z, T1, T2):

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

    y=yorig[:]
    yerr=yerrorig[:]

    y[2] = yorig[2] * f250(T1) # KEcorr factor dependant on model SED temperature
    y[3] = yorig[3] * f350(T1)
    y[4] = yorig[4] * f500(T1)
    yerr[2] = yerrorig[2] * f250(T1)
    yerr[3] = yerrorig[3] * f350(T1)
    yerr[4] = yerrorig[4] * f500(T1)
    y2 = func2(wa,D,z, Dust[0],T1,Dust[1],T2)

    som=0
    for i in range(len(wa)):
        som+=((y[i]-y2[i])/yerr[i])**2

    return som

# -----------------------------------------------------------------

def newdata(y,err):

    """
    # generate new SED based on observed SED and its errors
    """

    yb=y[:]
    for i in range(len(y)):
        yb[i]=random.gauss(y[i],err[i])
    return yb

# -----------------------------------------------------------------

class BlackBodyDustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(BlackBodyDustMapMaker, self).__init__()

        # -- Attributes --

        # The datacube
        self.datacube = None

        # The error datacube
        self.errorcube = None

        # The truncation mask
        self.mask = None

        # The pixel indices
        self.pixels = []

        # The pixel seds
        self.seds = []

        # The dust map
        self.map = None

        # The process pool
        self.pool = None

        # The SPIRE instance
        self.spire = SPIRE()

    # -----------------------------------------------------------------

    def run(self, method="grid"):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.create_datacube()

        # 3. Initialize the dust map to the pixel grid of the rebinned images
        self.initialize_map()

        # Create the SEDs
        self.create_seds()

        # 4. Make the dust map
        self.make_map(method)

        # 5. Make everything positive
        self.make_positive()

        # 6. Normalize the dust map
        self.normalize()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BlackBodyDustMapMaker, self).setup()

        # Initialize the process pool
        nprocesses = 8
        self.pool = Pool(processes=nprocesses)

    # -----------------------------------------------------------------

    def create_datacube(self):

        """
        This function ...
        :return:
        """

        # Exclude these bands
        exclude_filters = ["MIPS 70mu", "MIPS 160mu"]

        # Determine the minimum and maximum wavelength
        #min_wavelength = 23. * Unit("micron") for 'old' version
        min_wavelength = 50. * Unit("micron")
        max_wavelength = 1000. * Unit("micron")

        # Create the datacube
        self.datacube = self.dataset.create_datacube(min_wavelength, max_wavelength, exclude=exclude_filters)

        # Create the error cube
        self.errorcube = self.dataset.create_errorcube(min_wavelength, max_wavelength, exclude=exclude_filters)

        # Determine the conversion factor from MJy / sr to Jy/sr
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Determine the conversion factor to Jansky
        pixelscale = self.datacube.average_pixelscale
        pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
        conversion_factor /= pixel_factor

        # Convert the datacube
        self.datacube *= conversion_factor
        self.datacube.unit = "Jy"

        # Convert the errorcube
        self.errorcube *= conversion_factor
        self.errorcube.unit = "Jy"

    # -----------------------------------------------------------------

    def initialize_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the dust map ...")

        # Create a map with the same shape as all the data
        self.map = Frame.zeros(self.datacube.shape)

    # -----------------------------------------------------------------

    def create_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating pixel SEDs ...")

        # Create truncation mask for datacube
        truncation_mask = self.truncation_ellipse.to_pixel(self.datacube.wcs).to_mask(self.datacube.xsize, self.datacube.ysize)

        # Get list of x and y pixels
        pixels_y, pixels_x = np.where(truncation_mask)
        npixels = pixels_x.size

        # Loop over all pixels
        for index in range(npixels):

            # Debugging
            log.debug("Generating SED for pixel " + str(index + 1) + " of " + str(npixels) + " ...")

            # Get x and y index of the pixel
            x = pixels_x[index]
            y = pixels_y[index]

            # Get the FIR-submm SED for this pixel
            sed = self.datacube.pixel_sed(x, y, errorcube=self.errorcube)

            # Set the pixel coordinate
            self.pixels.append((x,y))

            # Add the sed
            self.seds.append(sed)

    # -----------------------------------------------------------------

    def make_map(self, method):

        """
        This function ...
        :param method:
        :return:
        """

        # Inform the user
        log.info("Making the dust map ...")

        # Get the list of wavelengths
        wavelengths = self.datacube.wavelengths(unit="micron", asarray=True) # wavelengths in um

        # Get galaxy distance and redshift
        distance = self.galaxy_distance.to("Mpc").value
        redshift = self.galaxy_redshift

        file1 = open("dustmasses_G15c.dat", "wb")  # write out data dust masses in solar masses (log), temp in K
        file1.write("#ID \t log_Md \t log_Mderr \t Tc \t Tcerr \t Tw \t Twerr \t log_Md_bestfit \t log_Md_c_bestfit \t log_Md_w_bestfit \t Tc_bestfit \t Tw_bestfit \t Chi2_bestfit \n")


        #data = np.genfromtxt("./G15_dustselect_DR1cat.dat")  # read Herschel fluxes
        #dataID = np.genfromtxt("./G15_dustselect_DR1cat.dat", usecols=0, dtype='string')  # read galaxy ID
        #dist, zhel = np.loadtxt("./G15zD_.txt", unpack=True, usecols=[3, 1])  # dist in Mpc

        results = []

        # Loop over the pixels
        for index in range(len(self.seds)):

            # Get the sed
            sed = self.seds[index]

            # Execute
            result = self.pool.apply_async(_fit_one_pixel, args=(index,)) # All simple types (strings)

            # Add the result
            results.append(result)

        #print([res.get(timeout=1) for res in multiple_results])



        # CLOSE AND JOIN THE PROCESS POOL
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------

    def make_map_old(self, method, plot=False):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust map ...")

        # Get the list of wavelengths
        wavelengths = self.datacube.wavelengths(unit="micron", asarray=True)

        # Create truncation mask for datacube
        truncation_mask = self.truncation_ellipse.to_pixel(self.datacube.wcs).to_mask(self.datacube.xsize, self.datacube.ysize)

        # Get list of x and y pixels
        pixels_y, pixels_x = np.where(truncation_mask)
        npixels = pixels_x.size

        #ellipse = self.truncation_ellipse.to_pixel(self.datacube.wcs)
        #center = ellipse.center

        # Loop over all pixels
        for index in range(npixels):
        #for index in range(1):

            # Debugging
            log.debug("Fitting to pixel " + str(index + 1) + " of " + str(npixels) + " ...")

            x = pixels_x[index]
            y = pixels_y[index]

            #x = int(center.x)
            #y = int(center.y)

            # Get the FIR-submm SED for this pixel
            sed = self.datacube.pixel_sed(x, y)

            # Get the fluxes
            fluxes = sed.fluxes(unit="Jy", asarray=True)

            if np.any(fluxes < 0): continue

            #print(fluxes)

            #from ....core.plot.sed import SEDPlotter
            #plotter = SEDPlotter()
            # plotter.add_observed_sed(sed, "pixel")
            #plotter.add_modeled_sed(sed, "pixel")
            #plotter.add_observed_sed(sed, "pixel")
            #plotter.run()

            #continue

            # The errors
            errors = fluxes * 0.0 + 1.0

            # The distance
            distance = 3.62 * Unit("Mpc")
            #distance = distance.to("pc").value
            distance = distance.value

            #D = data[i, 1] * 3 * 10 ** 5 / 67.30
            # The distance

            # Do the fit for this pixel
            if method == "grid": t1, t2, mdust, ratio, dust_mass_range, dust_mass_probs = self.fit_grid(wavelengths, fluxes, errors, distance)
            elif method == "genetic": t1, t2, mdust, ratio = self.fit_genetic(wavelengths, fluxes, errors, distance)
            else: raise ValueError("Invalid method (" + method + ")")

            # Set the dust mass in the dust mass map
            self.map[y, x] = mdust

            plot = True

            if plot:

                plt.plot(dust_mass_range, dust_mass_probs)
                #plt.title(dataID[i])
                plt.xlabel('T_cold (K)')
                plt.show()

                print(mdust, t1, t2)

                plt.errorbar(wavelengths, fluxes, yerr=errors, fmt='bo')
                x2 = np.arange(10., 600., 0.1)
                plt.loglog()
                #plt.title(dataID[i])
                plt.ylim(0.001, 1)
                plt.plot(x2, two_blackbodies(x2, distance, np.log10(ratio) + mdust, t1, np.log10(1.0 - ratio) + mdust, t2), 'r', label="best fit")
                plt.xlabel('Wavelength (microns)')
                plt.ylabel('Flux (Jy)')
                plt.plot(x2, blackbody(x2, distance, np.log10(ratio) + mdust, t1), ':', lw=2, label="cold dust: logMd = %s, Tc= %s K " % (np.log10(ratio) + mdust, t1))
                plt.plot(x2, blackbody(x2, distance, np.log10(1.0 - ratio) + mdust, t2), ':', lw=2, label="warm dust: logMd = %s, Tc= %s K " % (np.log10(1.0 - ratio) + mdust, t2))
                plt.legend(loc=4)
                plt.show()

        #plt.legend(frameon=False)
        #plt.savefig('fit_bb.png')

    # -----------------------------------------------------------------

    def make_positive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Replacing NaNs and negative pixels by zeros ...")

        # Set negatives and NaNs to zero
        self.map.replace_nans(0.0)
        self.map.replace_negatives(0.0)

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the dust map ...")

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

    # -----------------------------------------------------------------

    def fit_grid(self, wa, ydata, yerr, D):

        """
        This function ...
        :return:
        """

        chi2 = float('inf')

        # Parameter ranges
        cold_temp_range = np.arange(t1_min, t1_max, t1_step)
        warm_temp_range = np.arange(t2_min, t2_max, t2_step)
        dust_mass_range = np.arange(md_min, md_max, md_step)

        ncombinations = len(cold_temp_range) * len(warm_temp_range) * len(dust_mass_range)

        dust_mass_probs = np.zeros(len(dust_mass_range))

        # Best values
        cold_temp_best = None
        warm_temp_best = None
        dust_mass_best = None
        ratio_best = None

        ptot = 0
        # Loop over the temperatures and dust masses
        index = 0
        dust_mass_index = 0
        for dust_mass in dust_mass_range:
            for warm_temp in warm_temp_range:
                for cold_temp in cold_temp_range:

                    # Debugging
                    log.debug("Fitting dust ratio for a dust mass of " + str(dust_mass) + ", a warm temperature of "
                              + str(warm_temp) + ", and a cold temperature of " + str(cold_temp) + " (" + str(index+1)
                              + " of " + str(ncombinations) + ") ...")

                    # Optimize
                    popt = minimize(leastsq_old, [ratio_guess], args=(dust_mass, wa, ydata, yerr, D, cold_temp, warm_temp), method='Nelder-Mead', options={'maxiter': 200})

                    Mdr_new = popt.x
                    chi2_new = leastsq_old(Mdr_new, dust_mass, wa, ydata, yerr, D, cold_temp, warm_temp)

                    if chi2 > chi2_new:

                        chi2 = chi2_new

                        # Set best parameters
                        ratio_best = Mdr_new
                        cold_temp_best = cold_temp
                        warm_temp_best = warm_temp
                        dust_mass_best = dust_mass

                        #print(Mddist[ii], Mdr_new, cold_temp_best, warm_temp_best)

                    prob = np.exp(-0.5 * chi2_new)
                    dust_mass_probs[dust_mass_index] += prob
                    ptot += prob

                    index += 1

            dust_mass_index += 1

        # Normalize probabilities
        dust_mass_probs = dust_mass_probs / ptot

        #print("tcold", Mddist[np.where(Mdprob == max(Mdprob))], percentiles(Mddist, Mdprob, 16), percentiles(Mddist, Mdprob, 50), percentiles(Mddist, Mdprob, 84))

        return cold_temp_best, warm_temp_best, dust_mass_best, ratio_best, dust_mass_range, dust_mass_probs

    # -----------------------------------------------------------------

    def fit_genetic(self, wavelengths, ydata, yerr, D):

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

def blackbody_base(lam, T):

    """
    Blackbody as a function of wavelength (um) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """

    from scipy.constants import h, k, c
    lam = 1e-6 * lam  # convert to metres
    return 2. * h * c / (lam ** 3. * (np.exp(h * c / (lam * k * T)) - 1))

# -----------------------------------------------------------------

def blackbody(wa, D, Md, T1):

    """
    This function ...
    :param wa:
    :param D:
    :param Md:
    :param T1:
    :return:
    """

    kv = k850 * (850./wa)**2
    flux = kv * 10**Md/D**2 * blackbody_base(wa, T1)*2.08*10**11 #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def two_blackbodies(wa, D, Md, T1, Md2, T2):

    """
    This function ...
    :param wa:
    :param D:
    :param Md:
    :param T1:
    :param Md2:
    :param T2:
    :return:
    """

    kv = k850 * (850./wa)**2
    flux = kv*10**Md/D**2*blackbody_base(wa, T1)*2.08*10**11+kv*10**Md2/D**2*blackbody_base(wa, T2)*2.08*10**11 #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def leastsq_old(ratio, Md, wa, y, yerr, D, T1, T2):

    """
    This function ...
    :param ratio:
    :param Md:
    :param wa:
    :param y:
    :param yerr:
    :param D:
    :param T1:
    :param T2:
    :return:
    """

    som = 0
    y2 = two_blackbodies(wa, D, np.log10(ratio) + Md, T1, np.log10(1.- ratio) + Md, T2)

    for i in range(len(wa)):
        if i!=0 or y[i]<y2[i]:
            som += ((y[i]-y2[i])/yerr[i])**2

    return som

# -----------------------------------------------------------------

def percentiles(T1dist, T1prob, percentile):

    """
    This function ...
    :param T1dist:
    :param T1prob:
    :param percentile:
    :return:
    """

    percentile = percentile/100.
    perc = 0

    for ii in range(len(T1dist)):

        perc += T1prob[ii]
        if perc > percentile:
            return T1dist[ii]

# -----------------------------------------------------------------

def _fit_one_pixel(wavelengths, ydata, yerr, distance, redshift):

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

                initial_guess = np.array([7., 6.])

                # Do the minimization
                bnds = ((1, None), (0, None))
                popt = minimize(chi_squared, initial_guess, args=(wavelengths, ydatab, yerr, distance, redshift, T1dist[ii], T2), bounds=bnds, method='TNC', options={'maxiter': 20})

                Md_new = popt.x
                chi2_new = chi_squared(Md_new, wavelengths, ydatab, yerr, distance, redshift, T1dist[ii], T2)

                # If chi squared is lower
                if chi2_new < chi2:  # store parameters if chi2<previous chi2

                    chi2 = chi2_new
                    Mds = Md_new
                    T1_sav = T1dist[ii]
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

            # Initial guess
            initial_guess = np.array([7., 6.])

            # Do minimization
            bnds = ((1, None), (0, None))
            popt = minimize(chi_squared, initial_guess, args=(wavelengths, ydatab, yerr, distance, redshift, T1dist[ii], T2), bounds=bnds, method='TNC', options={'maxiter': 20})
            #  minimize(function, [initial guesses], args=(function arguments),bounds=bnds,method='TNC',options={'maxiter':20})  maxiter limited to reduce runtime

            # Get the chi squared value
            Md_new = popt.x
            chi2_new = chi_squared(Md_new, wavelengths, ydatab, yerr, distance, redshift, T1dist[ii], T2)

            # Check if chi squared value is better
            if chi2 > chi2_new:

                # Update the values
                chi2 = chi2_new
                Mds = Md_new
                T1_sav = T1dist[ii]
                T2_sav = T2

        # line = "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n" %
        # (dataID[i], np.mean(Mdslist), np.std(Mdslist), np.mean(Tclist), np.std(Tclist), np.mean(Twlist),
        # , np.log10(10 ** Mds[0] + 10 ** Mds[1]), Mds[0], Mds[1], T1_sav, T2_sav, chi2)

        # file1.write("#ID \t log_Md \t log_Mderr \t Tc \t Tcerr \t Tw \t Twerr \t log_Md_bestfit \t log_Md_c_bestfit \t log_Md_w_bestfit \t Tc_bestfit \t Tw_bestfit \t Chi2_bestfit \n")

        # print(line)
        # file1.write(line)

        # BOOTSTRAPPING

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

        # BEST FIT

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
