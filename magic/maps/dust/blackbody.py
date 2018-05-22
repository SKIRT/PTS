#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.blackbody Contains the BlackBodyDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....magic.core.frame import Frame
from ....core.launch.pts import PTSRemoteLauncher
from .fitter import GridBlackBodyFitter, GeneticBlackBodyFitter
from ....magic.tools import wavelengths
from ....magic.basics.vector import Pixel
from ....core.units.parsing import parse_unit as u
from ....core.basics.configurable import Configurable

# PTS evolution classes and modules
from ....evolve.core.engine import GeneticEngine, RawScoreCriteria
from ....evolve.genomes.list1d import G1DList
from ....evolve.core import mutators
from ....evolve.core import initializators
from ....evolve.core import constants

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

k850 = 0.077  # Maarten (and SKIRT) use different value you so you might want to change to be consistent

# -----------------------------------------------------------------

class BlackBodyDustMapsMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param config:
        :param interactive
        :return:
        """

        # Call the constructor of the base class
        super(BlackBodyDustMapsMaker, self).__init__(*args, **kwargs)

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

        # The results
        self.maps = dict()
        self.origins = dict()
        self.methods = dict()

        # The dust mass error map
        self.error_map = None

        # The process pool
        self.pool = None

        # Create the PTS remote launcher
        self.launcher = PTSRemoteLauncher()

        # The dust masses for each pixel
        self.dust_masses = None

        # The error on the dust mas for each pixel
        self.dust_mass_errors = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the necessary input
        self.create_datacube()

        # 3. Initialize the dust map to the pixel grid of the rebinned images
        self.initialize_map()

        # 4. Create the SEDs
        self.create_seds()

        # 5. Do the fitting
        self.fit()

        # 6. Make the dust map
        self.make_map()

        # 7. Make everything positive
        self.make_positive()

        # 8. Normalize the dust map
        self.normalize()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BlackBodyDustMapsMaker, self).setup()

        # Setup the remote PTS launcher
        if self.config.remote is not None: self.launcher.setup(self.config.remote)

    # -----------------------------------------------------------------

    def create_datacube(self):

        """
        This function ...
        :return:
        """

        # Exclude these bands
        exclude_filters = ["MIPS 70mu", "MIPS 160mu"]

        # Determine the minimum and maximum wavelength
        #min_wavelength = 23. * u("micron") for 'old' version

        # Get the wavelength range
        wavelength_range = wavelengths.black_body_wavelength_range

        # Create the datacube
        self.datacube = self.dataset.create_datacube(wavelength_range.min, wavelength_range.max, exclude=exclude_filters)

        # Create the error cube
        self.errorcube = self.dataset.create_errorcube(wavelength_range.min, wavelength_range.max, exclude=exclude_filters)

        #print(self.datacube.wavelengths())
        #print(self.datacube.frames.keys())
        #print(self.errorcube.wavelengths())
        #print(self.errorcube.frames.keys())

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

        # Initialize the error map
        self.error_map = Frame.zeros(self.datacube.shape)

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

            # Create the pixel ID
            pixel = Pixel(x, y)

            # Get the FIR-submm SED for this pixel
            sed = self.datacube.pixel_sed(x, y, errorcube=self.errorcube)

            # Check values
            if np.any(sed.photometry(asarray=True) < 0): continue

            # Set the pixel
            self.pixels.append(pixel)

            # Add the sed
            self.seds.append(sed)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return len(self.pixels)

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting ...")

        # Make remotely or locally
        if self.config.remote is not None: self.fit_remote()
        else: self.fit_local()

    # -----------------------------------------------------------------

    def fit_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust map locally ...")

        # file1 = open("dustmasses_G15c.dat", "wb")  # write out data dust masses in solar masses (log), temp in K
        # file1.write("#ID \t log_Md \t log_Mderr \t Tc \t Tcerr \t Tw \t Twerr \t log_Md_bestfit \t log_Md_c_bestfit \t log_Md_w_bestfit \t Tc_bestfit \t Tw_bestfit \t Chi2_bestfit \n")

        # data = np.genfromtxt("./G15_dustselect_DR1cat.dat")  # read Herschel fluxes
        # dataID = np.genfromtxt("./G15_dustselect_DR1cat.dat", usecols=0, dtype='string')  # read galaxy ID
        # dist, zhel = np.loadtxt("./G15zD_.txt", unpack=True, usecols=[3, 1])  # dist in Mpc

        # Get the list of wavelengths
        #wavelengths = self.datacube.wavelengths(unit="micron", asarray=True)  # wavelengths in um

        # Get galaxy distance and redshift
        distance = self.galaxy_distance.to("Mpc").value
        redshift = self.galaxy_redshift

        # Create the GridBlackBodyFitter instance
        if self.config.method == "grid":
            fitter = GridBlackBodyFitter()
            #fitter.config.nprocesses = ...

        # Create the GeneticBlackBodyFitter
        elif self.config.method == "genetic": fitter = GeneticBlackBodyFitter()

        # Invalid value
        else: raise ValueError("Invalid value for 'method'")

        # Don't write
        fitter.config.write = False

        # Run the fitter
        fitter.run(seds=self.seds, distances=distance, redshifts=redshift)

        # Get the dust masses
        self.dust_masses = fitter.dust_masses

    # -----------------------------------------------------------------

    def fit_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust map remotely on host '" + self.config.remote + "' ...")

        # INPUT DICTIONARY
        input_dict = dict()

        # CONFIGURATION DICTIONARY
        config_dict = dict()

        # Set input
        input_dict["seds"] = self.seds
        input_dict["distances"] = self.galaxy_distance.to("Mpc").value
        input_dict["redshift"] = self.galaxy_redshift

        # Set configuration
        #config_dict["nprocesses"] = ...

        # Determine PTS command name
        if self.config.method == "grid": command_name = "fit_blackbody_grid"
        elif self.config.method == "genetic": command_name = "fit_blackbody_genetic"
        else: raise ValueError("Invalid value for 'method'")

        # Run PTS remotely and get the output
        dust_masses, dust_mass_errors = self.launcher.run_attached(command_name, config_dict, input_dict, return_output_names=["dust_masses", "dust_mass_errors"], unpack=True)

        # Set the dust masses
        self.dust_masses = dust_masses

        # Set the dust mass errors
        self.dust_mass_errors = dust_mass_errors

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust mass map ...")

        # Loop over the pixels
        for index in range(self.npixels):

            # Get the pixel
            pixel = self.pixels[index]

            # Set the dust mass in the dust mass map
            self.map[pixel] = self.dust_masses[index]

            # Set the dust mass error
            self.error_map[pixel] = self.dust_mass_errors[index]

    # -----------------------------------------------------------------

    def make_map_old(self, plot=False):

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

            # Debugging
            log.debug("Fitting to pixel " + str(index + 1) + " of " + str(npixels) + " ...")

            x = pixels_x[index]
            y = pixels_y[index]

            #x = int(center.x)
            #y = int(center.y)

            # Get the FIR-submm SED for this pixel
            sed = self.datacube.pixel_sed(x, y)

            # Get the fluxes
            fluxes = sed.photometry(unit="Jy", asarray=True)

            if np.any(fluxes < 0): continue

            #print(fluxes)

            #from ....core.plot.sed import SEDPlotter
            #plotter = SEDPlotter()
            # plotter.add_sed(sed, "pixel")
            #plotter.add_sed(sed, "pixel")
            #plotter.add_sed(sed, "pixel")
            #plotter.run()

            #continue

            # The errors
            errors = fluxes * 0.0 + 1.0

            # The distance
            distance = 3.62 * u("Mpc")
            #distance = distance.to("pc").value
            distance = distance.value

            #D = data[i, 1] * 3 * 10 ** 5 / 67.30
            # The distance

            # Do the fit for this pixel
            if self.config.method == "grid": t1, t2, mdust, ratio, dust_mass_range, dust_mass_probs = self.fit_grid(wavelengths, fluxes, errors, distance)
            elif self.config.method == "genetic": t1, t2, mdust, ratio = self.fit_genetic(wavelengths, fluxes, errors, distance)
            else: raise ValueError("Invalid method (" + self.config.method + ")")

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

        # Normalize the error map to relative errors
        self.error_map /= self.map

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
        engine = GeneticEngine(genome)

        # Set options for the engine
        engine.terminationCriteria.set(RawScoreCriteria)
        engine.setMinimax("minimize")
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

        pass

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
