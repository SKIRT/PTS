#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.playground Contains the MappingsPlayground and BruzualCharlot playground classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import tempfile
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.simulation.execute import SkirtExec
from ...core.data.sed import SED
from ...core.prep.templates import get_oneparticle_template
from ...core.basics.log import log
from ...core.units.unit import parse_unit as u

# -----------------------------------------------------------------

k = 1.3806488e-23    # Boltmann constant in J/K
Zsun = 0.0122        # solar metallicity according to Groves et al. 2008

# -----------------------------------------------------------------

colors = ('c','b','r','m')

# -----------------------------------------------------------------

def plot_mappings_examples(path, pressure=True, compactness=True, fpdr=True, metallicity=True, format="pdf",
                           temp_path=None):

    """
    This function ...
    :param path:
    :param pressure:
    :param covering_factor:
    :param fpdr:
    :param metallicity:
    :param format:
    :param temp_path:
    :return:
    """

    # Create the playground
    playground = MappingsPlayground(temp_path=temp_path)

    # Define plot paths
    pressure_path = fs.join(path, "hiiregion-logP." + format)
    compactness_path = fs.join(path, "hiiregion-logC." + format)
    fpdr_path = fs.join(path, "hiiregion-fpdr.pdf")
    metallicity_path = fs.join(path, "hiiregion-z-a.pdf")
    metallicity2_path = fs.join(path, "hiiregion-z-b.pdf")

    # Produce plot of varying pressure
    if pressure: playground.vary_pressure(4., 8., 1., Zsun, 5., 1., path=pressure_path)

    # Produce plot of varying compactness
    if compactness: playground.vary_compactness(4., 6.5, 1., Zsun, 5., 1., path=compactness_path)

    # Produce plot of varying fPDR (cloud covering factor)
    if fpdr: playground.vary_fpdr(0., 1., 1., Zsun, 5., 5., path=fpdr_path)

    # Produce plot of varying metallicity
    if metallicity: playground.vary_metallicity(0.05 * Zsun, 2. * Zsun, 1., 5., 5., 0., path=metallicity_path)

    # Another with varying metallicity
    if metallicity: playground.vary_metallicity(0.05 * Zsun, 2. * Zsun, 1., 5., 5., 1., path=metallicity2_path)

# -----------------------------------------------------------------

class MappingsPlayground(object):
    
    """
    This class...
    """

    def __init__(self, temp_path=None):

        """
        The constructor ...
        :param temp_path:
        :return:
        """

        # Call the constructor of the base class
        super(MappingsPlayground, self).__init__()

        # -- Attributes --

        # Determine the path to a temporary directory
        self.temp_path = tempfile.gettempdir() if temp_path is None else temp_path

        # Load the template ski file
        self.ski = get_oneparticle_template()

        # Set the number of wavelengths
        self.ski.set_nwavelengths(1000)

    # -----------------------------------------------------------------

    def vary_pressure(self, min_logp, max_logp, sfr, met, logc, fpdr, nvalues=3, path=None):

        """
        This function ...
        :param min_logp:
        :param max_logp:
        :param sfr:
        :param met:
        :param logc:
        :param fpdr:
        :param nvalues:
        :param path:
        :return:
        """

        # Calculate the range of logp values
        logp_range = np.linspace(min_logp, max_logp, nvalues)

        seds = dict()

        # Loop over the range of logp values
        for logp in logp_range:

            # Simulate the SED and add it to the dictionary
            sed = self.simulate_sed(logp, sfr, met, logc, fpdr)
            seds[(logp, sfr, met, logc, fpdr)] = sed

        # Plot the seds
        self.plot_seds(seds, path)

    # -----------------------------------------------------------------

    def vary_compactness(self, min_logc, max_logc, sfr, met, logp, fpdr, nvalues=3, path=None):

        """
        This function ...
        :param min_logc:
        :param max_logc:
        :param sfr:
        :param met:
        :param logp:
        :param fpdr:
        :param nvalues:
        :param path:
        :return:
        """

        # Calculate the range of logc values
        logc_range = np.linspace(min_logc, max_logc, nvalues)

        seds = dict()

        for logc in logc_range:

            # Simulate the SED and add it to the dictionary
            sed = self.simulate_sed(logp, sfr, met, logc, fpdr)
            seds[(logp, sfr, met, logc, fpdr)] = sed

        # Plot the seds
        self.plot_seds(seds, path)

    # -----------------------------------------------------------------

    def vary_fpdr(self, min_fpdr, max_fpdr, sfr, met, logc, logp, nvalues=3, path=None):

        """
        This function ...
        :param min_fpdr:
        :param max_fpdr:
        :param sfr:
        :param met:
        :param logc:
        :param logp:
        :param nvalues:
        :param path:
        :return:
        """

        # Calculate the range of fpdr values
        fpdr_range = np.linspace(min_fpdr, max_fpdr, nvalues)

        seds = dict()

        for fpdr in fpdr_range:

            # Simulate the SED and add it to the dictionary
            sed = self.simulate_sed(logp, sfr, met, logc, fpdr)
            seds[(logp, sfr, met, logc, fpdr)] = sed

        # Plot the seds
        self.plot_seds(seds, path)

    # -----------------------------------------------------------------

    def vary_metallicity(self, min_met, max_met, sfr, logc, logp, fpdr, nvalues=3, path=None):

        """
        This function ...
        :param min_met:
        :param max_met:
        :param sfr:
        :param logc:
        :param logp:
        :param fpdr:
        :param nvalues:
        :param path:
        :return:
        """

        # Calculate the range of metallicity values
        metallicity_range = np.linspace(min_met, max_met, nvalues)

        seds = dict()

        for met in metallicity_range:

            # Simulate the SED and add it to the dictionary
            sed = self.simulate_sed(logp, sfr, met, logc, fpdr)
            seds[(logp, sfr, met, logc, fpdr)] = sed

        # Plot the seds
        self.plot_seds(seds, path)

    # -----------------------------------------------------------------

    def simulate_sed(self, logp, sfr, met, logc, fpdr, silent=True, output_path=None):

        """
        This function ...
        :param logp:
        :param sfr:
        :param met:
        :param logc:
        :param fpdr:
        :param silent:
        :param output_path:
        :return:
        """

        # Parameter string
        parameter_string = "SFR{} Z{} logC{} logP{} fPDR{}".format(sfr, met / Zsun, logc, logp, fpdr)

        # Create output path
        if output_path is None:
            # Create a directory within the temporary directory
            output_path = fs.join(self.temp_path, parameter_string)
            fs.create_directory(output_path)

        # Log debug
        log.info("Simulating SED in temporary directory '" + output_path + "' ...")

        # Convert pressure to Pascal
        pressure = 10. ** logp * k * 1e6

        # Determine the path to the temporary data file
        particle_path = fs.join(output_path, "oneparticle.dat")

        # Output an appropriate particle data file
        datafile = open(particle_path, 'w')
        datafile.write("0 0 0 1 {} {} {} {} {}\n".format(sfr, met, logc, pressure, fpdr))
        datafile.close()

        # Determine the path to the ski file
        ski_path = fs.join(output_path, "oneparticle.ski")

        # Save the ski file to the temporary directory
        self.ski.saveto(ski_path)

        # Perform the SKIRT simulation
        simulation = SkirtExec().execute(ski_path, brief=True, inpath=output_path, outpath=output_path, silent=silent)[0]

        # Load the fluxes, convert them to luminosities in erg/s
        sedpath = simulation.seddatpaths()[0]

        # Get distance
        distance = simulation.instrumentdistance(unit='m') * u("m")

        # Simpler
        sed = SED.from_skirt(sedpath, distance=distance)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def plot_seds(self, seds, path=None):

        """
        This function ...
        :param seds:
        :param path:
        :return:
        """

        # setup the figure
        figure = plt.figure(figsize=(10, 6))
        plt.xscale('log')
        plt.yscale('log')

        counter = 0

        for (logp, sfr, met, logc, fpdr), sed in seds.items():

            wavelengths = sed.wavelengths(unit="micron", add_unit=False)
            luminosities = sed.photometry(unit="erg/s", density=True, add_unit=False)

            # plot the SED
            plt.plot(wavelengths, luminosities, color=colors[counter], label="SFR={} Z={} logC={} logP={} fPDR={}".format(sfr, met/Zsun, logc, logp, fpdr))

            counter += 1

        # add axis labels, legend and title
        plt.grid('on')
        plt.xlabel(r"$\lambda\,(\mu \mathrm{m})$", fontsize='medium')
        plt.ylabel(r"$\nu F_\nu\,(\mathrm{erg}\,\mathrm{s}^{-1})$", fontsize='medium')
        plt.xlim(0.1, 1000)
        # ymax = lambdaLlambdav.max()
        # plt.ylim(ymax*1.1e-3, ymax*1.1)
        plt.ylim(10 ** 40., 10 ** 43.5)
        plt.legend(loc='best', prop={'size': 'small'})

        # Save the figure
        if path is not None:
            plt.savefig(path, bbox_inches='tight', pad_inches=0.25)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------

class BruzualCharlotPlayground(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        pass

# -----------------------------------------------------------------
