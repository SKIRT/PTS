#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.simulationseds Contains the SimulationSEDPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from .sed import SEDPlotter
from ..data.sed import SED, ObservedSED
from ..tools import filesystem as fs
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

contributions = ["total", "direct", "scattered", "dust", "dustscattered", "transparent"]

# -----------------------------------------------------------------

class SimulationSEDPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SimulationSEDPlotter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The SKIRT simulation object
        self.simulation = None

        # Plot limits for the different instruments
        self.min_wavelengths = dict()
        self.max_wavelengths = dict()
        self.min_fluxes = dict()
        self.max_fluxes = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param simulation:
        :return:
        """

        # 2. Make SED plots for each instrument seperately
        self.plot_single()

        # 2. Plot SED comparing different instruments
        if self.config.instruments: self.plot_instruments()

        # 3. Plot various contributions to the SEDs
        if self.config.contributions: self.plot_contributions()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        For ALL instruments
        :return:
        """

        min_wavelength = None
        for instr_name in self.instrument_names:
            if min_wavelength is None: min_wavelength = self.min_wavelengths[instr_name]
            elif self.min_wavelengths[instr_name] < min_wavelength: min_wavelength = self.min_wavelengths[instr_name]
        return min_wavelength

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        For ALL instruments
        :return:
        """

        max_wavelength = None
        for instr_name in self.instrument_names:
            if max_wavelength is None: max_wavelength = self.max_wavelengths[instr_name]
            elif self.max_wavelengths[instr_name] > max_wavelength: max_wavelength = self.max_wavelengths[instr_name]
        return max_wavelength

    # -----------------------------------------------------------------

    @property
    def min_flux(self):

        """
        For ALL instruments
        :return:
        """

        min_flux = None
        for instr_name in self.instrument_names:
            if min_flux is None: min_flux = self.min_fluxes[instr_name]
            elif self.min_fluxes[instr_name] < min_flux: min_flux = self.min_fluxes[instr_name]
        return min_flux

    # -----------------------------------------------------------------

    @property
    def max_flux(self):

        """
        For ALL instruments
        :return:
        """

        max_flux = None
        for instr_name in self.instrument_names:
            if max_flux is None: max_flux = self.max_fluxes[instr_name]
            elif self.max_fluxes[instr_name] > max_flux: max_flux = self.max_fluxes[instr_name]
        return max_flux

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the simulation SED plotter ...")

        # Set the attributes to default values
        self.simulation = None

        # Clear the ranges
        self.min_wavelengths = dict()
        self.max_wavelengths = dict()
        self.min_fluxes = dict()
        self.max_fluxes = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationSEDPlotter, self).setup()

        # Make a local reference to the simulation object
        self.simulation = kwargs.pop("simulation")

    # -----------------------------------------------------------------

    @lazyproperty
    def instrument_names(self):

        """
        This function ...
        :return:
        """

        names = []

        # Loop over the simulated SED files and add the
        for sed_path in self.simulation.seddatpaths():

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Add the name
            names.append(instr_name)

        return names

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        return self.simulation.prefix()

    # -----------------------------------------------------------------

    def plot_single(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs for the different instruments separately ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Set plotting options
        plotter.config.library = self.config.library
        plotter.config.ignore_filters = self.config.ignore_filters

        # Loop over the simulated SED files and add the
        for sed_path in self.simulation.seddatpaths():

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Load the SED
            sed = SED.from_skirt(sed_path) # contribution = "total"

            # Add the simulated SED to the plotter
            plotter.add_sed(sed, instr_name)

            # Check if reference SED is defined, and add it
            if self.config.reference_seds is not None:

                # Add the reference SEDs
                for reference_sed_path in self.config.reference_seds:

                    # Determine name
                    reference_sed_name = fs.strip_extension(fs.name(reference_sed_path))

                    # Add the reference SED
                    reference_sed = ObservedSED.from_file(reference_sed_path)
                    plotter.add_sed(reference_sed, reference_sed_name)

            # Determine the path to the plot file
            path = self.output_path_file("sed_" + instr_name + "." + self.config.format)

            # Make the plot
            plotter.run(title=self.simulation.name, output=path)

            # Get the axis limits
            min_wavelength = plotter.min_wavelength
            max_wavelength = plotter.max_wavelength
            min_flux = plotter.min_flux
            max_flux = plotter.max_flux

            # Reset the plotter
            plotter.clear()

            # Set the axis limits
            self.min_wavelengths[instr_name] = min_wavelength
            self.max_wavelengths[instr_name] = max_wavelength
            self.min_fluxes[instr_name] = min_flux
            self.max_fluxes[instr_name] = max_flux

    # -----------------------------------------------------------------

    def plot_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting comparison of the SEDs of the different instruments ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Set plotting options
        plotter.config.library = self.config.library
        plotter.config.ignore_filters = self.config.ignore_filters

        # Loop over the simulated SED files and add the SEDs to the SEDPlotter
        for sed_path in self.simulation.seddatpaths():

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Load the SED
            sed = SED.from_skirt(sed_path) # contribution = "total"

            # Add the simulated SED to the plotter
            plotter.add_sed(sed, instr_name)

        # Check if reference SED is defined, and add it
        if self.config.reference_seds is not None:

            # Add the reference SEDs
            for reference_sed_path in self.config.reference_seds:

                # Determine name
                reference_sed_name = fs.strip_extension(fs.name(reference_sed_path))

                # Add the reference SED
                reference_sed = ObservedSED.from_file(reference_sed_path)
                plotter.add_sed(reference_sed, reference_sed_name)

        # Determine the path to the plot file
        path = self.output_path_file("sed." + self.config.format)
        plotter.run(title=self.simulation.name, output=path, min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength, min_flux=self.min_flux, max_flux=self.max_flux)

    # -----------------------------------------------------------------

    def plot_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs with different contributions to the total flux for each instrument ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Set plotting options
        plotter.config.library = self.config.library
        plotter.config.ignore_filters = self.config.ignore_filters

        # Check which SED files are produced by a FullInstrument (these files also contain the full SED of the various contributions)
        for sed_path in self.simulation.seddatpaths():

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Check how many columns the SED file contains
            ncols = number_of_columns(sed_path)

            # Check the type of the Instrument / SED
            if ncols == 2: continue  # SEDInstrument

            # Loop over the different contributions
            for contribution in contributions:

                # Load the SED contribution
                sed = SED.from_skirt(sed_path, contribution=contribution)

                # Add the SED to the plotter
                plotter.add_sed(sed, contribution, residuals=(contribution == "total"))

            # Check whether reference SEDs must be added
            if self.config.reference_seds is not None:

                # Add the reference SEDs
                for reference_sed_path in self.config.reference_seds:

                    # Determine name
                    reference_sed_name = fs.strip_extension(fs.name(reference_sed_path))

                    # Add the reference SED
                    reference_sed = ObservedSED.from_file(reference_sed_path)
                    plotter.add_sed(reference_sed, reference_sed_name)

            # Determine the path to the plot file
            path = self.output_path_file("sed_" + instr_name + "_contributions." + self.config.format)

            # Get the axis limits
            min_wavelength = self.min_wavelengths[instr_name]
            max_wavelength = self.max_wavelengths[instr_name]
            min_flux = self.min_fluxes[instr_name]
            max_flux = self.max_fluxes[instr_name]

            # Plot
            plotter.run(output=path, min_wavelength=min_wavelength, max_wavelength=max_wavelength, min_flux=min_flux, max_flux=max_flux)

            # Clear the SED plotter
            # Doesn't clear the config (e.g. ignore filters)
            plotter.clear()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------

def instrument_name(sed_path, prefix):

    """
    This function ...
    :param sed_path:
    :param prefix:
    :return:
    """

    return fs.name(sed_path).split("_sed.dat")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------

def number_of_columns(sed_path):

    """
    This function ...
    :param sed_path:
    :return:
    """

    with open(sed_path, 'r') as f:

        ncols = 0
        for line in f:

            if "# column" not in line: break
            else: ncols = int(line.split("column ")[1].split(": ")[0])

    return ncols

# -----------------------------------------------------------------
