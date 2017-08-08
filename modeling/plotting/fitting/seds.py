#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.seds Contains the SEDsPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.data.sed import SED, ObservedSED
from ....core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

class SEDsPLotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SEDsPLotter, self).__init__(*args, **kwargs)

        # The simulated SEDs of all models, as lists for each generation
        self.seds = dict()

        # The SEDs of the different stellar contributions (total, old, young, ionizing)
        self.sed_contributions = OrderedDict()

        # The observed SED
        self.observed_sed = None

    # -----------------------------------------------------------------

    def load(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs ...")

        # Load the observed SEDs
        self.load_observed_sed()

        # Load the simulated SEDs
        self.load_model_seds()

        # Load the contributions to the SEDs
        self.load_sed_contributions()

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed SED ...")

        # Load the observed SED
        self.observed_sed = ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    def load_model_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs of all the fit models ...")

        # Loop over all finished generations
        for generation_name in self.fitting_run.finished_generations:

            # Initialize dictionary to contain SEDs for this generation
            #seds_generation = dict()
            seds_generation = []

            # Loop over all simulations in this generation
            for simulation_name in self.fitting_run.get_simulations_in_generation(generation_name):

                # Determine the path to the 'plot' directory for this simulation
                #plot_path = fs.join(self.fitting_run.generations_path, generation_name, simulation_name, "plot")

                out_path = fs.join(self.fitting_run.generations_path, generation_name, simulation_name, "out")

                # Determine the path to the SED plot
                #sed_path = fs.join(plot_path, self.galaxy_name + "_earth_sed.dat")

                # Determine the path to the SED data file
                sed_path = fs.join(out_path, self.object_name + "_earth_sed.dat")

                # Check whether the SED file is present
                if not fs.is_file(sed_path):

                    # Give warning and continue
                    log.warning("The SED file for simulation " + simulation_name + " of generation " + generation_name + " is missing")
                    continue

                # Load the SED
                sed = SED.from_skirt(sed_path)

                # Add the SED
                seds_generation.append(sed)

            # Add the sed to the list of SEDs
            self.seds[generation_name] = seds_generation

    # -----------------------------------------------------------------

    def load_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs of the various stellar contributions for the best models of each generation ...")

        # Loop over the directories in the fit_best directory
        for path, generation_name in fs.directories_in_path(self.fitting_run.best_path, returns=["path", "name"]):

            # Initialize ...
            seds_generation = dict()

            # Loop over the contributions that have been simulated
            for contribution in fs.directories_in_path(path, returns="name"):

                # Determine the path to the output directory
                out_path = fs.join(path, contribution, "out")

                # Determine the path to the SED file
                sed_path = fs.join(out_path, self.object_name + "_earth_sed.dat")

                # Load the SED
                sed = SED.from_skirt(sed_path)

                # Add the SED
                seds_generation[contribution] = sed

            # Add the SEDS
            self.sed_contributions[generation_name] = seds_generation

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs ...")

        # Plot model SEDs
        self.plot_model_seds()

        # Plot contributions to the SED
        self.plot_sed_contributions()

    # -----------------------------------------------------------------

    def plot_model_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs of all the models ...")

        # Loop over the generations
        for generation_name in self.seds:

            # Create the SEDPlotter object
            plotter = SEDPlotter()

            # Add all model SEDs (these have to be plotted in gray)
            counter = 0
            for sed in self.seds[generation_name]:

                # Add the model SEDs
                plotter.add_sed(sed, str(counter), ghost=True)
                counter += 1

            # Add the 'best' model total SED
            if generation_name in self.sed_contributions: plotter.add_sed(self.sed_contributions[generation_name]["total"], "best") # this SED has to be plotted in black

            # Add the observed SED to the plotter
            plotter.add_sed(self.observed_sed, "observation")

            # Determine the path to the SED plot file
            path = fs.join(self.plot_fitting_seds_path, "model_seds_" + generation_name + ".pdf")

            # Run the plotter
            plotter.run(title=self.object_name, output=path)

    # -----------------------------------------------------------------

    def plot_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs of various contributions ...")

        # Loop over the generations
        for generation_name in self.sed_contributions:

            # Create the SEDPlotter object
            plotter = SEDPlotter()

            # Loop over the contributions
            for contribution in self.sed_contributions[generation_name]:

                # Add the simulated SED to the plotter
                plotter.add_sed(self.sed_contributions[generation_name][contribution], contribution,
                                residuals=(contribution == "total"))

            # Add the observed SED to the plotter
            plotter.add_sed(self.observed_sed, "observation")

            # Determine the path to the SED plot file
            path = fs.join(self.plot_fitting_seds_path, "sed_contributions_" + generation_name + ".pdf")

            # Run the plotter
            plotter.run(output=path, title=self.object_name)

    # -----------------------------------------------------------------

    def create_animation(self):

        """
        This function ...
        :return:
        """

        return

        ## LOAD IMAGES:

        # Inform the user
        log.info("Loading the SED plot files ...")

        # Find all PNG files within the the fit/plot directory
        for path in fs.files_in_path(self.fit_plot_path, extension="png", recursive=True):
            # Load the image (as a NumPy array)
            image = imageio.imread(path)

            # Add the image to the list of frames
            self.frames.append(image)

        ## CREATE ANIMATION:

        # Create the animated GIF instance
        self.animation = Animation(self.frames)

        ## WRITE:

        # Inform the user
        log.info("Writing the GIF animation ...")

        # Determine the path to the animation file
        path = fs.join(self.plot_fitting_path, "fitting.gif")

        # Save the animation as a GIF file
        self.animation.saveto(path)

# -----------------------------------------------------------------
