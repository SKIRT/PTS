#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.bestmodellauncher Contains the BestModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.launch.batchlauncher import BatchLauncher
from ...magic.misc.kernels import AnianoKernels
from ...core.basics.filter import Filter

# -----------------------------------------------------------------

class BestModelLauncher(FittingComponent):
    
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
        super(BestModelLauncher, self).__init__(config)

        # -- Attributes --

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski files for simulating the contributions of the various stellar components
        self.ski_contributions = dict()

        # The ski file for generating simulated images for the total model
        self.ski_total = None

        # The paths to the ski files
        self.ski_paths = dict()

        # The Pacs 160 micron filter
        self.pacs_160 = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Adjust the ski template
        self.adjust_ski()

        # 3. Write first, then launch the simulations
        self.write()

        # 4. Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BestModelLauncher, self).setup()

        # Create the Pacs 160 micron filter
        self.pacs_160 = Filter.from_string("Pacs 160")

        # Set options for the batch launcher
        self.set_launcher_options()

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Basic options
        self.launcher.config.shared_input = True  # The input directories (or files) for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job
        self.launcher.config.remotes = [self.config.remote]  # the remote host(s) on which to run the simulations
        #self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        #self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file

        # Simulation analysis options

        ## General
        self.launcher.config.analysis.relative = True

        ## Extraction
        self.launcher.config.analysis.extraction.path = "extr"     # name for the extraction directory
        self.launcher.config.analysis.extraction.progress = True   # extract progress information
        self.launcher.config.analysis.extraction.timeline = True   # extract the simulation timeline
        self.launcher.config.analysis.extraction.memory = True     # extract memory usage information

        ## Plotting
        self.launcher.config.analysis.plotting.path = "plot"
        self.launcher.config.analysis.plotting.progress = True    # plot progress information
        self.launcher.config.analysis.plotting.timeline = True    # plot timeline
        self.launcher.config.analysis.plotting.memory = True      # plot memory usage
        self.launcher.config.analysis.plotting.seds = True        # plot the simulated SEDs
        self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_path  # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.plotting.format = "png"     # plot in PNG format so that an animation can be made from the fit SEDs

        # Set the paths to the kernel for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path(self.pacs_160)
        for filter_name in self.observed_filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        ## Miscellaneous
        self.launcher.config.analysis.misc.path = "misc"          # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.misc.images = True          # create observed images
        self.launcher.config.analysis.misc.fluxes = True          # calculate observed fluxes
        self.launcher.config.analysis.misc.observation_filters = self.observed_filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.misc.make_images_remote = self.config.images_remote
        self.launcher.config.analysis.misc.images_wcs = self.reference_wcs_path
        self.launcher.config.analysis.misc.images_unit = "MJy/sr"
        self.launcher.config.analysis.misc.images_kernels = kernel_paths

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski template ...")

        # 1. Adjust the ski files for simulating the contributions of the various stellar components
        self.adjust_ski_contributions()

        # 2. Adjust the ski file for generating simulated images
        self.adjust_ski_images()
    
    # -----------------------------------------------------------------

    def adjust_ski_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for simulating the contribution of the various stellar components ...")

        # Loop over the different contributions, create seperate ski file instance
        contributions = ["old", "young", "ionizing"]
        component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                           "young": "Young stars",
                           "ionizing": "Ionizing stars"}
        for contribution in contributions:

            # Create a copy of the ski file instance
            ski = self.ski_template.copy()

            # Remove other stellar components
            ski.remove_stellar_components_except(component_names[contribution])

            # Add the ski file to the dictionary
            self.ski_contributions[contribution] = ski

            # Set the ski path
            ski_path = fs.join(self.fit_best_contribution_paths[contribution], self.galaxy_name + ".ski")
            self.ski_paths[contribution] = ski_path

        # Set the path to the ski file for the total model
        self.ski_paths["total"] = fs.join(self.fit_best_total_path, self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    def adjust_ski_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for generating simulated images ...")

        # Create a copy of the ski file instance
        self.ski_total = self.ski_template.copy()

        # Remove all instruments
        self.ski_total.remove_all_instruments()

        # Add the simple instrument
        self.ski_total.add_instrument("earth", self.simple_instrument)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")



    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski files for simulating the contributions of the various stellar components
        self.write_ski_files_contributions()

        # Write the ski file for generating simulated images
        self.write_ski_file_total()

    # -----------------------------------------------------------------

    def write_ski_files_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files for simulating the contribution of the various stellar components ...")

        # Loop over the ski files
        for contribution in self.ski_contributions: self.ski_contributions[contribution].saveto(self.ski_paths[contribution])

    # -----------------------------------------------------------------

    def write_ski_file_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file for creating simulated images for the total model ...")

        # Write the ski file
        self.ski_total.saveto(self.ski_paths["total"])

# -----------------------------------------------------------------
