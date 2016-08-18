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

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.launch.batchlauncher import BatchLauncher
from ...magic.misc.kernels import AnianoKernels
from ...core.basics.filter import Filter
from .wavelengthgrids import create_one_wavelength_grid
from .dustgrids import create_one_dust_grid
from ..core.emissionlines import EmissionLines

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

        # The wavelength grid and dust grid
        self.wavelength_grid = None
        self.dust_grid = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Create the dust grid
        self.create_dust_grid()

        # 4. Adjust the ski template
        self.adjust_ski()

        # 5. Write first, then launch the simulations
        self.write()

        # 6. Launch the simulations
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

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the emission lines instance
        emission_lines = EmissionLines()

        # Fixed wavelengths in the grid
        fixed = [self.i1_filter.pivotwavelength(), self.fuv_filter.pivotwavelength()]

        # Create the grid
        grid, subgrid_npoints, emission_npoints, fixed_npoints = create_one_wavelength_grid(self.config.nwavelengths, emission_lines, fixed)

        # Set the grid
        self.wavelength_grid = grid

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.truncation_ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        x_radius = radius_physical
        y_radius = radius_physical
        z_radius = 3. * Unit("kpc")

        x_min = - x_radius
        x_max = x_radius
        y_min = - y_radius
        y_max = y_radius
        z_min = - z_radius
        z_max = z_radius

        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min

        # Set the scale
        scale = self.config.dg.rel_scale * pixelscale

        # Create the grid
        grid = create_one_dust_grid(self.config.dg.grid_type, scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, self.config.dg.min_level, self.config.dg.max_mass_fraction)

        # Set the grid
        self.dust_grid = grid

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the dust grid
        self.write_dust_grid()

        # Write the ski files for simulating the contributions of the various stellar components
        self.write_ski_files_contributions()

        # Write the ski file for generating simulated images
        self.write_ski_file_total()

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid ...")

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

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

# -----------------------------------------------------------------
