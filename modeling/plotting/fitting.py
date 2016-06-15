#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting Contains the FittingPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.distribution import Distribution
from ...core.launch.timing import TimingTable
from ..core.transmission import TransmissionCurve
from ...core.plot.transmission import TransmissionPlotter
from ...core.plot.wavelengthgrid import WavelengthGridPlotter
from ...core.plot.grids import plotgrids
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.logfile import LogFile
from ...core.simulation.skifile import SkiFile
from ..core.emissionlines import EmissionLines
from ..core.sed import load_example_mappings_sed, load_example_bruzualcharlot_sed, load_example_zubko_sed
from ..core.sed import SED, ObservedSED
from ...magic.plot.imagegrid import ResidualImageGridPlotter
from ...magic.core.frame import Frame
from ...core.plot.sed import SEDPlotter
from ...magic.basics.skyregion import SkyRegion
from ..basics.models import load_model
from ..misc.geometryplotter import GeometryPlotter

# -----------------------------------------------------------------

class FittingPlotter(PlottingComponent):
    
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
        super(FittingPlotter, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski = None

        # The wavelength grids
        self.lowres_wavelength_grid = None
        self.highres_wavelength_grid = None

        # The transmission curves
        self.transmission_curves = dict()

        # The distribution of dust cells per level
        self.cell_distribution_lowres = None
        self.cell_distribution_highres = None

        # The runtimes
        self.runtimes = None

        # The SEDs of all the models
        self.seds = []

        # The SEDs of the different stellar contributions (total, old, young, ionizing)
        self.sed_contributions = OrderedDict()

        # The observed SED
        self.observed_sed = None

        # The simulated images
        self.simulated_images = dict()

        # The observed imags
        self.observed_images = dict()

        # The geometries
        self.geometries = dict()

    # -----------------------------------------------------------------

    def run(self, features=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file
        self.load_ski_file()

        # 3. Load the wavelength grids
        self.load_wavelength_grids()

        # 4. Load the transmission curves
        self.load_transmission_curves()

        # 5. Load the dust cell tree data
        self.load_dust_cell_tree()

        # 6. Load the runtimes
        self.load_runtimes()

        # 7. Load the observed SED
        self.load_observed_sed()

        # 8. Load the SEDs of the various models
        self.load_model_seds()

        # 9. Load the SEDs for the various contribution
        self.load_sed_contributions()

        # 10. Load the simulated images
        self.load_images()

        # 11. Load the geometries
        self.load_geometries()

        # 12. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_ski_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Determine the path to the initialized ski file
        path = fs.join(self.fit_path, self.galaxy_name + ".ski")

        # Load the ski file
        if fs.is_file(path): self.ski = SkiFile(path)

    # -----------------------------------------------------------------

    def load_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grids used for the fitting ...")

        # Load the low-resolution wavelength grid
        self.load_low_res_wavelength_grid()

        # Load the high-resolution wavelength grid
        self.load_high_res_wavelength_grid()

    # -----------------------------------------------------------------

    def load_low_res_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the low-resolution wavelength grid ...")

        # Determine the path to the wavelength grid file
        path = fs.join(self.fit_path, "in", "wavelengths_lowres.txt")

        # Load the wavelength grid
        if fs.is_file(path): self.lowres_wavelength_grid = WavelengthGrid.from_skirt_input(path)

    # -----------------------------------------------------------------

    def load_high_res_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the high-resolution wavelength grid ...")

        # Determine the path to the wavelength grid file
        path = fs.join(self.fit_path, "in", "wavelengths_highres.txt")

        # Load the wavelength grid
        if fs.is_file(path): self.highres_wavelength_grid = WavelengthGrid.from_skirt_input(path)

    # -----------------------------------------------------------------

    def load_transmission_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the transmission curves for the observation filters ...")

        # Load the observed SED
        sed = ObservedSED.from_file(self.observed_sed_path)
        for fltr in sed.filters():

            # Create the transmission curve
            transmission = TransmissionCurve.from_filter(fltr)

            # Normalize the transmission curve
            transmission.normalize(value=1.0, method="max")

            # Add the transmission curve to the dictionary
            self.transmission_curves[str(fltr)] = transmission

    # -----------------------------------------------------------------

    def load_dust_cell_tree(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust cell tree information ...")

        # Load the low-resolution dust cell tree information
        self.load_low_res_dust_cell_tree()

        # Load the high-resolution dust cell tree information
        self.load_high_res_dust_cell_tree()

    # -----------------------------------------------------------------

    def load_low_res_dust_cell_tree(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the log file of the dust grid generating simulation
        fit_grid_lowres_path = fs.join(self.fit_path, "grid", "low-res")
        log_file_path = fs.join(fit_grid_lowres_path, self.galaxy_name + "_log.txt")

        # If the log file does not exist, break
        if not fs.is_file(log_file_path): return

        # Open the log file
        log_file = LogFile(log_file_path)

        # Get the distribution of cells per level of the tree
        self.cell_distribution_lowres = log_file.tree_leaf_distribution

    # -----------------------------------------------------------------

    def load_high_res_dust_cell_tree(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the log file of the dust grid generating simulation
        fit_grid_highres_path = fs.join(self.fit_path, "grid", "high-res")
        log_file_path = fs.join(fit_grid_highres_path, self.galaxy_name + "_log.txt")

        # If the log file does not exist, break
        if not fs.is_file(log_file_path): return

        # Open the log file
        log_file = LogFile(log_file_path)

        # Get the distribution of cells per level of the tree
        self.cell_distribution_highres = log_file.tree_leaf_distribution

    # -----------------------------------------------------------------

    def load_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the runtimes ...")

        # Determine the path to the timing table
        timing_table_path = fs.join(self.fit_path, "timing.dat")

        try:
            # Load the timing table
            timing_table = TimingTable.read(timing_table_path)
        except ValueError: return # ValueError: Column Timestamp failed to convert

        # Keep a list of all the runtimes recorded for a certain remote host
        self.runtimes = defaultdict(lambda: defaultdict(list))

        # Loop over the entries in the timing table
        for i in range(len(timing_table)):

            # Get the ID of the host and the cluster name for this particular simulation
            host_id = timing_table["Host id"][i]
            cluster_name = timing_table["Cluster name"][i]

            # Get the parallelization properties for this particular simulation
            cores = timing_table["Cores"][i]
            threads_per_core = timing_table["Threads per core"][i]
            processes = timing_table["Processes"][i]

            # Get the number of photon packages (per wavelength) used for this simulation
            packages = timing_table["Packages"][i]

            # Get the total runtime
            runtime = timing_table["Total runtime"][i]

            parallelization = (cores, threads_per_core, processes)

            # Add the runtime to the list of runtimes
            self.runtimes[host_id][packages, parallelization].append(runtime)

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

        # Determine the path to the fit/out directory
        fit_out_path = fs.join(self.fit_path, "out")

        # Loop over all directories in the fit/out directory
        for path in fs.directories_in_path(fit_out_path):

            # Check whether the SED file is present in the simulation's output directory
            sed_path = fs.join(path, "out", self.galaxy_name + "_earth_sed.dat")
            if not fs.is_file(sed_path): continue

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the sed to the list of SEDs
            self.seds.append(sed)

    # -----------------------------------------------------------------

    def load_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs of the various stellar contributions ...")

        # Determine the path to the fit/best directory
        fit_best_path = fs.join(self.fit_path, "best")

        # Determine the path to the fit/best/images directory
        images_path = fs.join(fit_best_path, "images")

        # Determine the path to the SED file from the 'images' simulation
        sed_path = fs.join(images_path, self.galaxy_name + "_earth_sed.dat")

        # Load the SED if the file exists
        if fs.is_file(sed_path):

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the total SED to the dictionary of SEDs
            self.sed_contributions["total"] = sed

        # Add the SEDs of the simulations with the individual stellar populations
        contributions = ["old", "young", "ionizing"]
        for contribution in contributions:

            # Determine the output path for this simulation
            out_path = fs.join(fit_best_path, contribution)

            # Determine the path to the SED file
            sed_path = fs.join(out_path, self.galaxy_name + "_earth_sed.dat")

            # Skip if the file does not exist
            if not fs.is_file(sed_path): continue

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the SED to the dictionary of SEDs
            self.sed_contributions[contribution] = sed

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Load simulated images
        self.load_simulated_images()

        # Load observed images
        self.load_observed_images()

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Determine the path to the fit/best directory
        fit_best_path = fs.join(self.fit_path, "best")

        # Determine the path to the fit/best/images directory
        out_path = fs.join(fit_best_path, "images")

        # Loop over all FITS files found in the fit/best/images directory
        for path, name in fs.files_in_path(out_path, extension="fits", returns=["path", "name"], contains="__"):

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # Get the filter name
            #filter_name = name.split("__")[1]

            # Open the image
            frame = Frame.from_file(path)

            # Get the filter name
            filter_name = str(frame.filter)

            # Add the image frame to the dictionary
            self.simulated_images[filter_name] = frame

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over all FITS files found in the 'truncated' directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Check whether a simulated image exists for this band
            if name not in self.simulated_images:
                log.warning("The simulated version of the " + name + " image could not be found, skipping " + name + " data ...")
                continue

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # The filter name is the image name
            filter_name = name

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.observed_images[filter_name] = frame

    # -----------------------------------------------------------------

    def load_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the geometries ...")

        # Determine the path to the fit/geometries directory
        geometries_path = fs.join(self.fit_path, "geometries")

        # Load all geometries in the directory
        for path, name in fs.files_in_path(geometries_path, extension="mod", returns=["path", "name"]):

            # Load the geometry
            model = load_model(path)

            # Add the geometry
            self.geometries[name] = model

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the wavelength grid used for the fitting
        self.plot_wavelengths()

        # Plot the dust grid
        self.plot_dust_grids()

        # Plot the distribution of dust cells for the different tree levels
        self.plot_dust_cell_distribution()

        # Plot the distributions of the runtimes on different remote systems
        if self.runtimes is not None: self.plot_runtimes()

        # Plot the SEDs of the models
        if len(self.seds) > 0: self.plot_model_seds()

        # Plot the SEDs
        if len(self.sed_contributions) > 0: self.plot_sed_contributions()

        # Plot the images
        if len(self.simulated_images) > 0: self.plot_images()

        # Plot the geometries
        if len(self.geometries) > 0: self.plot_geometries()

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grid ...")

        # Plot the low-resolution and high-resolution wavelength grids with the observation filters
        self.plot_wavelengths_filters(self.lowres_wavelength_grid, "wavelengths_filters_lowres.pdf")
        self.plot_wavelengths_filters(self.highres_wavelength_grid, "wavelengths_filters_highres.pdf")

        # Plot the wavelengths with the SEDs of stars and dust
        self.plot_wavelengths_seds(self.lowres_wavelength_grid, "wavelengths_seds_lowres.pdf")
        self.plot_wavelengths_seds(self.highres_wavelength_grid, "wavelengths_seds_highres.pdf")

    # -----------------------------------------------------------------

    def plot_wavelengths_filters(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with")

        # Create the transmission plotter
        plotter = TransmissionPlotter()

        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the transmission curves
        for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label], label)

        # Add the wavelength points
        for wavelength in grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength, min_transmission=0.0, max_transmission=1.05)

    # -----------------------------------------------------------------

    def plot_wavelengths_seds(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with SEDs of stars and dust ...")

        # Create the wavelength grid plotter
        plotter = WavelengthGridPlotter()

        # Set title
        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the wavelength grid
        plotter.add_wavelength_grid(grid, "fitting simulations")

        # Add MAPPINGS SFR SED
        mappings_sed = load_example_mappings_sed()
        plotter.add_sed(mappings_sed, "MAPPINGS")

        # Add Bruzual-Charlot stellar SED
        bc_sed = load_example_bruzualcharlot_sed()
        plotter.add_sed(bc_sed, "Bruzual-Charlot")

        # Add Zubko dust emission SED
        zubko_sed = load_example_zubko_sed()
        plotter.add_sed(zubko_sed, "Zubko")

        # Add emission lines
        emission_lines = EmissionLines()
        for line in emission_lines: plotter.add_emission_line(line)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength)

    # -----------------------------------------------------------------

    def plot_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust grid ...")

        # Plot the low-resolution dust grid
        self.plot_low_res_dust_grid()

        # Plot the high-resolution dust grid
        self.plot_high_res_dust_grid()

    # -----------------------------------------------------------------

    def plot_low_res_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Create a SkirtSimulation instance for the grid generating simulation
        fit_grid_lowres_path = fs.join(self.fit_path, "grid", "low-res")
        ski_path = fs.join(fit_grid_lowres_path, self.galaxy_name + ".ski")
        simulation = SkirtSimulation(ski_path=ski_path, outpath=fit_grid_lowres_path)

        # Plot the grid
        plotgrids(simulation, output_path=self.plot_fitting_path, silent=True, prefix="lowres")

    # -----------------------------------------------------------------

    def plot_high_res_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Create a SkirtSimulation instance for the grid generating simulation
        fit_grid_highres_path = fs.join(self.fit_path, "grid", "high-res")
        ski_path = fs.join(fit_grid_highres_path, self.galaxy_name + ".ski")
        simulation = SkirtSimulation(ski_path=ski_path, outpath=fit_grid_highres_path)

        # Plot the grid
        plotgrids(simulation, output_path=self.plot_fitting_path, silent=True, prefix="highres")

    # -----------------------------------------------------------------

    def plot_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust cell distribution ...")

        # Plot low-resolution dust cell distribution
        if self.cell_distribution_lowres is not None: self.plot_low_res_dust_cell_distribution()

        # Plot high-resolution dust cell distribution
        if self.cell_distribution_highres is not None: self.plot_high_res_dust_cell_distribution()

    # -----------------------------------------------------------------

    def plot_low_res_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Set plot title and path
        title = "Dust cells in each tree level"
        path = fs.join(self.plot_fitting_path, "cells_tree_lowres.pdf")

        # Make the plot
        self.cell_distribution_lowres.plot(title=title, path=path)

    # -----------------------------------------------------------------

    def plot_high_res_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Set plot title and path
        title = "Dust cells in each tree level"
        path = fs.join(self.plot_fitting_path, "cells_tree_highres.pdf")

        # Make the plot
        self.cell_distribution_highres.plot(title=title, path=path)

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        bins = 15

        # Loop over the different remote hosts
        for host_id in self.runtimes:

            # Loop over the different configurations (packages, parallelization)
            for packages, parallelization in self.runtimes[host_id]:

                runtimes = self.runtimes[host_id][packages, parallelization]

                distribution = Distribution.from_values(runtimes, bins)

                path = fs.join(self.plot_fitting_path, "runtimes_" + host_id + "_" + str(parallelization) + "_" + str(packages) + ".pdf")
                title = host_id + " " + str(parallelization) + " " + str(packages)

                # Plot the distribution
                distribution.plot(title=title, path=path)

    # -----------------------------------------------------------------

    def plot_model_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs of all the models ...")

        # Create the SEDPlotter object
        plotter = SEDPlotter(self.galaxy_name)

        # Add all model SEDs (these have to be plotted in gray)
        counter = 0
        for sed in self.seds:
            plotter.add_modeled_sed(sed, str(counter), ghost=True)
            counter += 1

        # Add the 'best' model total SED
        plotter.add_modeled_sed(self.sed_contributions["total"], "best") # this SED has to be plotted in black

        # Add the observed SED to the plotter
        plotter.add_observed_sed(self.observed_sed, "observation")

        # Determine the path to the SED plot file
        path = fs.join(self.plot_fitting_path, "model_seds.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs of various contributions ...")

        # Create the SEDPlotter object
        plotter = SEDPlotter(self.galaxy_name)

        # Loop over the simulated SEDs of the various stellar contributions
        for label in self.sed_contributions:

            # Add the simulated SED to the plotter
            plotter.add_modeled_sed(self.sed_contributions[label], label, residuals=(label=="total"))

        # Add the observed SED to the plotter
        plotter.add_observed_sed(self.observed_sed, "observation")

        # Determine the path to the SED plot file
        path = fs.join(self.plot_fitting_path, "sed_contributions.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual images ...")

        # Get the ellipse
        path = fs.join(self.truncation_path, "ellipse.reg")
        region = SkyRegion.from_file(path)
        ellipse = region[0]

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Image residuals")

        # Create list of filter names sorted by increasing wavelength
        sorted_filter_names = sorted(self.observed_images.keys(), key=lambda key: self.observed_images[key].filter.pivotwavelength())

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in sorted_filter_names:

            observed = self.observed_images[filter_name]
            simulated = self.simulated_images[filter_name]

            plotter.add_row(observed, simulated, filter_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(ellipse.bounding_box)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "images.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the geometries of the model ...")

        # Create the geometry plotter
        plotter = GeometryPlotter()

        # Add the geometries
        for label in self.geometries: plotter.add_geometry(self.geometries[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "geometries.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
