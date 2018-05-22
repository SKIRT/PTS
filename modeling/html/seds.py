#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.seds Contains the SEDsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import imageio

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class, hover_table_class
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ...core.tools.utils import lazyproperty
from ...core.plot.sed import SEDPlotter
from ...core.tools import filesystem as fs
from ...core.data.sed import SED
from ...core.plot.simulationseds import instrument_name, number_of_columns, contributions
from pts.core.basics.rgbimage import invert_colors

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

class SEDsPageGenerator(HTMLPageComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SEDsPageGenerator, self).__init__(*args, **kwargs)

        # Axis limits
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_flux = None
        self.max_flux = None

        # Plots
        self.instruments_plot = None
        self.contributions_plots = dict() # per instrument
        self.stellar_contributions_plot = None
        self.stellar_contributions_noionizing_plot = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make tables
        self.make_tables()

        # 3. Make plots
        self.make_plots()

        # 4. Generate the html
        self.generate()

        # 5. Write
        self.write()

        # 6. Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDsPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.analysis_runs.load(self.config.analysis_run)

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_output_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.analysis_run.output_path

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # 1. Plot a comparison between the different instruments
        self.plot_instruments()

        # 2. Plot the contributions to the total SEDs
        self.plot_contributions()

        # 3. Plot stellar contributions to the total SED
        self.plot_stellar_contributions()

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_name

    # -----------------------------------------------------------------

    @property
    def reference_sed_name(self):

        """
        Thisfunction ...
        :return:
        """

        return "Observation"

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.analysis_output_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @property
    def image_width(self):

        """
        Thisf unction ...
        :return:
        """

        return 600

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

        # Set ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Loop over the simulated SED files and add the SEDs to the SEDPlotter
        for sed_path in self.sed_paths:

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Load the SED
            sed = SED.from_skirt(sed_path)

            if instr_name == "earth": residuals = True
            else: residuals = False

            # Add the simulated SED to the plotter
            plotter.add_sed(sed, instr_name, residuals=residuals)

        # Add the reference SED
        plotter.add_sed(self.observed_sed, self.reference_sed_name)

        # Determine the path to the plot file
        path = fs.join(self.images_seds_path, "instruments.png")

        # Plot
        plotter.run(output=path)

        # Get the axis limits
        self.min_wavelength = plotter.min_wavelength
        self.max_wavelength = plotter.max_wavelength
        self.min_flux = plotter.min_flux
        self.max_flux = plotter.max_flux

        # Determine path for the dark version
        dark_path = fs.join(self.images_seds_path, "instruments_dark.png")
        create_inverted(path, dark_path)

        # Create image
        self.instruments_plot = html.image(path, width=self.image_width)

    # -----------------------------------------------------------------

    def plot_contributions(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs with different contributions to the total flux for each instrument ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Set ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Check which SED files are produced by a FullInstrument (these files also contain the full SED of the various contributions)
        #for sed_path in self.simulation.seddatpaths():
        for sed_path in self.sed_paths:

            # Determine the name of the corresponding instrument
            instr_name = instrument_name(sed_path, self.prefix)

            # Check how many columns the SED file contains
            ncols = number_of_columns(sed_path)

            # Check the type of the Instrument / SED
            if ncols == 2: continue # SEDInstrument

            # Loop over the different contributions
            for contribution in contributions:

                # Load the SED contribution
                sed = SED.from_skirt(sed_path, contribution=contribution)

                # Add the SED to the plotter
                plotter.add_sed(sed, contribution, residuals=(contribution == "total"))

            # Add the reference SED
            plotter.add_sed(self.observed_sed, self.reference_sed_name)

            # Determine the path to the plot file
            path = fs.join(self.images_seds_path, "contributions_" + instr_name + ".png")

            # Plot
            plotter.run(output=path, min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength, min_flux=self.min_flux, max_flux=self.max_flux)

            # Clear the SED plotter
            plotter.clear()

            # Determine path for the dark version
            dark_path = fs.join(self.images_seds_path, "contributions_" + instr_name + "_dark.png")
            create_inverted(path, dark_path)

            # Create plot
            plot = html.image(path, width=self.image_width)

            # Add the plot
            self.contributions_plots[instr_name] = plot

    # -----------------------------------------------------------------

    def plot_stellar_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the different stellar contributions to the SED ...")

        ## THIS IS TEMPORARY

        playground_run_path = fs.join(self.playground_path, "best_contributions")
        simulations_path = fs.join(playground_run_path, "simulations")

        total_path = fs.join(simulations_path, "total")
        old_path = fs.join(simulations_path, "old")
        young_path = fs.join(simulations_path, "young")
        ionizing_path = fs.join(simulations_path, "ionizing")

        total_out_path = fs.join(total_path, "out")
        old_out_path = fs.join(old_path, "out")
        young_out_path = fs.join(young_path, "out")
        ionizing_out_path = fs.join(ionizing_path, "out")

        total_sed_path = fs.join(total_out_path, "M81_earth_sed.dat")
        old_sed_path = fs.join(old_out_path, "M81_earth_sed.dat")
        young_sed_path = fs.join(young_out_path, "M81_earth_sed.dat")
        ionizing_sed_path = fs.join(ionizing_out_path, "M81_earth_sed.dat")

        # Load the SEDs
        total_sed = SED.from_skirt(total_sed_path)
        old_sed = SED.from_skirt(old_sed_path)
        young_sed = SED.from_skirt(young_sed_path)
        ionizing_sed = SED.from_skirt(ionizing_sed_path)

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Set ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the different stellar contribution SEDs
        plotter.add_sed(total_sed, "total", residuals=True)
        plotter.add_sed(old_sed, "old", residuals=False)
        plotter.add_sed(young_sed, "young", residuals=False)
        plotter.add_sed(ionizing_sed, "ionizing", residuals=False)

        # Add the reference SED
        plotter.add_sed(self.observed_sed, self.reference_sed_name)

        # Determine the path to the plot file
        path = fs.join(self.images_seds_path, "stellar_contributions.png")

        # Plot
        plotter.run(output=path)

        # Determine path for the dark version
        dark_path = fs.join(self.images_seds_path, "stellar_contributions_dark.png")
        create_inverted(path, dark_path)

        # Create image
        self.stellar_contributions_plot = html.image(path, width=self.image_width)

        ## NOW WITHOUT IONIZING

        # Clear
        plotter.clear()

        # Add the different stellar contribution SEDs BUT NOT IONIZING
        plotter.add_sed(total_sed, "total", residuals=True)
        plotter.add_sed(old_sed, "old", residuals=False)
        plotter.add_sed(young_sed, "young", residuals=False)

        # Add the reference SED
        plotter.add_sed(self.observed_sed, self.reference_sed_name)

        # Determine the path to the plot file
        path = fs.join(self.images_seds_path, "stellar_contributions_noionizing.png")

        # Plot
        plotter.run(output=path)

        # Determine path for the dark version
        dark_path = fs.join(self.images_seds_path, "stellar_contributions_noionizing_dark.png")
        create_inverted(path, dark_path)

        # Create image
        self.stellar_contributions_noionizing_plot = html.image(path, width=self.image_width)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the status
        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the status page ...")

        # Create list of css scripts
        #css_paths = css_scripts[:]
        css_paths = []
        css_paths.append(stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Make javascripts urls
        #javascript_paths = javascripts[:]
        javascript_paths = []
        javascript_paths.append(sortable_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        # Theme button
        self.page += html.center(html.make_theme_button(images=True))
        self.page += html.newline + html.newline

        self.page += "INSTRUMENTS:"
        self.page += html.newline + html.newline

        self.page += self.instruments_plot

        self.page += html.newline + html.newline
        self.page += html.line
        self.page += html.newline + html.newline

        self.page += "CONTRIBUTIONS:"
        self.page += html.newline + html.newline

        for instr_name in self.contributions_plots:

            self.page += "  " + instr_name.upper()
            self.page += html.newline + html.newline

            self.page += self.contributions_plots[instr_name]

            self.page += html.newline + html.newline

        self.page += html.line
        self.page += html.newline + html.newline

        self.page += "STELLAR CONTRIBUTIONS:"
        self.page += html.newline + html.newline

        self.page += self.stellar_contributions_noionizing_plot
        self.page += html.newline + html.newline

        self.page += self.stellar_contributions_plot
        self.page += html.newline + html.newline

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write status page
        self.write_page()

    # -----------------------------------------------------------------

    @property
    def page_path(self):

        """
        This function ...
        :return:
        """

        return self.seds_page_path

# -----------------------------------------------------------------

def create_inverted(original_path, new_path):

    """
    This function ...
    :param original_path:
    :param new_path:
    :return:
    """

    # Open the original image
    image = imageio.imread(original_path)

    # Invert the colours
    invert_colors(image)

    # Write the inverted image
    imageio.imwrite(new_path, image)

# -----------------------------------------------------------------
