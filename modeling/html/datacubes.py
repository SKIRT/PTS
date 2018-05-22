#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.datacubes Contains the DatacubesPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, sortable_table_class
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...core.tools.html import HTMLPage, updated_footing, make_page_width
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.misc.images import is_total_datacube, instrument_name
from ...magic.core.datacube import DataCube
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

page_width = 600
thumbnail_title = "Thumbnail"

# -----------------------------------------------------------------

class DatacubesPageGenerator(HTMLPageComponent):

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
        super(DatacubesPageGenerator, self).__init__(*args, **kwargs)

        # The plot paths per instrument
        self.plot_paths = dict()

        # The slider for each instrument
        self.sliders = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make tables
        self.make_tables()

        # Make plots
        self.make_plots()

        # Make the sliders
        self.make_sliders()

        # 2. Generate the html
        self.generate()

        # 3. Write
        self.write()

        # 4. Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DatacubesPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

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

    @lazyproperty
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.wavelengths()

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_strings(self):

        """
        This function ...
        :return:
        """

        return [tostr(wavelength) for wavelength in self.wavelengths]

    # -----------------------------------------------------------------

    @property
    def first_wavelength_string(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_strings[0]

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_name

    # -----------------------------------------------------------------

    @lazyproperty
    def total_datacube_paths(self):

        """
        Thisn function ...
        :return:
        """

        paths = []
        for path in fs.files_in_path(self.analysis_output_path, extension="fits", startswith=self.prefix):
            if not is_total_datacube(path): continue
            paths.append(path)
        return paths

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Loop over the total datacube paths
        for path in self.total_datacube_paths:

            # Get the instrument name
            instr_name = instrument_name(path, self.prefix)

            # Create a plot directory for this instrument
            plot_path = fs.create_directory_in(self.images_datacubes_path, instr_name)

            # Load the datacube
            datacube = DataCube.from_file(path, self.wavelength_grid)

            # Plot paths for this instrument
            plot_paths = []

            # Loop over the frames, create PNG image for each frame
            for index in range(datacube.nframes):

                # Get the frame
                frame = datacube.frames[index]

                # Determine frame plot path
                frame_plot_path = fs.join(plot_path, str(index) + ".png")

                # Add the path
                plot_paths.append(frame_plot_path)

                # Check if already present
                if fs.is_file(frame_plot_path): continue

                # Save as PNG
                vmin, vmax = frame.saveto_png(frame_plot_path, colours=self.config.colours,
                                                      interval=self.config.interval,
                                                      scale=self.config.scale, alpha=self.config.alpha_method,
                                                      peak_alpha=self.config.peak_alpha)

            # Set plot paths for this instrument
            self.plot_paths[instr_name] = plot_paths

    # -----------------------------------------------------------------

    @property
    def image_width(self):

        """
        This function ...
        :return:
        """

        return 600

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        Thisf unction ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def make_sliders(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the image sliders ...")

        # Loop over the instruments
        for instr_name in self.plot_paths:

            # Create the slider
            slider = html.make_image_slider(instr_name, self.plot_paths[instr_name], self.wavelength_strings, self.first_wavelength_string,
                                            width=self.image_width, height=self.image_height, basic=True,
                                            img_class="pixelated", extra_img_class="pixelated")

            # Set the slider
            self.sliders[instr_name] = slider

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the page
        self.generate_page()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

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
        #javascript_paths.append(preview_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        #classes = dict()
        #classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        self.page += html.newline

        self.page += html.line
        self.page += html.newline

        # Add the sliders
        for instr_name in self.sliders:

            # Add the name
            self.page += instr_name.upper() + " INSTRUMENT"
            self.page += html.newline

            # Add the image slider
            self.page += self.sliders[instr_name]

            # Add line
            self.page += html.line
            self.page += html.newline

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

        return self.datacubes_page_path

# -----------------------------------------------------------------
