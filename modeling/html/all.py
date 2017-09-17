#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.all Contains the AllPagesGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.html.index import IndexPageGenerator
from pts.modeling.html.status import StatusPageGenerator
from pts.modeling.html.data import DataPageGenerator
from pts.modeling.html.preparation import PreparationPageGenerator
from pts.modeling.html.components import ComponentsPageGenerator
from pts.modeling.html.photometry import PhotometryPageGenerator
from pts.modeling.html.maps import MapsPageGenerator
from pts.modeling.html.model import ModelPageGenerator
from pts.modeling.html.fitting import FittingPageGenerator
from pts.modeling.html.datacubes import DatacubesPageGenerator
from pts.modeling.html.fluxes import FluxesPageGenerator
from pts.modeling.html.images import ImagePageGenerator
from pts.modeling.html.attenuation import AttenuationPageGenerator
from pts.modeling.html.colours import ColoursPageGenerator
from pts.modeling.html.heating import HeatingPageGenerator
from pts.core.basics.log import log
from ..component.galaxy import GalaxyModelingComponent
from ..core.progression import create_modeling_progression
from ...core.tools import browser

# -----------------------------------------------------------------

class AllPagesGenerator(GalaxyModelingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(AllPagesGenerator, self).__init__(*args, **kwargs)

        # The modeling progression to use to generate the pages
        self.progression = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Generate the index page
        if self.has_properties: self.generate_index()

        # 3. Generate the status page
        self.generate_status()

        # 4. Generate the data page
        if self.has_images: self.generate_data()

        # 5. Generate the preparation page
        if self.has_prepared: self.generate_preparation()

        # 6. Generate the components page
        if self.has_components: self.generate_components()

        # 7. Generate the photometry page
        if self.has_photometry: self.generate_photometry()

        # 8. Generate the maps page, if maps are chosen to construct a model
        if self.has_model: self.generate_maps()

        # 9. GEnerate the model page
        if self.has_fitting_run: self.generate_model()

        # 10. Generate the fitting page
        if self.has_generation: self.generate_fitting()

        # 11. Generate the datacubes page
        if self.has_datacubes: self.generate_datacubes()

        # 12. Generate the fluxes page
        if self.has_fluxes: self.generate_fluxes()

        # 13. Generate the images page
        if self.has_model_images: self.generate_images()

        # 14. Generate the attenuation page
        if self.has_attenuation: self.generate_attenuation()

        # 15. Generate the colours page
        if self.has_colours: self.generate_colours()

        # 16. Generate the heating page
        if self.has_heating: self.generate_heating()

        # 17. Write
        self.write()

        # 18. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    @property
    def has_properties(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("fetch_properties")

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("fetch_images")

    # -----------------------------------------------------------------

    @property
    def has_prepared(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("prepare_data")

    # -----------------------------------------------------------------

    @property
    def has_components(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("decompose")

    # -----------------------------------------------------------------

    @property
    def has_photometry(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("photometry")

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("build_model")

    # -----------------------------------------------------------------

    @property
    def has_fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("configure_fit")

    # -----------------------------------------------------------------

    @property
    def has_generation(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("fit_sed")

    # -----------------------------------------------------------------

    @property
    def has_datacubes(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_model_images(self):

        """
        This function ....
        :return:
        """

        return self.history.finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.history.finished_any("analyse_attenuation_map", "analyse_attenuation_curve")

    # -----------------------------------------------------------------

    @property
    def has_colours(self):

        """
        This function ...
        :return:
        """

        return self.history.finished("analyse_colours")

    # -----------------------------------------------------------------

    @property
    def has_heating(self):

        """
        This function ...
        :return:
        """

        return self.history.finished_any("analyse_cell_heating", "analyse_projected_heating")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AllPagesGenerator, self).setup(**kwargs)

        # Create the progression
        if "progression" in kwargs: self.progression = kwargs.pop("progression")
        else: self.create_progression()

    # -----------------------------------------------------------------

    def create_progression(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling progression ...")

        # Create
        self.progression = create_modeling_progression(self.config.path)

    # -----------------------------------------------------------------

    def generate_index(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the index page ...")

        # Generate
        # 'generate_index_page'
        generator = IndexPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the status page ...")

        # Generate
        # 'generate_status_page'
        generator = StatusPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the data page ...")

        # Generate
        # 'generate_data_page'
        generator = DataPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the preparation page ...")

        # Generate
        # 'generate_preparation_page'
        generator = PreparationPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the components page ...")

        # Generate
        # 'generate_components_page'
        generator = ComponentsPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the photometry page ...")

        # Generate
        # 'generate_photometry_page'
        generator = PhotometryPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the maps page ...")

        # Generate the maps page
        # 'generate_maps_page'
        generator = MapsPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.model_name = self.progression.model_name
        generator.run()

    # -----------------------------------------------------------------

    def generate_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the model page ...")

        # Generate the model page
        # 'generate_model_page'
        generator = ModelPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.fitting_run = self.progression.fitting_run_name
        generator.run()

    # -----------------------------------------------------------------

    def generate_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the fitting page ...")

        # Generate the fitting page
        # 'generate_fitting_page'
        generator = FittingPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.fitting_run = self.progression.fitting_run_name
        generator.run()

    # -----------------------------------------------------------------

    def generate_datacubes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the datacubes page ...")

        # Generate the datacubes page
        # 'generate_datacubes_page'
        generator = DatacubesPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_fluxes(self):

        """
        This function ...
        :return:
        """

        # Generate the fluxes page
        # 'generate_fluxes_page'
        generator = FluxesPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_images(self):

        """
        This function ...
        :return:
        """

        # Generate the images page
        # 'generate_images_page'
        generator = ImagePageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.run()

    # -----------------------------------------------------------------

    def generate_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the attenuation page ...")

        # Generate the attenuation page
        # 'generate_attenuation_page'
        generator = AttenuationPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.analysis_run = self.progression.analysis_run_name
        generator.run()

    # -----------------------------------------------------------------

    def generate_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the colours page ...")

        # Generate the colours page
        # 'generate_colours_page'
        generator = ColoursPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.analysis_run = self.progression.analysis_run_name
        generator.run()

    # -----------------------------------------------------------------

    def generate_heating(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the heating page ...")

        # Generate the heating page
        # 'generate_heating_page'
        generator = HeatingPageGenerator()
        generator.config.path = self.config.path
        generator.config.replot = self.config.replot
        generator.config.analysis_run = self.progression.analysis_run_name
        generator.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the pages ...")

        # Open
        browser.open_path(self.environment.html_status_path)

# -----------------------------------------------------------------
