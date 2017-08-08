#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.model Contains the ModelPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, stylesheet_url, table_class
from ...core.tools import html
from ..fitting.run import load_fitting_run
from ..tests.data import M81TestData

# -----------------------------------------------------------------

widget_width = 400
widget_height = 500
widget_style = "minimal"

# -----------------------------------------------------------------

data = M81TestData()

# -----------------------------------------------------------------

class ModelPageGenerator(HTMLPageComponent):

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
        super(ModelPageGenerator, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

        # Models
        self.old_model = None
        self.young_model = None
        self.ionizing_model = None
        self.dust_model = None

        # Tables
        self.old_bulge_table = None
        self.old_disk_table = None
        self.young_table = None
        self.ionizing_table = None
        self.dust_table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

        # Make tables
        self.make_tables()

        # Make plots
        self.make_plots()

        # Generaet the html
        self.generate()

        # Write
        self.write()

        # Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelPageGenerator, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = load_fitting_run(self.config.path, self.config.fitting_run)

    # -----------------------------------------------------------------

    @property
    def definition(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_definition

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make
        self.make_geometry_tables()

    # -----------------------------------------------------------------

    def make_geometry_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making geometry tables ...")

        # Old bulge
        self.make_old_bulge_geometry_table()

        # Old disk
        self.make_old_disk_geometry_table()

        # Young stars
        self.make_young_geometry_table()

        # Ionizing stars
        self.make_ionizing_geometry_table()

        # Dust
        self.make_dust_geometry_table()

    # -----------------------------------------------------------------

    def make_old_bulge_geometry_table(self):

        """
        This function ...
        :return:
        """

        self.old_bulge_table = html.SimpleTable(data.bulge.as_tuples(), header_row=["Property", "Value"], css_class=table_class)

    # -----------------------------------------------------------------

    def make_old_disk_geometry_table(self):

        """
        This function ...
        :return: 
        """

        self.old_disk_table = html.SimpleTable(data.old_deprojection.as_tuples(), header_row=["Property", "Value"], css_class=table_class)

    # -----------------------------------------------------------------

    def make_young_geometry_table(self):

        """
        This function ...
        :return:
        """

        self.young_table = html.SimpleTable(data.young_deprojection.as_tuples(), header_row=["Property", "Value"], css_class=table_class)

    # -----------------------------------------------------------------

    def make_ionizing_geometry_table(self):

        """
        This function ...
        :return:
        """

        self.ionizing_table = html.SimpleTable(data.ionizing_deprojection.as_tuples(), header_row=["Property", "Value"], css_class=table_class)

    # -----------------------------------------------------------------

    def make_dust_geometry_table(self):

        """
        This function ...
        :return:
        """

        self.dust_table = html.SimpleTable(data.dust_deprojection.as_tuples(), header_row=["Property", "Value"], css_class=table_class)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Model
        self.make_model_plots()

    # -----------------------------------------------------------------

    def make_model_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making model plots ...")

        # TEMPORARY: LOAD TEST COMPONENTS
        components = data.components

        # Group the components
        old_components = {"disk": components["old"], "bulge": components["bulge"]}
        young_components = {"young": components["young"]}
        ionizing_components = {"ionizing": components["ionizing"]}
        dust_components = {"dust": components["dust"]}

        # New
        #old = self.definition.old_stars_deprojection
        #young = self.definition.young_stars_deprojection
        #ionizing = self.definition.ionizing_stars_deprojection
        #dust = self.definition.dust_deprojection

        # Generate HTML
        self.old_model = data.render_components_html(old_components, only_body=True, width=widget_width, height=widget_height, style=widget_style)
        self.young_model = data.render_components_html(young_components, only_body=True, width=widget_width, height=widget_height, style=widget_style)
        self.ionizing_model = data.render_components_html(ionizing_components, only_body=True, width=widget_width, height=widget_height, style=widget_style)
        self.dust_model = data.render_components_html(dust_components, only_body=True, width=widget_width, height=widget_height, style=widget_style)

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
        log.info("Generating the status page ...")

        body = self.heading

        # Create titles
        title_model = html.underline_template.format(text="3D MODEL GEOMETRY")

        body += title_model + html.newline + html.newline
        body += html.bold_template.format(text="Old stars") + html.newline + html.newline
        body += self.old_model + html.newline + html.newline

        body += html.bold_template.format(text="Young stars") + html.newline + html.newline
        body += self.young_model + html.newline + html.newline

        body += html.bold_template.format(text="Ionizing stars") + html.newline + html.newline
        body += self.ionizing_model + html.newline + html.newline

        body += html.bold_template.format(text="Dust") + html.newline + html.newline
        body += self.dust_model + html.newline + html.newline

        # Create contents
        contents = dict()
        contents["title"] = self.title
        contents["head"] = html.link_stylesheet_header_template.format(url=stylesheet_url)
        contents["body"] = body
        contents["style"] = self.style

        # Create the page
        self.page = html.page_template.format(**contents)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write models page
        self.write_page()

    # -----------------------------------------------------------------

    @property
    def page_path(self):

        """
        This function ...
        :return:
        """

        return self.model_page_path

# -----------------------------------------------------------------
