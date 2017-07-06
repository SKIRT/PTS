#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.generator Contains the HTMLGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import html
from ...core.tools import filesystem as fs
from ..plotting.model import load_test_components, render_components_html

# -----------------------------------------------------------------

page_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
{head}
</head>
<body>
{body}
</body>
</html>
"""

# -----------------------------------------------------------------

underline_template = "<span style='text-decoration: underline;'>{text}</span>"
bold_template = "<span style='font-weight:bold'>{text}</span>"
fontsize_template = "<span style='font-size:{size}px'>{text}</span>"

# -----------------------------------------------------------------

newline = "<br>"

# -----------------------------------------------------------------

class HTMLGenerator(GalaxyModelingComponent):

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
        super(HTMLGenerator, self).__init__(*args, **kwargs)

        # Tables
        self.info_table = None
        self.properties_table = None
        self.status_table = None

        # Models
        self.old_model = None
        self.young_model = None
        self.ionizing_model = None
        self.dust_model = None

        # Pages
        self.status_page = None

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

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(HTMLGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make info table
        self.make_info_table()

        # Make properties table
        self.make_properties_table()

        # Make status table
        self.make_status_table()

    # -----------------------------------------------------------------

    def make_info_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the info table ...")

        # Create the table
        self.info_table = html.SimpleTable(self.galaxy_info.items(), header_row=["Property", "Value"])

    # -----------------------------------------------------------------

    def make_properties_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the properties table ...")

        # Create the table
        self.properties_table = html.SimpleTable(self.galaxy_properties.as_tuples(), header_row=["Property", "Value"])

    # -----------------------------------------------------------------

    def make_status_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the status table ...")

        # Set the background colours of the cells
        bgcolors = [(None, color) for color in self.status.colors]

        # Create the table
        self.status_table = html.SimpleTable(self.status, header_row=["Step", "Status"], bgcolors=bgcolors)

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
        components = load_test_components()

        old_components = {"disk": components["old"], "bulge": components["bulge"]}
        #old_components = {"disk": components["disk"]}
        young_components = {"young": components["young"]}
        ionizing_components = {"ionizing": components["ionizing"]}
        dust_components = {"dust": components["dust"]}

        # Generate HTML
        self.old_model = render_components_html(old_components, only_body=True, width=400, height=500)
        self.young_model = render_components_html(young_components, only_body=True, width=400, height=500)
        self.ionizing_model = render_components_html(ionizing_components, only_body=True, width=400, height=500)
        self.dust_model = render_components_html(dust_components, only_body=True, width=400, height=500)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the status
        self.generate_status()

    # -----------------------------------------------------------------

    def generate_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the status page ...")

        # Create title
        title = fontsize_template.format(size=20, text=underline_template.format(text="Modeling of " + self.galaxy_name))

        # Create titles
        title_info = underline_template.format(text="GALAXY INFO")
        title_properties = underline_template.format(text="GALAXY PROPERTIES")
        title_status = underline_template.format(text="MODELING STATUS")
        title_model = underline_template.format(text="3D MODEL GEOMETRY")

        body = title + newline + newline + title_info + newline + newline + str(self.info_table) + newline + newline
        body += title_properties + newline + newline + str(self.properties_table) + newline + newline
        body += title_status + newline + newline + str(self.status_table) + newline + newline

        body += title_model + newline + newline
        body += bold_template.format(text="Old stars") + newline + newline
        body += self.old_model + newline + newline

        body += bold_template.format(text="Young stars") + newline + newline
        body += self.young_model + newline + newline

        body += bold_template.format(text="Ionizing stars") + newline + newline
        body += self.ionizing_model + newline + newline

        body += bold_template.format(text="Dust") + newline + newline
        body += self.dust_model + newline + newline

        # Create contents
        contents = dict()
        contents["title"] = "Modeling of " + self.galaxy_name
        contents["head"] = ""
        contents["body"] = body

        # Create the status page
        self.status_page = page_template.format(**contents)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write status page
        self.write_status_page()

    # -----------------------------------------------------------------

    @property
    def status_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_status_path

    # -----------------------------------------------------------------

    def write_status_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing status page ...")

        # Write
        fs.write_text(self.status_page_path, self.status_page)

# -----------------------------------------------------------------
