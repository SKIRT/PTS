#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.component Contains the HTMLPageComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod, abstractproperty

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import html
from ...core.tools import filesystem as fs
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import time
from ...core.tools import browser

# -----------------------------------------------------------------

# Stylesheets
stylesheet_url = "https://users.ugent.be/~sjversto/stylesheet.css"
slider_stylesheet_url = "https://users.ugent.be/~sjversto/slider.css"

# Scripts
sortable_url = "https://users.ugent.be/~sjversto/sorttable.js"
preview_url = "https://users.ugent.be/~sjversto/preview.js"
slider_url = "https://users.ugent.be/~sjversto/slider.js"

# -----------------------------------------------------------------

table_class = "realtable"
hover_table_class = "hovertable"
sortable_table_class = "sortable"

# -----------------------------------------------------------------

page_style = "ugentstyle"

# -----------------------------------------------------------------

top_title_size = 24
title_size = 20

# -----------------------------------------------------------------

data_page_filename = "data.html"
photometry_page_filename = "photometry.html"
components_page_filename = "components.html"
preparation_page_filename = "preparation.html"
model_page_filename = "model.html"
maps_page_filename = "maps.html"
fitting_page_filename = "fitting.html"
datacubes_page_filename = "datacubes.html"
fluxes_page_filename = "fluxes.html"
images_page_filename = "images.html"
heating_page_filename = "heating.html"
colours_page_filename = "colours.html"
attenuation_page_filename = "attenuation.html"

# -----------------------------------------------------------------

class HTMLPageComponent(GalaxyModelingComponent):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HTMLPageComponent, self).__init__(*args, **kwargs)

        # Page
        self.page = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(HTMLPageComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @abstractmethod
    def make_tables(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def make_plots(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def images_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_images_path

    # -----------------------------------------------------------------

    @property
    def index_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_index_path

    # -----------------------------------------------------------------

    @property
    def status_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_status_path

    # -----------------------------------------------------------------

    @property
    def has_index_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.index_page_path)

    # -----------------------------------------------------------------

    @property
    def has_status_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.status_page_path)

    # -----------------------------------------------------------------

    @property
    def status_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.status_page_path)

    # -----------------------------------------------------------------

    @property
    def data_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, data_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_data_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.data_page_path)

    # -----------------------------------------------------------------

    @property
    def data_page_name(self):

        """
        THisn function ...
        :return:
        """

        return fs.name(self.data_page_path)

    # -----------------------------------------------------------------

    @property
    def photometry_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, photometry_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_photometry_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.photometry_page_path)

    # -----------------------------------------------------------------

    @property
    def photometry_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.photometry_page_path)

    # -----------------------------------------------------------------

    @property
    def components_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, components_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_components_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.components_page_path)

    # -----------------------------------------------------------------

    @property
    def components_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.components_page_path)

    # -----------------------------------------------------------------

    @property
    def preparation_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, preparation_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_preparation_page(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.preparation_page_path)

    # -----------------------------------------------------------------

    @property
    def preparation_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.preparation_page_path)

    # -----------------------------------------------------------------

    @property
    def model_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, model_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_model_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.model_page_path)

    # -----------------------------------------------------------------

    @property
    def model_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.model_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, maps_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_maps_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.maps_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_page_path)

    # -----------------------------------------------------------------

    @property
    def fitting_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, fitting_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_fitting_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fitting_page_path)

    # -----------------------------------------------------------------

    @property
    def fitting_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.fitting_page_path)

    # -----------------------------------------------------------------

    @property
    def datacubes_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, datacubes_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_datacubes_page(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.datacubes_page_path)

    # -----------------------------------------------------------------

    @property
    def datacubes_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.datacubes_page_path)

    # -----------------------------------------------------------------

    @property
    def fluxes_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, fluxes_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_fluxes_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_page_path)

    # -----------------------------------------------------------------

    @property
    def fluxes_page_name(self):

        """
        This fnc
        :return:
        """

        return fs.name(self.fluxes_page_path)

    # -----------------------------------------------------------------

    @property
    def images_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, images_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_images_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.images_page_path)

    # -----------------------------------------------------------------

    @property
    def images_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.images_page_path)

    # -----------------------------------------------------------------

    @property
    def heating_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, heating_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_heating_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.heating_page_path)

    # -----------------------------------------------------------------

    @property
    def heating_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.heating_page_path)

    # -----------------------------------------------------------------

    @property
    def colours_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, colours_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_colours_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.colours_page_path)

    # -----------------------------------------------------------------

    @property
    def colours_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.colours_page_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, attenuation_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_attenuation_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.attenuation_page_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.attenuation_page_path)

    # -----------------------------------------------------------------

    @property
    def top_title(self):

        """
        This function ...
        :return:
        """

        # Create title
        return "High resolution 3D radiative transfer modeling of DustPedia galaxies"

    # -----------------------------------------------------------------

    @property
    def html_top_title(self):

        """
        This function ...
        :return:
        """

        return html.fontsize_template.format(size=top_title_size, text=html.underline_template.format(text=self.top_title))

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        #return "Modeling of " + self.galaxy_name
        return self.galaxy_name + " (" + self.ngc_name_nospaces + ")"

    # -----------------------------------------------------------------

    @property
    def html_title(self):

        """
        This function ...
        :return:
        """

        # Create title
        return html.fontsize_template.format(size=title_size, text=html.underline_template.format(text=self.title))

    # -----------------------------------------------------------------

    @property
    def heading(self):

        """
        This function ...
        :return:
        """

        heading = ""
        heading += self.html_top_title + html.newline + html.newline + html.line + html.newline
        heading += self.html_title + html.newline
        return heading

    # -----------------------------------------------------------------

    @property
    def footing(self):

        """
        This function ...
        :return:
        """

        text = html.small_template.format(text="Last updated " + time.pretty_time())

        footing = ""
        footing += html.newline + html.make_line("heavy") + html.center_template.format(text=text)
        return footing

    # -----------------------------------------------------------------

    @property
    def style(self):

        """
        This function ...
        :return:
        """

        return page_style

    # -----------------------------------------------------------------

    @property
    def tostr_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["round"] = True
        # kwargs["scientific"] = True  NO: let tostr decide
        kwargs["ndigits"] = 3
        kwargs["delimiter"] = ", "
        kwargs["html"] = True
        return kwargs

    # -----------------------------------------------------------------

    def make_page(self, body):

        """
        This function ...
        :param body:
        :return:
        """

        # Create contents
        contents = dict()
        contents["title"] = self.title
        contents["head"] = html.link_stylesheet_header_template.format(url=stylesheet_url)
        contents["body"] = body
        contents["style"] = self.style

        # Make page
        self.page = html.page_template.format(**contents)

    # -----------------------------------------------------------------

    @abstractproperty
    def page_path(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Write
        #fs.write_text(self.page_path, self.page)
        self.page.saveto(self.page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open
        browser.open_path(self.page_path)

# -----------------------------------------------------------------
