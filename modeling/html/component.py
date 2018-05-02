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
from ...core.tools.utils import lazyproperty
from ..truncation.component import ellipse_page_filename, significance_page_filename
from ..core.environment import all_maps_filename, maps_summary_filename, old_maps_filename, young_maps_filename, ionizing_maps_filename, dust_maps_filename, clip_maps_filename

# -----------------------------------------------------------------

# Stylesheets
stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"
slider_stylesheet_url = "http://users.ugent.be/~sjversto/slider.css"

# Scripts
sortable_url = "http://users.ugent.be/~sjversto/sorttable.js"
preview_url = "http://users.ugent.be/~sjversto/preview.js"
slider_url = "http://users.ugent.be/~sjversto/slider.js"

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
seds_page_filename = "seds.html"
datacubes_page_filename = "datacubes.html"
fluxes_page_filename = "fluxes.html"
images_page_filename = "images.html"
heating_page_filename = "heating.html"
colours_page_filename = "colours.html"
attenuation_page_filename = "attenuation.html"

# -----------------------------------------------------------------

class HTMLComponent(GalaxyModelingComponent):


    def __init__(self, *args, **kwargs):

        """
        Thisf ucntion ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HTMLComponent, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def has_properties(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("fetch_properties")

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("fetch_images")

    # -----------------------------------------------------------------

    @property
    def has_prepared(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("prepare_data")

    # -----------------------------------------------------------------

    @property
    def has_components(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("decompose")

    # -----------------------------------------------------------------

    @property
    def has_photometry(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("photometry")

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("build_model")

    # -----------------------------------------------------------------

    @property
    def has_fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("configure_fit")

    # -----------------------------------------------------------------

    @property
    def has_generation(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("fit_sed")

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        This function ...
        :return:
        """

        return "launch_analysis" in self.history
        #return self.history.is_finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_datacubes(self):

        """
        This function ...
        :return:
        """

        return "launch_analysis" in self.history
        #return self.history.is_finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_model_images(self):

        """
        This function ....
        :return:
        """

        return self.history.is_finished("launch_analysis")

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.history.has_finished_any("analyse_attenuation_map", "analyse_attenuation_curve")

    # -----------------------------------------------------------------

    @property
    def has_colours(self):

        """
        This function ...
        :return:
        """

        return self.history.is_finished("analyse_colours")

    # -----------------------------------------------------------------

    @property
    def has_heating(self):

        """
        This function ...
        :return:
        """

        return self.history.has_finished_any("analyse_cell_heating", "analyse_projected_heating")

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
    def seds_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, seds_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_seds_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.seds_page_path)

    # -----------------------------------------------------------------

    @property
    def seds_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, seds_page_filename)

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

    @lazyproperty
    def truncation_html_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.truncation_path, "html")

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.truncation_html_path, ellipse_page_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_significance_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.truncation_html_path, significance_page_filename)

    # -----------------------------------------------------------------

    @property
    def has_truncation_ellipse_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.truncation_ellipse_page_path)

    # -----------------------------------------------------------------

    @property
    def has_truncation_significance_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.truncation_significance_page_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_html_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "html")

    # -----------------------------------------------------------------

    @lazyproperty
    def all_maps_page_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.maps_html_path, all_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_summary_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, maps_summary_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_maps_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, old_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_maps_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, young_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_maps_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, ionizing_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_maps_page_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.maps_html_path, dust_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def clip_maps_page_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.maps_html_path, clip_maps_filename)

    # -----------------------------------------------------------------

    @property
    def has_all_maps_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.all_maps_page_path)

    # -----------------------------------------------------------------

    @property
    def has_maps_summary_page(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.maps_summary_page_path)

    # -----------------------------------------------------------------

    @property
    def has_old_maps_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_maps_page_path)

    # -----------------------------------------------------------------

    @property
    def has_young_maps_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_maps_page_path)

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps_page(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_maps_page_path)

    # -----------------------------------------------------------------

    @property
    def has_dust_maps_page(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.dust_maps_page_path)

    # -----------------------------------------------------------------

    @property
    def has_clip_maps_page(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.clip_maps_page_path)

# -----------------------------------------------------------------

class HTMLPageComponent(HTMLComponent):

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
    def scripts_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.environment.html_scripts_path

    # -----------------------------------------------------------------

    @lazyproperty
    def scripts_model_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.scripts_path, "model")

    # -----------------------------------------------------------------

    @lazyproperty
    def scripts_components_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.scripts_path, "components")

    # -----------------------------------------------------------------

    @property
    def images_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_images_path

    # -----------------------------------------------------------------

    @lazyproperty
    def images_data_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_path, "data")

    # -----------------------------------------------------------------

    @lazyproperty
    def images_maps_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_path, "maps")

    # -----------------------------------------------------------------

    @lazyproperty
    def images_seds_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_path, "seds")

    # -----------------------------------------------------------------

    @lazyproperty
    def images_datacubes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_path, "datacubes")

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

        # Fix
        fix_local_paths(self.page_path, self.html_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open
        with browser.serve_local_host(): browser.open_path(self.page_path)

# -----------------------------------------------------------------

def fix_local_paths(filepath, html_path):

    """
    This function ...
    :param filepath:
    :param html_path:
    :return:
    """

    filename = fs.name(filepath)

    # Debugging
    log.debug("Fixing local paths in the '" + filename + "' file ...")

    fixed = False
    newlines = []

    html_path = html_path + "/"

    # Loop over the lines in the file
    for line in fs.read_lines(filepath):

        if html_path in line:

            #print(line)
            line = line.replace(html_path, "")
            fixed = True
            #print("NEWLINE", line)

        # Add the line
        newlines.append(line)

    # Write, only if fixed
    if not fixed: return

    # Debugging
    log.debug("Overwriting the '" + filename + "' file with the new image paths ...")

    # Remove original file
    fs.remove_file(filepath)

    # Write the new lines
    fs.write_lines(filepath, newlines)

# -----------------------------------------------------------------
