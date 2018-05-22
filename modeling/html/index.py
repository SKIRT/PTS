#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.index Contains the IndexPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class, hover_table_class
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ..core.progression import create_modeling_progression

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

class IndexPageGenerator(HTMLPageComponent):

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
        super(IndexPageGenerator, self).__init__(*args, **kwargs)

        # The modeling progression
        self.progression = None

        # Tables
        self.info_table = None

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
        super(IndexPageGenerator, self).setup(**kwargs)

        # Get the modeling progression
        if "progression" in kwargs: self.progression = kwargs.pop("progression")
        else: self.progression = create_modeling_progression(self.config.path)

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

    # -----------------------------------------------------------------

    def make_info_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the info table ...")

        # Create list of tuples
        tuples = self.galaxy_info.items()
        tuples.append(("Distance", self.galaxy_properties.distance))
        tuples.append(("Redshift", self.galaxy_properties.redshift))
        tuples.append(("Center coordinate", self.galaxy_properties.center))
        tuples.append(("Major axis length", self.galaxy_properties.major))
        tuples.append(("Ellipticity", self.galaxy_properties.ellipticity))
        tuples.append(("Position angle", self.galaxy_properties.position_angle))
        tuples.append(("Inclination", self.galaxy_properties.inclination))

        # Create the table
        self.info_table = html.SimpleTable(tuples, header=["Property", "Value"], css_class=hover_table_class, tostr_kwargs=self.tostr_kwargs)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

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
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Make javascripts urls
        javascript_paths = javascripts[:]
        javascript_paths.append(sortable_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        #classes = dict()
        #classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        #body = self.heading

        #body += html.newline + html.line
        self.page += html.newline + "Pages:" + html.newline

        items = []

        # Status page
        items.append(html.hyperlink(self.status_page_name, "status") + ": status of the modeling")

        # Data page
        if self.has_images: items.append(html.hyperlink(self.data_page_name, "data") + ": summary of the observational dataset")

        # Preparation page
        if self.has_prepared: items.append(html.hyperlink(self.preparation_page_name, "preparation") + ": image preparation statistics")

        # Components page
        if self.has_components: items.append(html.hyperlink(self.components_page_name, "components") + ": decomposition details")

        # Photometry page
        if self.has_photometry: items.append(html.hyperlink(self.photometry_page_name, "photometry") + ": observed SED")

        # Maps page
        #if self.has_model: items.append(html.hyperlink(self.maps_page_name, "maps"))
        if self.progression.model_name is not None: items.append(html.hyperlink(self.maps_page_name, "maps") + ": distributions used as input geometries for the model")

        # Model page
        #if self.has_fitting_run: items.append(html.hyperlink(self.model_page_name, "model"))
        if self.progression.model_name is not None: items.append(html.hyperlink(self.model_page_name, "model") + ": components and properties of the model")

        # Fitting page
        if self.has_generation: items.append(html.hyperlink(self.fitting_page_name, "fitting") + ": fitting details")

        # SEDs page
        if self.has_seds: items.append(html.hyperlink(self.seds_page_name, "SEDs") + ": output SEDs and their different contributions")

        # Datacubes page
        if self.has_datacubes: items.append(html.hyperlink(self.datacubes_page_name, "datacubes") + ": output datacubes")

        # Fluxes page
        if self.has_fluxes: items.append(html.hyperlink(self.fluxes_page_name, "fluxes") + ": calculated observed fluxes from simulation")

        # Images page
        if self.has_model_images: items.append(html.hyperlink(self.images_page_name, "images") + ": observed images from simulation")

        # Attenuation page
        if self.has_attenuation: items.append(html.hyperlink(self.attenuation_page_name, "attenuation") + ": dust attenuation analysis output")

        # Colours page
        if self.has_colours: items.append(html.hyperlink(self.colours_page_name, "colours") + ": colours calculated based on simulated images")

        # Heating page
        if self.has_heating: items.append(html.hyperlink(self.heating_page_name, "heating") + ": dust heating analysis output")

        # Add the list
        self.page += html.unordered_list(items, css_class="b")
        self.page += html.line + html.newline

        # Add links to detailed pages
        if self.config.details:

            items = []

            self.page += html.newline + "Detailed pages:" + html.newline

            # TRUNCATION

            if self.has_truncation_ellipse_page: items.append(html.hyperlink(self.truncation_ellipse_page_path, "truncation ellipse") + ": to determine the boundary at which all image data is truncated")
            if self.has_truncation_significance_page: items.append(html.hyperlink(self.truncation_significance_page_path, "significance levels") + ": to determine the minimal signal-to-noise for each image to compare with simulation")

            # MAPS

            if self.has_all_maps_page: items.append(html.hyperlink(self.all_maps_page_path, "all generated maps") + ": collection of all possible maps that have been generated for visual inspection")
            if self.has_maps_summary_page: items.append(html.hyperlink(self.maps_summary_page_path, "maps summary") + ": summary of map properties")
            if self.has_old_maps_page: items.append(html.hyperlink(self.old_maps_page_path, "old stellar maps") + ": summary of the old stellar maps")
            if self.has_young_maps_page: items.append(html.hyperlink(self.young_maps_page_path, "young stellar maps") + ": summary of the young stellar maps")
            if self.has_ionizing_maps_page: items.append(html.hyperlink(self.ionizing_maps_page_path, "ionizing stellar maps") + ": summary of the ionizing stellar maps")
            if self.has_dust_maps_page: items.append(html.hyperlink(self.dust_maps_page_path, "dust maps") + ": summary of the dust maps")
            if self.has_clip_maps_page: items.append(html.hyperlink(self.clip_maps_page_path, "map clipping") + ": to determine the minimal signal-to-noise levels at which to clip each input image to obtain the final map pixels")

            # Add the list
            self.page += html.unordered_list(items, css_class="b")
            self.page += html.line + html.newline

        # Add the galaxy info table
        title_info = html.underline_template.format(text="Galaxy info")
        self.page += title_info + html.newline + html.newline + str(self.info_table) + html.newline

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

        return self.index_page_path

# -----------------------------------------------------------------
