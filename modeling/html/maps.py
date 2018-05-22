#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.maps Contains the MapsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class
from ...core.tools.utils import lazyproperty
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.tools.html import SimpleTable, HTMLPage
from ...core.basics.composite import SimplePropertyComposite
from ...magic.tools.info import get_image_info_strings, get_image_info
from ...core.tools import filesystem as fs
from ...core.tools import html
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

class MapsPageGenerator(HTMLPageComponent):

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
        super(MapsPageGenerator, self).__init__(*args, **kwargs)

        # The infos
        self.old_info = None
        self.young_info = None
        self.ionizing_info = None
        self.dust_info = None

        # The tables
        self.old_table = None
        self.young_table = None
        self.ionizing_table = None
        self.dust_table = None

        # Plots of the old stellar map
        self.old_plot = None
        self.old_deprojected_plot = None
        self.old_edgeon_plot = None

        # Plots of the young stellar map
        self.young_plot = None
        self.young_deprojected_plot = None
        self.young_edgeon_plot = None

        # Plots of the ionizing stellar map
        self.ionizing_plot = None
        self.ionizing_deprojected_plot = None
        self.ionizing_edgeon_plot = None

        # Plots of the dust map
        self.dust_plot = None
        self.dust_deprojected_plot = None
        self.dust_edgeon_plot = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get info
        self.get_info()

        # 2. Make tables
        self.make_tables()

        # 3. Make plots
        self.make_plots()

        # 4. Generaet the html
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
        super(MapsPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_model_definition(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_name(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_old_map_name_for_model(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_old_map_path(self.old_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.old_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojected_old_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_deprojected_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_old_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_edgeon_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_deprojected_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_old_skirt_deprojected_map_path(self.old_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_edgeon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_old_edgeon_map_path(self.old_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_name(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_young_map_name_for_model(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_young_map_path(self.young_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map(self):

        """
        Thisf unction ...
        :return:
        """

        return Frame.from_file(self.young_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.young_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojected_young_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.young_deprojected_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_young_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.young_edgeon_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_deprojected_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_young_skirt_deprojected_map_path(self.young_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_edgeon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_young_edgeon_map_path(self.young_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_name(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_ionizing_map_name_for_model(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_ionizing_map_path(self.ionizing_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.ionizing_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.ionizing_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojected_ionizing_map(self):

        """
        Thisf unction ...
        :return:
        """

        return Frame.from_file(self.ionizing_deprojected_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_ionizing_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.ionizing_edgeon_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_deprojected_map_path(self):

        """
        THis function ...
        :return:
        """

        return self.static_maps_selection.get_ionizing_skirt_deprojected_map_path(self.ionizing_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_edgeon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_ionizing_edgeon_map_path(self.ionizing_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_name(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.get_dust_map_name_for_model(self.config.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_dust_map_path(self.dust_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojected_dust_map(self):

        """
        Thisf unction ...
        :return:
        """

        return Frame.from_file(self.dust_deprojected_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_dust_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.dust_edgeon_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_deprojected_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_dust_skirt_deprojected_map_path(self.dust_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_edgeon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.static_maps_selection.get_dust_edgeon_map_path(self.dust_map_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def info_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["path"] = False
        kwargs["name"] = True
        kwargs["xsize"] = True
        kwargs["ysize"] = True
        kwargs["psf_filter"] = False
        kwargs["filesize"] = True
        kwargs["filter"] = False
        kwargs["wavelength"] = False
        kwargs["unit"] = False
        return kwargs

    # -----------------------------------------------------------------

    def get_info(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Getting map info ...")

        # Old stellar map
        self.get_old_info()

        # Young stellar map
        self.get_young_info()

        # Ionizing stellar map
        self.get_ionizing_info()

        # Dust map
        self.get_dust_info()

    # -----------------------------------------------------------------

    def get_old_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Get old stellar map info ...")

        # Get info
        info = get_image_info(self.old_map_name, self.old_map, **self.info_kwargs)

        # Make property composite
        self.old_info = SimplePropertyComposite.from_dict(info)

    # -----------------------------------------------------------------

    def get_young_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Get young stellar map info ...")

        # Get info
        info = get_image_info(self.young_map_name, self.young_map, **self.info_kwargs)

        # Make property composite
        self.young_info = SimplePropertyComposite.from_dict(info)

    # -----------------------------------------------------------------

    def get_ionizing_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Get ionizing stellar map info ...")

        # Get info
        info = get_image_info(self.ionizing_map_name, self.ionizing_map, **self.info_kwargs)

        # Make property composite
        self.ionizing_info = SimplePropertyComposite.from_dict(info)

    # -----------------------------------------------------------------

    def get_dust_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Get dust map info ...")

        # Get info
        info = get_image_info(self.dust_map_name, self.dust_map, **self.info_kwargs)

        # Make property composite
        self.dust_info = SimplePropertyComposite.from_dict(info)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make table of old stellar map
        self.make_old_table()

        # Make table of young stellar map
        self.make_young_table()

        # Make table of ionizing stellar map
        self.make_ionizing_table()

        # Make table of dust map
        self.make_dust_table()

    # -----------------------------------------------------------------

    def make_old_table(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making table for the old stellar map ...")

        # Make the table
        self.old_table = SimpleTable.from_composite(self.old_info, css_class=table_class)

    # -----------------------------------------------------------------

    def make_young_table(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making table of the young stellar map ...")

        # Make the table
        self.young_table = SimpleTable.from_composite(self.young_info, css_class=table_class)

    # -----------------------------------------------------------------

    def make_ionizing_table(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making table of the ionizing stellar map ...")

        # Make the table
        self.ionizing_table = SimpleTable.from_composite(self.ionizing_info, css_class=table_class)

    # -----------------------------------------------------------------

    def make_dust_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making table of the dust map ...")

        # Make the table
        self.dust_table = SimpleTable.from_composite(self.dust_info, css_class=table_class)

    # -----------------------------------------------------------------

    @property
    def image_width(self):

        """
        This function ...
        :return:
        """

        return 400

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make old stellar map plots
        self.make_old_plots()

        # Make young stellar map plots
        self.make_young_plots()

        # Make ionizing stellar map plots
        self.make_ionizing_plots()

        # Make dust map plots
        self.make_dust_plots()

    # -----------------------------------------------------------------

    def make_old_plots(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the old stellar map ...")

        # Determine the paths
        map_path = fs.join(self.images_maps_path, "old_map.png")
        deprojected_map_path = fs.join(self.images_maps_path, "old_map_deprojected.png")
        edgeon_map_path = fs.join(self.images_maps_path, "old_map_edgeon.png")

        # Plot the map
        if not fs.is_file(map_path):
            vmin, vmax = self.old_map.saveto_png(map_path, colours=self.config.colours, interval=self.config.interval,
                                    scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)
        #interval = [vmin, vmax]
        interval = self.config.interval

        if not fs.is_file(deprojected_map_path):
            # Plot the deprojected map
            self.deprojected_old_map.saveto_png(deprojected_map_path, colours=self.config.colours, interval=interval,
                                                scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        if not fs.is_file(edgeon_map_path):
            # Plot the edgeon map
            self.edgeon_old_map.saveto_png(edgeon_map_path, colours=self.config.colours, interval=interval,
                                           scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Make images
        self.old_plot = html.image(map_path, alttext=self.old_map_name, width=self.image_width)
        self.old_deprojected_plot = html.image(deprojected_map_path, alttext=self.old_map_name, width=self.image_width)
        self.old_edgeon_plot = html.image(edgeon_map_path, alttext=self.old_map_name, width=self.image_width)

    # -----------------------------------------------------------------

    def make_young_plots(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the young stellar map ...")

        # Determine the paths
        map_path = fs.join(self.images_maps_path, "young_map.png")
        deprojected_map_path = fs.join(self.images_maps_path, "young_map_deprojected.png")
        edgeon_map_path = fs.join(self.images_maps_path, "young_map_edgeon.png")

        # Plot the map
        if not fs.is_file(map_path):
            vmin, vmax = self.young_map.saveto_png(map_path, colours=self.config.colours, interval=self.config.interval,
                                                   scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)
        #interval = [vmin, vmax]
        interval = self.config.interval

        if not fs.is_file(deprojected_map_path):
            # Plot the deprojected map
            self.deprojected_young_map.saveto_png(deprojected_map_path, colours=self.config.colours, interval=interval,
                                                  scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        if not fs.is_file(edgeon_map_path):
            # Plot the edgeon map
            self.edgeon_young_map.saveto_png(edgeon_map_path, colours=self.config.colours, interval=interval,
                                             scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Make images
        self.young_plot = html.image(map_path, alttext=self.young_map_name, width=self.image_width)
        self.young_deprojected_plot = html.image(deprojected_map_path, alttext=self.young_map_name, width=self.image_width)
        self.young_edgeon_plot = html.image(edgeon_map_path, alttext=self.young_map_name, width=self.image_width)

    # -----------------------------------------------------------------

    def make_ionizing_plots(self):

        """
        Thisnfunction ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the ionizing stellar map ...")

        # Determine the paths
        map_path = fs.join(self.images_maps_path, "ionizing_map.png")
        deprojected_map_path = fs.join(self.images_maps_path, "ionizing_map_deprojected.png")
        edgeon_map_path = fs.join(self.images_maps_path, "ionizing_map_edgeon.png")

        # Plot the map
        if not fs.is_file(map_path):
            vmin, vmax = self.ionizing_map.saveto_png(map_path, colours=self.config.colours, interval=self.config.interval,
                                                      scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)
        #interval = [vmin, vmax]
        interval = self.config.interval

        if not fs.is_file(deprojected_map_path):
            # Plot the deprojected map
            self.deprojected_ionizing_map.saveto_png(deprojected_map_path, colours=self.config.colours, interval=interval,
                                                     scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        if not fs.is_file(edgeon_map_path):
            # Plot the edgeon map
            self.edgeon_ionizing_map.saveto_png(edgeon_map_path, colours=self.config.colours, interval=interval,
                                                scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Make images
        self.ionizing_plot = html.image(map_path, alttext=self.ionizing_map_name, width=self.image_width)
        self.ionizing_deprojected_plot = html.image(deprojected_map_path, alttext=self.ionizing_map_name, width=self.image_width)
        self.ionizing_edgeon_plot = html.image(edgeon_map_path, alttext=self.ionizing_map_name, width=self.image_width)

    # -----------------------------------------------------------------

    def make_dust_plots(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the dust map ...")

        # Determine the paths
        map_path = fs.join(self.images_maps_path, "dust_map.png")
        deprojected_map_path = fs.join(self.images_maps_path, "dust_map_deprojected.png")
        edgeon_map_path = fs.join(self.images_maps_path, "dust_map_edgeon.png")

        # Plot the map
        if not fs.is_file(map_path):
            vmin, vmax = self.dust_map.saveto_png(map_path, colours=self.config.colours, interval=self.config.interval,
                                              scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)
        #interval = [vmin, vmax]
        interval = self.config.interval

        if not fs.is_file(deprojected_map_path):
            # Plot the deprojected map
            self.deprojected_dust_map.saveto_png(deprojected_map_path, colours=self.config.colours, interval=interval,
                                                 scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        if not fs.is_file(edgeon_map_path):
            # Plot the edgeon map
            self.edgeon_dust_map.saveto_png(edgeon_map_path, colours=self.config.colours, interval=interval,
                                            scale=self.config.scale, alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Make images
        self.dust_plot = html.image(map_path, alttext=self.dust_map_name, width=self.image_width)
        self.dust_deprojected_plot = html.image(deprojected_map_path, alttext=self.dust_map_name, width=self.image_width)
        self.dust_edgeon_plot = html.image(edgeon_map_path, alttext=self.dust_map_name, width=self.image_width)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Genreate the page
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
        #javascript_paths.append(sortable_url)
        #javascript_paths.append(preview_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths, footing=updated_footing())

        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        self.page += "OLD STELLAR MAP"
        self.page += html.newline + html.newline
        self.page += self.old_table
        self.page += html.newline + html.newline
        self.page += self.old_plot
        self.page += html.newline + html.newline
        self.page += "Deprojected:" + html.newline
        self.page += self.old_deprojected_plot
        self.page += html.newline + html.newline
        self.page += "Edge-on:" + html.newline
        self.page += self.old_edgeon_plot
        self.page += html.newline + html.newline

        self.page += html.line
        self.page += html.newline + html.newline

        self.page += "YOUNG STELLAR MAP"
        self.page += html.newline + html.newline
        self.page += self.young_table
        self.page += html.newline + html.newline
        self.page += self.young_plot
        self.page += html.newline + html.newline
        self.page += "Deprojected:" + html.newline
        self.page += self.young_deprojected_plot
        self.page += html.newline + html.newline
        self.page += "Edge-on:" + html.newline
        self.page += self.young_edgeon_plot
        self.page += html.newline + html.newline

        self.page += html.line
        self.page += html.newline + html.newline

        self.page += "IONIZING STELLAR MAP"
        self.page += html.newline + html.newline
        self.page += self.ionizing_table
        self.page += html.newline + html.newline
        self.page += self.ionizing_plot
        self.page += html.newline + html.newline
        self.page += "Deprojected:" + html.newline
        self.page += self.ionizing_deprojected_plot
        self.page += html.newline + html.newline
        self.page += "Edge-on:" + html.newline
        self.page += self.ionizing_edgeon_plot
        self.page += html.newline + html.newline

        self.page += html.line
        self.page += html.newline + html.newline

        self.page += "DUST MAP"
        self.page += html.newline + html.newline
        self.page += self.dust_table
        self.page += html.newline + html.newline
        self.page += self.dust_plot
        self.page += html.newline + html.newline
        self.page += "Deprojected:" + html.newline
        self.page += self.dust_deprojected_plot
        self.page += html.newline + html.newline
        self.page += "Edge-on:" + html.newline
        self.page += self.dust_edgeon_plot
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

        return self.maps_page_path

# -----------------------------------------------------------------
