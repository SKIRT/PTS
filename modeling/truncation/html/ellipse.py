#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.html.ellipse Contains the TruncationEllipsePageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc

# Import the relevant PTS classes and modules
from ..component import TruncationComponent
from ....core.tools.html import HTMLPage, SimpleTable, newline, updated_footing, make_theme_button, center, make_script_button, dictionary, make_usable
from ....core.basics.log import log
from ....magic.view.html import JS9Preloader, body_settings, javascripts, css_scripts, JS9Loader, JS9Spawner, JS9Window, make_load_region, make_synchronize_regions_script
from ....magic.view.html import make_spawn_code, add_to_div
from ....core.tools import filesystem as fs
from ....core.tools.utils import lazyproperty
from ....core.filter.filter import parse_filter
from ..analytics import mask_names
from ....core.tools import browser
from ....magic.core.frame import Frame
from ....magic.tools.info import get_image_info_from_header
from ....magic.region.list import PixelRegionList

# -----------------------------------------------------------------

ncolumns = 3
image_width = 300
image_height = 300
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

base_url = "http://users.ugent.be/~sjversto"
stylesheet_filename = "stylesheet.css"
stylesheet_url = fs.join(base_url, stylesheet_filename)

style = "ugentstyle"

# -----------------------------------------------------------------

class TruncationEllipsePageGenerator(TruncationComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(TruncationEllipsePageGenerator, self).__init__(*args, **kwargs)

        # --- Attributes ---

        # The truncation factor indicator
        self.indicator = None

        # The preloader
        self.preloader = None

        # The loaders
        self.loaders = dict()

        # The region loaders
        self.region_loaders = dict()

        # The windows
        self.windows = dict()

        # The page
        self.page = None

        # The plots directory
        self.plots_path = None

        # Plot paths
        self.plots_paths = dict()

        # The masks directory
        self.masks_path = None

        # The ellipses path
        self.ellipses_path = None

        # Mask paths
        self.mask_paths = dict()

        # Info
        self.info = dict()

        # The coordinate systems
        self.coordinate_systems = dict()

        # The display IDs
        self.display_ids = dict()

        # The ellipses
        self.ellipses = dict()

        # All images loader
        self.all_loader = None

        # The table
        self.table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the masks
        if self.config.mask: self.load_masks()

        # 3. Make the eindicator
        self.make_indicator()

        # 4. Make info
        if self.config.info: self.get_info()

        # 5. Set the coordinate systems
        self.set_coordinate_systems()

        # 6. Make plots
        self.make_plots()

        # 7. Create the regions, for the coordinate systems
        if not self.has_ellipses: self.create_ellipses()
        else: self.load_ellipses()

        # 8. Make the views
        self.make_views()

        # 9. Make table
        self.make_table()

        # 10. Generate the page
        self.generate_page()

        # 11. Writing
        self.write()

        # 12. Showing
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TruncationEllipsePageGenerator, self).setup(**kwargs)

        # Set the preloader
        #if self.config.preload: self.preloader = JS9Preloader()
        self.preloader = JS9Preloader()

        # Make directory to contain the plots
        self.plots_path = fs.join(self.truncation_html_path, "plots")
        if fs.is_directory(self.plots_path):
            if self.config.replot: fs.clear_directory(self.plots_path)
        else: fs.create_directory(self.plots_path)

        # Make directory to contain the mask plots
        self.masks_path  = fs.join(self.truncation_html_path, "masks")
        if fs.is_directory(self.masks_path):
            if self.config.replot: fs.clear_directory(self.masks_path)
        else: fs.create_directory(self.masks_path)

        # Make directory to contain the ellipses in pixel coordinates
        self.ellipses_path = fs.join(self.truncation_html_path, "ellipses")
        if fs.is_directory(self.ellipses_path):
            if self.config.replot: fs.clear_directory(self.ellipses_path)
        else: fs.create_directory(self.ellipses_path)

        # Check
        if self.config.reproject:
            if not self.config.png: raise ValueError("The 'png' option has to be enabled for the reprojection")

        # Cannot enable reprojection AND downsampling (one or the other)
        if self.config.downsample:
            if self.config.reproject: raise ValueError("Cannot use both reprojection and downsampling")
            if not self.config.png: raise ValueError("The 'png' option has to be enabled for downsampling")

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        # Sorted on wavelength!
        if self.config.filters is not None: return sorted(self.config.filters, key=lambda fltr: fltr.wavelength.to("micron").value)
        #else: return sorted(self.dataset.filters, key=lambda fltr: fltr.wavelength) # TOO SLOW
        else: return sorted((parse_filter(name) for name in self.dataset.names), key=lambda fltr: fltr.wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Truncation"

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        #return [self.dataset.get_name_for_filter(fltr) for fltr in self.filters]
        return self.dataset.get_names_for_filters(self.filters)

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Loop over the images
        for name, fltr in zip(self.names, self.filters):

            # Get the mask
            mask = self.dataset.get_image_masks_union(name, mask_names, strict=False)

            # Save the mask as a PNG image
            frame = Frame(mask.data.astype(int), wcs=mask.wcs)

            # Determine path
            filepath = fs.join(self.masks_path, name + ".png")

            # Set the path
            self.mask_paths[name] = filepath

            # Collect
            gc.collect()

            if fs.is_file(filepath): continue  # plot already there

            # Save as PNG
            frame.saveto_png(filepath, colours="grey")

    # -----------------------------------------------------------------

    @property
    def indicator_id(self):

        """
        This function ...
        :return:
        """

        return "indicator"

    # -----------------------------------------------------------------

    @property
    def default_factor(self):

        """
        This function ...
        :return:
        """

        # There is already a truncation ellipse
        if self.has_truncation_ellipse: return self.truncation_factor

        # There is not yet a truncation ellipse
        else: return 1.0

    # -----------------------------------------------------------------

    def make_indicator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the truncation factor indicator ...")

        # Display factor
        self.indicator = "<div id='" + self.indicator_id + "'>\n"
        self.indicator += "Factor: " + str(self.default_factor) + "\n"
        self.indicator += "</div>"

    # -----------------------------------------------------------------

    def get_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info for the images ...")

        # Loop over the images
        for name, fltr in zip(self.names, self.filters):

            # Get header
            header = self.dataset.get_header(name)

            # Get the image info
            path = self.dataset.get_frame_path(name)
            info = get_image_info_from_header(name, header, image_path=path, path=False, name=False)

            # Make list
            code = dictionary(info, key_color=key_color)

            # Add
            self.info[name] = code

    # -----------------------------------------------------------------

    def set_coordinate_systems(self):

        """
        This fnuction ...
        :return:
        """

        # Inform the suer
        log.info("Setting the coordinate systems ...")

        # Determine rebinning wcs
        if self.config.reproject:

            # Which WCS to use?
            if self.config.reproject_method == "max": rebin_name = self.dataset.max_pixelscale_name
            elif self.config.reproject_method == "median": rebin_name = self.dataset.median_pixelscale_name
            elif self.config.reproject_method == "largest": rebin_name = self.dataset.largest_wcs_name
            elif self.config.reproject_method == "largest_from_median": rebin_name = self.dataset.largest_wcs_below_median_pixelscale_name
            elif self.config.reproject_method == "closest_pixelscale": rebin_name = self.dataset.get_closest_pixelscale_name(self.config.reproject_pixelscale)
            else: raise ValueError("Invalid reproject method: '" + self.config.reproject_method + "'")

            # Debugging
            log.debug("Rebinning all images to the coordinate system of the '" + rebin_name + "' image ...")
            rebin_wcs = self.dataset.get_wcs(rebin_name)

        # Don't rebin
        else:
            log.debug("Not rebinning the images")
            rebin_wcs = None

        # Loop over the images
        # for name in self.dataset.names:
        for name, fltr in zip(self.names, self.filters):

            # Get frame
            frame = self.dataset.get_frame(name)

            # Set the WCS
            if rebin_wcs is not None: self.coordinate_systems[name] = rebin_wcs
            else: self.coordinate_systems[name] = frame.wcs

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make image plots
        if self.config.png: self.make_image_plots()

    # -----------------------------------------------------------------

    def make_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making image plots ...")

        # Loop over the images
        for name, fltr in zip(self.names, self.filters):

            # Determine path
            filepath = fs.join(self.plots_path, name + ".png")

            # Set the path
            self.plots_paths[name] = filepath

            # Is the plot already made?
            if fs.is_file(filepath):
                log.success("The PNG for the " + name + " frame is already created")
                continue  # plot already there

            # Debugging
            log.debug("Loading the '" + name + "' frame ...")

            # Get frame
            frame = self.dataset.get_frame(name)

            # Get the wcs
            rebin_wcs = self.coordinate_systems[name]

            # REPROJECT
            if rebin_wcs is not None: frame.rebin(rebin_wcs, exact=False, convert=True)

            # Downsample?
            if self.config.downsample:
                downsample_factor = frame.downsample_to_npixels(self.config.downsample_npixels_threshold)
                if downsample_factor != 1: self.coordinate_systems[name] = frame.wcs

            # Debugging
            log.debug("Making plot of the '" + name + "' image ...")

            # Save
            #frame.saveto_png(filepath, colours=self.config.colormap, alpha="absolute")
            frame.saveto_png(filepath, scale="linear", colours="gray", alpha="absolute") # PNG IN GRAYSCALE

            # Cleanup
            gc.collect()

    # -----------------------------------------------------------------

    @property
    def has_all_plots(self):

        """
        Thisf unction ...
        :return:
        """

        for name in self.names:

            # Determine path
            filepath = fs.join(self.plots_path, name + ".png")

            if not fs.is_file(filepath): return False

        return True

    # -----------------------------------------------------------------

    @property
    def has_some_plots(self):

        """
        This function ...
        :return:
        """

        for name in self.names:

            # Determine path
            filepath = fs.join(self.plots_path, name + ".png")

            if fs.is_file(filepath): return True

        return False

    # -----------------------------------------------------------------

    @property
    def has_ellipses(self):

        """
        Thsif unction ...
        :return:
        """

        # Regenerate?
        if self.config.regenerate_ellipses:
            fs.clear_directory(self.ellipses_path)
            return False

        # Loop over all prepared images, get the images
        for name in self.names:
            if not self.has_ellipse_for_name(name): return False

        if self.has_truncation_ellipse: log.warning("Ellipses present from previous run in ellipses/ directory are used but truncation factor is already determined: make sure the regions correspond to the truncation radius (add --regenerate_ellipses if necessary)")

        return True

    # -----------------------------------------------------------------

    def has_ellipse_for_name(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return fs.is_file(self.ellipse_path_for_name(name))

    # -----------------------------------------------------------------

    def ellipse_path_for_name(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return fs.join(self.ellipses_path, name + ".reg")

    # -----------------------------------------------------------------

    def create_ellipses(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the ellipses ...")

        # Warning
        if self.has_some_plots: log.warning("Ellipses are created but some plots are already made: coordinate systems may deviate. Run with --replot to solve issues.")

        # Loop over all prepared images, get the images
        for name, fltr in zip(self.names, self.filters):

            # Debugging
            log.debug("Creating the ellipse for the '" + name + "' image ...")

            # Get region in image coordinates
            region = self.disk_ellipse.to_pixel(self.coordinate_systems[name])

            # Add the region
            self.ellipses[name] = region

    # -----------------------------------------------------------------

    def load_ellipses(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Loading the ellipses ...")

        # Loop over the names
        for name in self.names:

            path = self.ellipse_path_for_name(name)
            self.ellipses[name] = PixelRegionList.from_file(path)[0]

    # -----------------------------------------------------------------

    def make_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views ...")

        load_info = dict()
        images = dict()
        region_loads = dict()
        placeholders = dict()
        views = dict()

        # Loop over all prepared images, get the images
        for name, fltr in zip(self.names, self.filters):

            # Get name
            display_name = make_usable(name)

            # Get path
            if self.config.png: path = self.plots_paths[name]
            else: path  = self.dataset.get_frame_path(name)

            settings = dict()
            settings["scale"] = self.config.scale
            settings["colormap"] = self.config.colormap
            settings["zoom"] = self.config.zoom
            #settings["fits2png"] = "true"

            #regions = "ellipse"

            # Get the angle
            #center = self.disk_ellipse.center  # in sky coordinates
            #semimajor = self.disk_ellipse.semimajor
            #semiminor = self.disk_ellipse.semiminor
            #angle = self.disk_ellipse.angle

            # Determine the ratio of semimajor and semiminor
            #ratio = semiminor / semimajor

            #region_string = str(self.disk_ellipse)

            # Get region
            region = self.ellipses[name]
            regions_for_loader = self.ellipses[name] if self.config.load_regions else None

            # Add preload
            if self.config.preload_all or (self.config.preload is not None and fltr in self.config.preload):

                # Add to preloader
                image = self.preloader.add_path(name, path, settings=settings, display=display_name, regions=regions_for_loader)

                # Create window
                self.windows[name] = JS9Window(display_name, width=image_width, height=image_height,
                                               background_color=background_color, menubar=self.config.menubar,
                                               colorbar=self.config.colorbar, resize=self.config.resize)
                display_id = display_name

                # Set load info
                views[display_id] = self.windows[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = image

            # Add dynamic load
            elif self.config.dynamic:

                # Create the loader
                self.loaders[name] = JS9Spawner.from_path("Load image", name, path, settings=settings, button=True,
                                                          menubar=self.config.menubar, colorbar=self.config.colorbar,
                                                          regions=regions_for_loader, add_placeholder=False, background_color=background_color)
                display_id = self.loaders[name].display_id

                self.windows[name] = self.loaders[name].placeholder

                # Set load info
                views[display_id] = self.loaders[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = self.loaders[name].image
                placeholders[display_id] = self.loaders[name].spawn_div_name

            # Regular button load in a pre-existing viewer
            else:

                # Create loader
                self.loaders[name] = JS9Loader.from_path("Load image", name, path, display=display_name,
                                                         settings=settings, button=True, regions=regions_for_loader)

                # Create window
                self.windows[name] = JS9Window(display_name, width=image_width, height=image_height, background_color=background_color, menubar=self.config.menubar, colorbar=self.config.colorbar, resize=self.config.resize)

                display_id = display_name

                # Set load info
                views[display_id] = self.windows[name].view
                load_info[display_id] = (name, path, regions_for_loader)
                images[display_id] = self.loaders[name].image

            # Set display ID
            self.display_ids[name] = display_id

            # Regions button
            region_button_id = display_name + "regionsbutton"
            load_region_function_name = "load_regions_" + display_name
            # load_region = make_load_region_function(load_region_function_name, regions, display=None)
            load_region = make_load_region(region, display=display_id, movable=False, rotatable=False,
                                           removable=False, resizable=True, quote_character="'")

            region_loads[display_id] = load_region

            # Create region loader
            self.region_loaders[name] = make_script_button(region_button_id, "Load regions", load_region,
                                                           load_region_function_name)

            # CREATE ALL IMAGES LOADER
            #buttonid = self.image.name + "Loader"
            #load_html = self.image.load(regions=self.regions)
            #return html.button(buttonid, self.text, load_html, quote_character=strings.other_quote_character(self.text, load_html))

        all_loader_name = "allimagesloaderbutton"
        all_loader_text = "Load all images"

        load_script = ""

        # Load over the images (displays)
        for display_id in load_info:

            #name, path, regions_for_loader = load_info[display_id]
            image = images[display_id]
            view = views[display_id]

            load_image = image.load()
            load_region = region_loads[display_id]

            if display_id in placeholders:

                spawn_div_name = placeholders[display_id]

                # Make spawn code
                spawn_code = make_spawn_code(view, image, menubar=self.config.menubar, colorbar=self.config.colorbar, width=image_width,
                                             background_color=background_color)

                # Add html code to DIV
                load_script += add_to_div(spawn_div_name, spawn_code)

                # Add DIV to JS9
                #load_script += "\n"
                load_script += "JS9.AddDivs('" + display_id + "');\n"

            # Load image and region code
            load_script += load_image
            load_script += "\n"
            load_script += load_region
            load_script += "\n\n"

        function_name = "loadAllImages"
        self.all_loader = make_script_button(all_loader_name, all_loader_text, load_script, function_name)

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        if self.config.dynamic:
            if self.config.info: return 2
            else: return 1
        else: return ncolumns

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "realtable"

    # -----------------------------------------------------------------

    def make_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making table ...")

        cells = []

        # Loop over the images
        for name, fltr in zip(self.names, self.filters):

            # Add name of the image
            string = center(name + newline)

            # Add loader if necessary
            if name in self.loaders:
                string += center(str(self.loaders[name]))
                string += newline

            # Add region loader if necessary
            if name in self.region_loaders:
                string += center(str(self.region_loaders[name]))
                string += newline

            # Add the window
            if name in self.windows: string += str(self.windows[name])

            # Add the cell
            cells.append(string)

            # Add info
            if self.config.info: cells.append(self.info[name])

        # Make the table
        self.table = SimpleTable.rasterize(cells, ncolumns=self.ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        settings_body = body_settings if self.preloader.has_images else None

        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = HTMLPage(self.title, body_settings=settings_body, javascript_path=javascripts, css_path=css_paths, style=style, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"

        ellipses = dict()
        for name in self.ellipses:
            display_id = self.display_ids[name]
            ellipses[display_id] = self.ellipses[name]

        #ellipse = self.ellipses[name]
        self.page += make_synchronize_regions_script(self.indicator_id, self.display_ids.values(), ellipses, start_factor=self.default_factor) + "\n"

        # Add theme button
        self.page += center(make_theme_button(classes=classes)) + newline

        # Add the indicator
        self.page += center(self.indicator) + newline

        # Add all images loader
        if self.all_loader is not None: self.page += center(str(self.all_loader)) + newline

        # Add the table
        self.page += self.table

        # Add preloader
        if self.preloader.has_images: self.page += self.preloader

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ellipses
        self.write_ellipses()

        # Write the page
        self.write_page()

    # -----------------------------------------------------------------

    def write_ellipses(self):

        """
        Thisf unction ...
        :return:
        """

        # INform the user
        log.info("Writing the ellipses ...")

        # Loop over the ellipses
        for name in self.ellipses:

            # Determine path
            path = self.ellipse_path_for_name(name)

            # Write
            self.ellipses[name].saveto(path)

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Save
        self.page.saveto(self.ellipse_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.ellipse_page_path)

# -----------------------------------------------------------------
