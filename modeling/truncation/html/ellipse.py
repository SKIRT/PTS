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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..component import TruncationComponent
from ....core.tools.html import newline, center, dictionary
from ....core.basics.log import log
from ....magic.view.html import make_synchronize_regions_script
from ....core.tools import filesystem as fs
from ....core.tools.utils import lazyproperty
from ..analytics import mask_names
from ....core.tools import browser
from ....magic.core.frame import Frame
from ....magic.tools.info import get_image_info_from_header
from ....magic.region.list import PixelRegionList
from ....magic.view.multi import MultiImageViewer
from ....core.filter.broad import BroadBandFilter

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

        # The ellipses
        self.ellipses = dict()

        # The viewer
        self.viewer = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # 9. Generate the page
        self.generate_page()

        # 10. Writing
        self.write()

        # 11. Showing
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
        else: return sorted((fltr for fltr in self.dataset.filters_from_names if not (isinstance(fltr, BroadBandFilter) and fltr.is_planck)), key=lambda fltr: fltr.wavelength.to("micron").value)

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
            #frame.saveto_png(filepath, scale="linear", colours="gray", alpha="absolute") # PNG IN GRAYSCALE
            frame.saveto_png(filepath, scale=self.config.scale, colours="gray", alpha="absolute") # USE DEDICATED SCALE, USE LINEAR IN THE VIEWER

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

    @lazyproperty
    def image_paths(self):

        """
        Thisf unction ...
        :return:
        """

        # Initialize
        paths = OrderedDict()

        # Loop over all prepared images, get the images
        for name, fltr in zip(self.names, self.filters):

            # Get path
            if self.config.png: path = self.plots_paths[name]
            else: path = self.dataset.get_frame_path(name)

            # Determine the relative path
            relpath = fs.relative_to(path, self.truncation_html_path)

            # Set the path
            paths[name] = relpath

        # Return the dictionary
        return paths

    # -----------------------------------------------------------------

    def make_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the views ...")

        # Initialize the viewer
        self.viewer = MultiImageViewer()

        # Set the scale
        if self.config.png: self.viewer.config.scale = "linear"
        else: self.viewer.config.scale = self.config.scale

        # Set colormap and zoom
        self.viewer.config.colormap = self.config.colormap
        self.viewer.config.zoom = self.config.zoom

        # Set advanced settings
        self.viewer.config.preload_all = self.config.preload_all
        self.viewer.config.preload = self.config.preload
        self.viewer.config.dynamic = self.config.dynamic
        self.viewer.config.load_regions = self.config.load_regions

        # Set regions settings
        self.viewer.config.regions.color = "white"
        self.viewer.config.regions.changeable = False
        self.viewer.config.regions.movable = False
        self.viewer.config.regions.rotatable = False
        self.viewer.config.regions.removable = False
        self.viewer.config.regions.resizable = True
        self.viewer.config.regions.changeable = True

        # Disable
        self.viewer.config.page = False
        self.viewer.config.show = False

        # Run the viewer
        self.viewer.run(paths=self.image_paths, regions=self.ellipses)

    # -----------------------------------------------------------------

    @lazyproperty
    def page(self):

        """
        This function ...
        :return:
        """

        # Initialize the page
        self.viewer._initialize_page()

        # Return the page
        return self.viewer.page

    # -----------------------------------------------------------------

    @property
    def display_ids(self):

        """
        Thisf unction ...
        :return:
        """

        return self.viewer.display_ids

    # -----------------------------------------------------------------

    @property
    def all_loader(self):

        """
        This function ...
        :return:
        """

        return self.viewer.all_loader

    # -----------------------------------------------------------------

    @property
    def table(self):

        """
        This function ...
        :return:
        """

        return self.viewer.table

    # -----------------------------------------------------------------

    @property
    def preloader(self):

        """
        Thsi function ...
        :return:
        """

        return self.viewer.preloader

    # -----------------------------------------------------------------

    def add_buttons(self):

        """
        This function ...
        :return:
        """

        # Theme button
        self.viewer._add_theme_button()

        # Infs and negatives
        self.viewer._add_infs_button()
        self.viewer._add_negatives_button()

    # -----------------------------------------------------------------

    def add_loader(self):

        """
        Thisn function ...
        :return:
        """

        # All images loader
        self.viewer._add_all_loader()

    # -----------------------------------------------------------------

    def add_preloader(self):

        """
        This function ...
        :return:
        """

        # Add preloader
        self.viewer._add_preloader()

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Set dictionary of ellipse per display
        ellipses = dict()
        for name in self.ellipses:
            display_id = self.display_ids[name]
            ellipses[display_id] = self.ellipses[name]

        # Add regions synchronization
        self.page += make_synchronize_regions_script(self.indicator_id, self.display_ids.values(), ellipses, start_factor=self.default_factor) + "\n"

        # Add buttons
        self.add_buttons()

        # Add the factor indicator
        self.page += center(self.indicator) + newline

        # Add loader
        self.add_loader()

        # Add the table
        self.page += self.table

        # Add preloader
        self.add_preloader()

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
        with browser.serve_local_host(): browser.open_path(self.ellipse_page_path)

# -----------------------------------------------------------------
