#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.html.significance Contains the SignificanceLevelsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import TruncationComponent
from ...html.component import stylesheet_url, page_style
from ....core.tools.html import HTMLPage, updated_footing
from ....core.tools import html
from ....core.tools import browser
from ....core.tools.utils import lazyproperty
from ....core.tools import filesystem as fs
from ....magic.core.mask import Mask
from ....magic.core.alpha import AlphaMask
from ....core.tools import sequences
from ....core.filter.filter import parse_filter
from ....core.basics.range import RealRange

# -----------------------------------------------------------------

significance_plots_name = "significance_plots"
ncolumns = 2
colour_map = "jet"

page_width = 600

# -----------------------------------------------------------------

class SignificanceLevelsPageGenerator(TruncationComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SignificanceLevelsPageGenerator, self).__init__(*args, **kwargs)

        # Plot paths for each filter
        self.image_plot_paths = dict()

        # The sliders
        self.sliders = dict()

        # The plots path
        self.plots_path = None

        # The page
        self.page = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make plots
        self.make_plots()

        # 3. Make the sliders
        self.make_sliders()

        # 4. Generate the page
        self.generate_page()

        # 5. Writing
        self.write()

        # 6. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SignificanceLevelsPageGenerator, self).setup(**kwargs)

        # Make directory to contain the plots
        self.plots_path = fs.join(self.truncation_html_path, "significance")
        if fs.is_directory(self.plots_path):
            if self.config.replot: fs.clear_directory(self.plots_path)
        else: fs.create_directory(self.plots_path)

        # Create directories for each filter
        for name in self.names:
            image_path = fs.join(self.plots_path, name)
            self.image_plot_paths[name] = image_path
            if not fs.is_directory(image_path): fs.create_directory(image_path)

        # Check
        if self.config.remove_dark_masks and self.config.dark_masks: raise ValueError("Cannot ask to create dark masks and to remove dark masks")

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        return self.dataset.filters

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        #return self.dataset.names
        # SORT BASED ON WAVELENGTH
        return sorted(self.dataset.names, key=lambda name: parse_filter(name).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Significance maps"

    # -----------------------------------------------------------------

    @property
    def image_width(self):

        """
        This fucntion ...
        :return:
        """

        return self.config.image_width

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        This function ...
        :return:
        """

        return self.config.image_height

    # -----------------------------------------------------------------

    def has_all_plots(self, name):

        """
        This function ...
        :return:
        """

        # Get plot path
        #plot_path = self.filter_plot_paths[parse_filter(name)]
        plot_path = self.image_plot_paths[name]

        #print(fs.files_in_path(plot_path, returns="name", extension="png", convert=float), self.config.sigma_levels)

        # Loop over the levels
        for level in self.config.sigma_levels:

            # Determine path
            path = fs.join(plot_path, str(level) + ".png")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def check_plots(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the plot path
        plot_path = self.image_plot_paths[name]

        #has_all = True

        # Loop over the levels
        for level in self.config.sigma_levels:

            # Determine path
            path = fs.join(plot_path, str(level) + ".png")

            # Determine mask path
            mask_path = fs.join(plot_path, str(level) + "_mask.png")

            # Determine dark mask path
            dark_mask_path = fs.join(plot_path, str(level) + "_mask_dark.png")

            # Remove dark masks
            if self.config.remove_dark_masks and fs.is_file(dark_mask_path): fs.remove_file(dark_mask_path)

            # Check
            #if fs.is_file(path) and fs.is_file(mask_path): pass
            #else: has_all = False
            if not fs.is_file(path) or not fs.is_file(mask_path) or (self.config.dark_masks and not fs.is_file(dark_mask_path)): return False

        # Return
        #return has_all
        return True

    # -----------------------------------------------------------------

    def create_mask_for_level(self, significance, level):

        """
        This function ...
        :param significance:
        :param level:
        :return:
        """

        # Create the mask
        mask = Mask(significance > level)

        # Only keep largest patch
        mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

        # Fill holes
        mask.fill_holes()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def create_fuzzy_mask_for_level(self, significance, level):

        """
        This function ...
        :param significance:
        :param level:
        :return:
        """

        # Construct value range
        lower_relative = 1. - self.config.fuzziness # example: 1. - 0.1
        upper_relative = 1. + self.config.fuzziness
        value_range = RealRange.around(level, lower_relative, upper_relative)

        # Create the mask
        mask = AlphaMask.between(significance, value_range)

        # Only keep largest patch
        mask = mask.largest(npixels=self.config.min_npixels, connectivity=self.config.connectivity)

        # Fill holes
        mask.fill_holes()

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @lazyproperty
    def normalization_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.normalization_ellipse_factor

    # -----------------------------------------------------------------

    def make_plots_for_level(self, name, frame, significance, plot_path, level):

        """
        This function ...
        :param name:
        :param frame:
        :param significance:
        :param plot_path:
        :param level:
        :return:
        """

        # Determine path
        path = fs.join(plot_path, str(level) + ".png")

        # Determine path for the mask
        mask_path = fs.join(plot_path, str(level) + "_mask.png")

        # Determine path for the dark mask
        mask_dark_path = fs.join(plot_path, str(level) + "_mask_dark.png")

        # Check
        if fs.is_file(path) and fs.is_file(mask_path) and (not self.config.dark_masks or fs.is_file(mask_dark_path)):
            log.success("The plots for the '" + name + "' image at a sigma level of '" + str(level) + "' are already present")
            return

        # Debugging
        log.debug("Making plots for the '" + name + "' image at a sigma level of '" + str(level) + "' ...")

        # Create the mask
        if self.config.fuzzy_mask: mask = self.create_fuzzy_mask_for_level(significance, level)
        else: mask = self.create_mask_for_level(significance, level)

        # Mask the frame
        if self.config.fuzzy_mask: frame = frame.applied_alpha_mask(mask)
        else: frame = frame.applied_mask(mask, invert=True)

        # Save the frame
        frame.saveto_png(path, colours=self.config.colours, interval=self.config.interval, scale=self.config.scale,
                         alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha,
                         normalize_in=self.normalization_ellipse)

        # Save the mask
        mask.saveto_png(mask_path, colour=self.config.mask_colour, alpha=self.config.mask_alpha)

        # Save the dark mask
        if self.config.dark_masks: mask.saveto_png(mask_dark_path, colour=self.config.dark_mask_colour, alpha=True) # background has to be transparent

    # -----------------------------------------------------------------

    @lazyproperty
    def levels(self):

        """
        This function ...
        :return:
        """

        # Levels have already been set
        if self.has_significance_levels:

            # Get configured levels
            levels = set(self.config.sigma_levels)

            # Loop over the image names
            for name in self.names:

                # Get the level
                level = self.significance_levels[name]

                # Add
                levels.add(level)

            # Return the levels, sorted
            return list(sorted(levels))

        # Levels are not yet defined
        else: return self.config.sigma_levels

    # -----------------------------------------------------------------

    def get_default_level(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_significance_levels: return self.significance_levels[name]
        else: return sequences.find_closest_value(self.levels, self.config.default_level)

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Loop over the frames
        for name in self.names:

            # Check whether not all plots are already present
            if self.check_plots(name): continue

            # Get plot path
            plot_path = self.image_plot_paths[name]

            # Get frame
            frame = self.dataset.get_frame(name)

            # Crop the frame
            frame.crop_to(self.truncation_box, factor=self.config.cropping_factor)

            # Get error map and crop as well
            errormap = self.dataset.get_errormap(name)
            errormap.rebin(frame.wcs, convert=True)

            # Create the significance map
            significance = frame / errormap

            # Create the plots
            for level in self.levels: self.make_plots_for_level(name, frame, significance, plot_path, level)

            # Clean
            gc.collect()

    # -----------------------------------------------------------------

    def make_sliders(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the image sliders ...")

        # Loop over the images
        for name in self.names:

            # Determine paths to the plots
            labels = self.config.sigma_levels
            paths = [fs.join(self.image_plot_paths[name], str(level) + ".png") for level in labels]
            mask_paths = [fs.join(self.image_plot_paths[name], str(level) + "_mask.png") for level in labels]

            # Create image ID
            image_id = name.replace("_", "").replace(" ", "")

            # Get the default level
            default_level = self.get_default_level(name)

            # Create the slider
            slider = html.make_image_slider(image_id, mask_paths, labels, default_level,
                                            width=self.image_width, height=self.image_height, basic=True,
                                            img_class="pixelated invertable", extra_urls=paths, extra_img_class="pixelated")

            # Set the slider
            self.sliders[name] = slider

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Create list of css scripts
        css_paths = [stylesheet_url]

        # Create CSS for the page width
        css = html.make_page_width(page_width)

        # Make javascripts urls
        javascript_paths = []

        # LATEST ATTEMPT
        javascript_paths_body = None

        # Create the page
        self.page = HTMLPage(self.title, style=page_style, css=css, css_path=css_paths, javascript_path=javascript_paths, javascript_path_body=javascript_paths_body, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes))

        self.page += html.newline

        self.page += html.line
        self.page += html.newline

        # Add the sliders
        for name in self.names:

            # Add the name
            self.page += name.upper()
            self.page += html.newline

            # Add the image slider
            self.page += self.sliders[name]

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

        # Write the page
        self.write_page()

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Save
        self.page.saveto(self.significance_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        with browser.serve_local_host(): browser.open_path(self.significance_page_path)

# -----------------------------------------------------------------
