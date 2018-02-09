#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.plot.imagegrid Contains the ImageGridPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from collections import OrderedDict
from abc import ABCMeta, abstractproperty
from collections import defaultdict
from matplotlib.colors import rgb2hex

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from ...core.tools import types
from ..basics.mask import MaskBase
from ..region.list import RegionList
from ..region.region import Region
from ..region.list import region_to_region_list
from ..core.image import Image
from ..core.frame import Frame, AllZeroError
from ...core.tools import filesystem as fs
from ...core.tools.numbers import nan
from ...core.basics.containers import ordered_by_value
from ...core.basics.plot import MPLFigure, BokehFigure, mpl, bokeh, dark_pretty_colors, pretty_colors, filled_markers
from ..tools import plotting
from ...core.tools.utils import lazyproperty
from ..region.list import load_region_list
from ..core.mask import Mask
from ...core.tools import numbers
from ...core.basics.map import Map
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.tools.stringify import tostr
from ...core.plot.distribution import plot_distribution
from ...core.tools import strings

# -----------------------------------------------------------------

class ImageGridPlotter(Configurable):

    """
    This class ...
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
        super(ImageGridPlotter, self).__init__(*args, **kwargs)

        # Set the title
        self.title = None

        # The figure
        self.figure = None

        # The plots
        self.plots = None

        # Output path
        self.out_path = None

        # Rebinning
        self.rebin_to = None

        # Cropping
        self.crop_to = None
        self.cropping_factor = 1.2

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImageGridPlotter, self).setup(**kwargs)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        if self.out_path is None and "output" in self.config and self.config.output is not None:
            full_output_path = fs.absolute_or_in(self.config.output, self.config.path)
            if fs.has_extension(full_output_path):
                directory_path = fs.directory_of(full_output_path)
                if not fs.is_directory(directory_path): fs.create_directory(directory_path)
            elif not fs.is_directory(full_output_path): fs.create_directory(full_output_path)
            self.out_path = full_output_path

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

        # Rebinning
        if "rebin_to" in kwargs: self.rebin_to = kwargs.pop("rebin_to")

        # Cropping
        if "crop_to" in kwargs: self.crop_to = kwargs.pop("crop_to")
        if "cropping_factor" in kwargs: self.cropping_factor = kwargs.pop("cropping_factor")

    # -----------------------------------------------------------------

    def initialize_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Initializing the figure with size " + str(self.figsize) + " ...")

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

    # -----------------------------------------------------------------

    @property
    def preserve_units(self):

        """
        This function ...
        :return:
        """

        return not self.config.normalize

    # -----------------------------------------------------------------

    def _process_frame(self, name, frame, copy=True, unit=None, return_normalization_sum=False, normalization_sum=None):

        """
        This function ...
        :param name:
        :param frame:
        :param copy:
        :param unit:
        :param return_normalization_sum:
        :param normalization_sum:
        :return:
        """

        # Initialize flags
        _rebinned = False
        _cropped = False
        _downsampled = False
        _normalized = False
        _converted = False

        # Initialize normalization sum
        norm_sum = None

        # Rebin the frame
        if self.rebin_to is not None:

            # Debugging
            log.debug("Rebinning the '" + name + "' frame ...")

            # Rebin
            if copy: frame = frame.rebinned(self.rebin_to, convert=self.preserve_units)
            else: frame.rebin(self.rebin_to, convert=self.preserve_units)

            # Set flag
            _rebinned = True

        # Crop the frame
        if self.crop_to is not None:

            # Debugging
            log.debug("Cropping the '" + name + "' frame ...")

            # Crop
            if _rebinned or not copy: frame.crop_to(self.crop_to, factor=self.cropping_factor, out_of_bounds="expand")
            else: frame = frame.cropped_to(self.crop_to, factor=self.cropping_factor, out_of_bounds="expand")

            # Set flag
            _cropped = True

        # Downsampling required?
        if self.config.downsample and frame.xsize > self.config.max_npixels or frame.ysize > self.config.max_npixels:

            # Determine the downsample factor
            downsample_factor = max(frame.xsize, frame.ysize) / float(self.config.max_npixels)

            # Debugging
            log.debug("Downsampling the '" + name + "' frame with a factor of " + str(downsample_factor) + "...")

            # Downsample
            if (_rebinned or _cropped) or not copy: frame.downsample(downsample_factor, convert=self.preserve_units, dilate_nans=False, dilate_infs=False)
            else: frame = frame.downsampled(downsample_factor, convert=self.preserve_units, dilate_nans=False, dilate_infs=False)

            # Set flag
            _downsampled = True

        # Normalize the frame
        if self.config.normalize:

            # Debugging
            log.debug("Normalizing the '" + name + "' frame ...")

            # Normalize this frame based on the sum of another frame
            if normalization_sum is not None:

                # Normalize
                if (_rebinned or _cropped or _downsampled) or not copy: frame /= normalization_sum
                else: frame = frame / normalization_sum

                # Remove the unit
                frame.unit = None

            # Try normalization this frame independently
            else:

                # Try
                try:

                    # Create normalized frame
                    if (_rebinned or _cropped or _downsampled) or not copy: norm_sum = frame.normalize()
                    else: frame, norm_sum = frame.normalized(return_sum=True)

                    # Set flag
                    _normalized = True

                except AllZeroError: log.warning("The '" + name + "' frame could not be normalized")

        # Convert the unit
        if unit is not None and frame.unit != unit:

            # Debugging
            log.debug("Converting the unit of the '" + name + "' frame to " + tostr(unit) + " ...")

            # Convert the unit
            if (_rebinned or _cropped or _downsampled or _normalized) or not copy: frame.convert_to(unit)
            else: frame = frame.converted_to(unit)

            # Set flag
            _converted = True

        # If no processing is performed, make a copy if necessary
        processed = _rebinned or _normalized or _cropped or _downsampled or _converted
        if not processed and copy: frame = frame.copy()

        # Return the new frame
        if return_normalization_sum: return frame, norm_sum
        else: return frame

    # -----------------------------------------------------------------

    @abstractproperty
    def nrows(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def ncolumns(self):

        """
        Thisn function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def width(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def height(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def figsize(self):

        """
        This function ...
        :return:
        """

        return (self.width, self.height)

    # -----------------------------------------------------------------

    def finish_plot(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Finishing the plot ...")

        # Tight layout
        self.figure.tight_layout()

        # Add title if requested
        if self.title is not None: self.figure.set_title(self.title)

        # Show the plot
        if self.config.show: self.figure.show()

        # Save the figure
        if self.out_path is not None: self.save_figure()

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Close the figure
        if self.figure is not None: self.figure.close()

    # -----------------------------------------------------------------

    def save_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Saving figure ...")

        # Determine path
        if types.is_string_type(self.out_path):
            if fs.is_directory(self.out_path): path = fs.join(self.out_path, "images." + self.config.format)
            else: path = self.out_path
        else: path = self.out_path

        # Save
        self.figure.saveto(path)

# -----------------------------------------------------------------

class StandardImageGridPlotter(ImageGridPlotter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(StandardImageGridPlotter, self).__init__(*args, **kwargs)

        # The frames to be plotted
        self.frames = OrderedDict()

        # Masks and regions to be overlayed on the images
        self.masks = defaultdict(dict)
        self.regions = defaultdict(dict)

        # The colorbar
        self.colorbar = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Uniformize the frames
        self.uniformize()

        # 2. Write
        if self.config.write: self.write()

        # 3. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    @lazyproperty
    def use_coordinates(self):

        """
        This function ...
        :return:
        """

        if self.config.coordinates is not None: return self.config.coordinates
        else: return self.all_have_wcs

    # -----------------------------------------------------------------

    @property
    def share_x(self):

        """
        This function ...
        :return:
        """

        #return False
        return self.use_coordinates

    # -----------------------------------------------------------------

    @property
    def share_y(self):

        """
        This function ...
        :return:
        """

        #return False
        return self.use_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def projections(self):

        """
        This function ...
        :return:
        """

        # Don't plot based on the coordinates
        if not self.use_coordinates: return None

        # Initialize list structure
        projections = sequences.create_nested_2d(self.nrows, self.ncolumns)

        # Loop over the images
        for index, name in enumerate(self.names):

            # Get row and col
            row, col = self.get_row_and_col(index)

            # Get the coordinate system
            wcs = self.get_wcs(name)

            # Loop over the column indices
            projections[row][col] = wcs.to_astropy()

        # Return the projections
        return projections

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(StandardImageGridPlotter, self).setup(**kwargs)

        # Load the images
        if self.no_frames:
            if self.config.from_data: self.load_data()
            else: self.load_images()

        # Still no images
        if self.no_frames: raise RuntimeError("No frames are loaded")

        # Initialize the figure
        self.initialize_figure()

        # Setup the plots
        #if self.use_coordinates: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white", projection=self.projection)
        #else: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white")

        # Create plots
        self.plots = self.figure.create_grid(self.nrows, self.ncolumns, sharex=self.share_x, sharey=self.share_y, projections=self.projections)

        # Sort the frames on filter
        if self.config.sort_filters: self.sort()

        # Show the grid
        if self.config.show_grid: self.show_grid()

    # -----------------------------------------------------------------

    def show_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the grid ...")

        # Create plots
        viz_figure = MPLFigure(size=self.figsize)
        viz_plots = viz_figure.create_grid(self.nrows, self.ncolumns, sharex=self.share_x, sharey=self.share_y)

        #for plots_row in viz_plots:
        #    for plot in plots_row:
        for plot in sequences.iterate_2d(viz_plots):

            axes = plot.axes
            # print(axes)

            adress = str(hex(id(axes)))

            print('<%s.%s object at %s>' % (
                axes.__class__.__module__,
                axes.__class__.__name__,
                hex(id(axes))
            ))

            print("x")
            for ax in axes.get_shared_x_axes(): print(ax)
            print("y")
            for ax in axes.get_shared_y_axes(): print(ax)
            # print(vars(axes))

            # Plot the text
            plot.text(0.5, 0.5, adress)

        # Show
        viz_figure.show()

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.frames[name].wcs

    # -----------------------------------------------------------------

    def has_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_wcs(name) is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def all_have_wcs(self):

        """
        This function ...
        :return:
        """

        for name in self.names:
            if not self.has_wcs(name): return False
        return True

    # -----------------------------------------------------------------

    @property
    def first_name(self):

        """
        This function ...
        :return:
        """

        return self.names[0]

    # -----------------------------------------------------------------

    @property
    def last_name(self):

        """
        This function ...
        :return:
        """

        return self.names[-1]

    # -----------------------------------------------------------------

    @property
    def first_wcs(self):

        """
        This function ...
        :return:
        """

        return self.get_wcs(self.first_name)

    # -----------------------------------------------------------------

    @property
    def projection(self):

        """
        This function ...
        :return:
        """

        return self.first_wcs

    # -----------------------------------------------------------------

    def sort(self):

        """
        This function ...
        :return:
        """

        # Sort on filter
        self.frames = ordered_by_value(self.frames, key=lambda frame: frame.wavelength)

    # -----------------------------------------------------------------

    @property
    def fixed_ncolumns(self):

        """
        This function ...
        :return:
        """

        return self.config.fixed == "columns"

    # -----------------------------------------------------------------

    @property
    def fixed_nrows(self):

        """
        This function ...
        :return:
        """

        return self.config.fixed == "rows"

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        # Determine the necessary number of rows, for a fixed number of columns
        if self.fixed_ncolumns: return int(math.ceil(float(self.nframes) / self.ncolumns))

        # Fixed number of rows
        else: return min(self.nframes, self.config.nrows)

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        Thisn function ...
        :return:
        """

        # Determine the necessary number of columns, for a fixed number of rows
        if self.fixed_nrows: return int(math.ceil(float(self.nframes) / self.nrows))

        # Fixed number of columns
        else: return min(self.nframes, self.config.ncolumns)

    # -----------------------------------------------------------------

    @property
    def ratio(self):

        """
        This function ...
        :return:
        """

        return float(self.nrows) / float(self.ncolumns)

    # -----------------------------------------------------------------

    @property
    def width(self):

        """
        This function ...
        :return:
        """

        return self.config.plot.xsize * self.ncolumns * self.mean_width_to_height

    # -----------------------------------------------------------------

    @property
    def height(self):

        """
        This function ...
        :return:
        """

        return self.config.plot.ysize * self.nrows

    # -----------------------------------------------------------------

    @property
    def mean_width_to_height(self):

        """
        This function ...
        :return:
        """

        values = []

        # Loop over the frames
        for name in self.names: values.append(self.get_width_to_height(name))

        # Return the mean
        return numbers.arithmetic_mean(*values)

    # -----------------------------------------------------------------

    def get_width_to_height(self, name):

        """
        This function ...
        :return:
        """

        width = self.get_width(name)
        height = self.get_height(name)
        return float(width) / float(height)

    # -----------------------------------------------------------------

    def get_width(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.frames[name].xsize

    # -----------------------------------------------------------------

    def get_height(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.frames[name].ysize

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading data from file ...")

        # Loop over the directories
        for path, name in fs.directories_in_path(self.config.path, returns=["path", "name"]):

            # Debugging
            log.debug("Loading the '" + name + "' frame ...")

            # Determine the frame path
            frame_path = fs.join(path, "frame.fits")
            if not fs.is_file(frame_path): raise IOError("Frame could not be found")

            # Load the frame
            frame = Frame.from_file(frame_path)

            # Determine the regions path
            regions_path = fs.join(path, "regions")
            if fs.is_directory(regions_path) and not fs.is_empty(regions_path):
                regions = dict()
                for regions_path, regions_name in fs.files_in_path(regions_path, extension="reg", returns=["path", "name"]):
                    region_list = load_region_list(regions_path)
                    regions[regions_name] = region_list
            else: regions = None

            # Determine the masks path
            masks_path = fs.join(path, "masks")
            if fs.is_directory(masks_path) and not fs.is_empty(masks_path):
                masks = dict()
                for mask_path, mask_name in fs.files_in_path(masks_path, extension="fits", returns=["path", "name"]):
                    mask = Mask.from_file(mask_path)
                    masks[mask_name] = mask
            else: masks = None

            # Add the frame, unprocessed
            self.add_frame(frame, name=name, mask=masks, regions=regions, process=False)

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading images from file ...")

        # Loop over the files in the directory
        for path, name in fs.files_in_path(self.config.path, extension="fits", contains=self.config.contains,
                                           not_contains=self.config.not_contains, exact_name=self.config.exact_name,
                                           exact_not_name=self.config.exact_not_name, startswith=self.config.startswith,
                                           endswith=self.config.endswith, returns=["path", "name"]):

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # Load as image?
            if self.config.planes: self.add_image_from_file(path, name, all_frames=True, masks=self.config.masks, regions=self.config.regions)

            # Load as frame
            else: self.add_frame_from_file(path, name)

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return len(self.frames)

    # -----------------------------------------------------------------

    @property
    def has_frames(self):

        """
        This function ...
        :return:
        """

        return self.nframes > 0

    # -----------------------------------------------------------------

    @property
    def no_frames(self):

        """
        This function ...
        :return:
        """

        return self.nframes == 0

    # -----------------------------------------------------------------

    @property
    def nmasks(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.masks)

    # -----------------------------------------------------------------

    @property
    def has_masks(self):

        """
        This function ...
        :return:
        """

        return self.nmasks > 0

    # -----------------------------------------------------------------

    @property
    def no_masks(self):

        """
        This function ...
        :return:
        """

        return self.nmasks == 0

    # -----------------------------------------------------------------

    @property
    def nregions(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.regions)

    # -----------------------------------------------------------------

    @property
    def has_regions(self):

        """
        This function ...
        :return:
        """

        return self.nregions > 0

    # -----------------------------------------------------------------

    @property
    def no_regions(self):

        """
        This function ...
        :return:
        """

        return self.nregions == 0

    # -----------------------------------------------------------------

    def nregions_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.regions: return 0
        else: return len(self.regions[name])

    # -----------------------------------------------------------------

    def has_regions_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nregions_for_frame(name) > 0

    # -----------------------------------------------------------------

    def no_regions_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nregions_for_frame(name) == 0

    # -----------------------------------------------------------------

    def nmasks_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.masks: return 0
        else: return len(self.masks[name])

    # -----------------------------------------------------------------

    def has_masks_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nmasks_for_frame(name) > 0

    # -----------------------------------------------------------------

    def no_masks_for_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nmasks_for_frame(name) == 0

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        Thisn function ...
        :return:
        """

        return self.frames.keys()

    # -----------------------------------------------------------------

    def has_frame(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return name in self.names

    # -----------------------------------------------------------------

    def add_image(self, image, name=None, all_frames=False, masks=True, regions=True, mask_name=None, regions_name=None,
                  copy=True, process=True):

        """
        This function ...
        :param image:
        :param name:
        :param all_frames:
        :param masks:
        :param regions:
        :param mask_name:
        :param regions_name:
        :param copy:
        :param process:
        :return:
        """

        # NO: don't do the downsampling here, why would we.. it is done for each frame (and mask) separately in add_frame,
        # by doing it here we are downsampling other frames and mask we don't actually need imported
        # # Downsampling required?
        # copied = False
        # if process and self.config.downsample and image.xsize > self.config.max_npixels or image.ysize > self.config.max_npixels:
        #
        #     # Determine the downsample factor
        #     downsample_factor = max(image.xsize, image.ysize) / float(self.config.max_npixels)
        #
        #     # Downsample
        #     if copy:
        #         image = image.copy()
        #         copied = True
        #     image.downsample(downsample_factor, convert=False, dilate_nans=False, dilate_infs=False)

        # Is copy still required?
        #copy = copy and not copied

        # Remember the names of the frames that we have added
        frame_names = []

        # Use all frames
        if all_frames:

            # Check
            if name is not None: raise ValueError("When using all frames, name cannot be specified")

            # Loop over the frames
            for frame_name in image.frame_names:

                # Get the frame
                frame = image.frames[frame_name]

                # Add the frame
                self.add_frame(frame, name=frame_name, copy=copy, process=process)

                # Add the frame name
                frame_names.append(frame_name)

        # Use only the primary frame
        else:

            # Set name
            if name is None: name = image.name
            if name is None: raise ValueError("Name of the image is not defined")

            # Get and add the primary frame
            frame = image.primary
            self.add_frame(frame, name, copy=copy, process=process)

            # Add the frame name
            frame_names = [name]

        # Add one mask
        if mask_name is not None:

            # Get the mask
            mask = image.masks[mask_name]

            # Loop over the frames
            for name in frame_names: self.add_mask(name, mask, label=mask_name)

        # Add all masks
        elif masks:

            # Loop over the frames
            for name in frame_names: self.add_masks(name, image.masks)

        # Add one region list
        if regions_name is not None:

            # Get the region list
            regions = image.regions[regions_name]

            # Loop over the frames
            for name in frame_names: self.add_region_list(name, regions, label=regions_name)

        # Add all region lists
        elif regions:

            # Loop over the frames
            for name in frame_names: self.add_region_lists(name, image.regions)

    # -----------------------------------------------------------------

    def add_image_from_file(self, path, name=None, all_frames=False, masks=True, regions=True, mask_name=None,
                            regions_name=None, copy=True):

        """
        This function ...
        :param path:
        :param name:
        :param all_frames:
        :param masks:
        :param regions:
        :param mask_name:
        :param regions_name:
        :param copy:
        :return:
        """

        # Load the image
        image = Image.from_file(path)

        # Add
        self.add_image(image, name=name, all_frames=all_frames, masks=masks, regions=regions, mask_name=mask_name,
                       regions_name=regions_name, copy=copy)

    # -----------------------------------------------------------------

    @property
    def preserve_units(self):

        """
        This function ...
        :return:
        """

        return not self.config.normalize

    # -----------------------------------------------------------------

    def add_frame(self, frame, name=None, mask=None, regions=None, copy=True, process=True):

        """
        This function ...
        :param frame:
        :param name:
        :param mask:
        :param regions:
        :param copy:
        :param process:
        :return:
        """

        # Check name
        if name is None: name = frame.name
        if name is None: raise ValueError("Frame has no name")

        # Check name
        if self.has_frame(name): raise ValueError("Name '" + name + "' is already in use")

        # Process the frame
        if process: frame = self._process_frame(name, frame, copy=copy)
        elif copy: frame = frame.copy()

        # Add the frame
        self.frames[name] = frame

        # Add masks
        if types.is_dictionary(mask): self.add_masks(name, mask)
        elif isinstance(mask, MaskBase): self.add_mask(name, mask)
        elif mask is None: pass
        else: raise ValueError("Invalid input")

        # Add regions
        if types.is_dictionary(regions): self.add_region_lists(name, regions)
        elif isinstance(regions, RegionList): self.add_region_list(name, regions)
        elif isinstance(regions, Region): self.add_region(name, regions)
        elif regions is None: pass
        else: raise ValueError("Invalid input")

    # -----------------------------------------------------------------

    def add_frame_from_file(self, path, name=None, process=True):

        """
        This function ...
        :param path:
        :param name:
        :param process:
        :return:
        """

        # Get the frame
        frame = Frame.from_file(path)

        # Add the frame
        self.add_frame(frame, name=name, process=process)

    # -----------------------------------------------------------------

    def has_mask(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        return name in self.masks and label in self.masks[name]

    # -----------------------------------------------------------------

    def add_mask(self, name, mask, label=None):

        """
        This function ...
        :param name:
        :param mask:
        :param label:
        :return:
        """

        # Set mask label
        if label is None: label = name

        # Check label
        if self.has_mask(name, label): raise ValueError("Label '" + label + "' for frame '" + name + "' is already in use")

        # Add
        self.masks[name][label] = mask

    # -----------------------------------------------------------------

    def add_masks(self, name, masks):

        """
        This function ...
        :param name:
        :param masks:
        :return:
        """

        # Loop over the masks
        for label in masks: self.add_mask(name, masks[label], label=label)

    # -----------------------------------------------------------------

    def has_region(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        return name in self.regions and label in self.regions[name]

    # -----------------------------------------------------------------

    def add_region(self, name, region, label=None):

        """
        This function ...
        :param name:
        :param region:
        :param label:
        :return:
        """

        # Create region list
        regions = region_to_region_list(region)

        # Add the region list
        self.add_region_list(name, regions, label=label)

    # -----------------------------------------------------------------

    def add_region_list(self, name, regions, label=None):

        """
        This function ...
        :param name:
        :param regions:
        :param label:
        :return:
        """

        # Set label
        if label is None: label = name

        # Check label
        if self.has_region(name, label): raise ValueError("Label '" + label + "' for frame '" + name + "' is already in use")

        # Add
        self.regions[name][label] = regions

    # -----------------------------------------------------------------

    def add_region_lists(self, name, regions):

        """
        This function ...
        :param name:
        :param regions:
        :return:
        """

        # Loop over the regions
        for label in regions:

            # Get the region(-list)
            region = regions[label]

            # Region or region list?
            if isinstance(region, RegionList): self.add_region_list(name, region, label=label)
            elif isinstance(region, Region): self.add_region(name, region, label=label)
            else: raise ValueError("Invalid input")

    # -----------------------------------------------------------------

    @lazyproperty
    def colormap(self):

        """
        This function ...
        :return:
        """

        return cm.get_cmap(self.config.colormap)

    # -----------------------------------------------------------------

    @lazyproperty
    def background_color(self):

        """
        This function ...
        :return:
        """

        return self.colormap(0.0)

    # -----------------------------------------------------------------

    def index_of_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.names.index(name)

    # -----------------------------------------------------------------

    def _plot_frame(self, name, row, col, vmin=None, vmax=None, return_image=False, return_normalization=False, add_label=True):

        """
        This function ...
        :param name:
        :param row:
        :param col:
        :param vmin:
        :param vmax:
        :param return_image:
        :param return_normalization:
        :param add_label:
        :return:
        """

        # Debugging
        log.debug("Adding the '" + name + "' frame to the plot ...")

        # Get the frame
        frame = self.frames[name]

        # Get the plot
        plot = self.plots[row][col]

        # Color spines
        #print(plot.axes.spines.keys())
        plot.axes.spines['bottom'].set_color("white")
        plot.axes.spines['top'].set_color("white")
        plot.axes.spines['left'].set_color("white")
        plot.axes.spines['right'].set_color("white")

        # Color ticks
        #plot.axes.xaxis.label.set_color("white")
        #plot.axes.yaxis.label.set_color("white")
        plot.axes.tick_params(axis='x', colors="white", direction="inout")
        plot.axes.tick_params(axis='y', colors="white", direction="inout")

        # Set background color: otherwise NaNs are not plotted (-> white/transparent)
        if self.config.background: plot.axes.set_facecolor(self.background_color)

        # Add mask if present
        if self.has_masks_for_frame(name):
            for label in self.masks[name]:
                mask = self.masks[name][label]
                frame[mask] = nan

        # Plot
        if self.config.share_scale and vmin is not None: interval = [vmin, vmax]
        else: interval = self.config.interval

        #from matplotlib.axes import subplot_class_factory
        # Set projection
        #projection_class, extra_kwargs = frame.wcs._as_mpl_axes()
        #projection_class, kwargs, key = process_projection_requirements(self.figure.figure, *args, **kwargs)
        #a = subplot_class_factory(projection_class)(self.figure.figure, **extra_kwargs)
        #plot.sca(a)

        # Plot
        plot.axes.set_adjustable('box-forced')
        vmin_image, vmax_image, image, normalization = plotting.plot_box(frame.data, axes=plot.axes, interval=interval, scale=self.config.scale, cmap=self.colormap, alpha=self.config.alpha, return_image=True, return_normalization=True)
        #plot.axes.set_adjustable('box-forced')
        #ax2.set_adjustable('box-forced')

        # Add region if present
        if self.has_regions_for_frame(name):
            for label in self.regions[name]:
                for patch in self.regions[name][label].to_mpl_patches(): plot.axes.add_patch(patch)

        # Add the label
        if add_label: plot.axes.text(0.95, 0.95, name, color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")  # fontweight='bold'

        # Return vmin and vmax
        if return_image:
            if return_normalization: return vmin_image, vmax_image, image, normalization
            else: return vmin_image, vmax_image, image
        else:
            if return_normalization: return vmin_image, vmax_image, normalization
            else: return vmin_image, vmax_image

    # -----------------------------------------------------------------

    def _plot_empty(self, row, col):

        """
        This function ...
        :param row:
        :param col:
        :return:
        """

        # Get the plot
        plot = self.plots[row][col]

        # Color spines
        plot.axes.spines['bottom'].set_color("white")
        plot.axes.spines['top'].set_color("white")
        plot.axes.spines['left'].set_color("white")
        plot.axes.spines['right'].set_color("white")

        # Color ticks
        #plot.axes.xaxis.label.set_color("white")
        #plot.axes.yaxis.label.set_color("white")
        plot.axes.tick_params(axis='x', colors="white", direction="inout")
        plot.axes.tick_params(axis='y', colors="white", direction="inout")

        # Set background color: otherwise NaNs are not plotted (-> white/transparent)
        if self.config.background: plot.axes.set_facecolor(self.background_color)

        # Create NaNs image
        #nans = np.full((100,100), np.nan)

        #vmin_image, vmax_image = plotting.plot_box(nans, axes=plot.axes, interval=interval, scale=self.config.scale, cmap=self.colormap, alpha=self.config.alpha)

        #extent = None
        #aspect = "equal"
        #plot.axes.imshow(nans, origin="lower", vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, aspect=aspect, extent=extent)

    # -----------------------------------------------------------------

    def uniformize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uniformizing the frames ...")


    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the frames
        self.write_frames()

        # Write the regions
        if self.has_regions: self.write_regions()

        # Write the masks
        if self.has_masks: self.write_masks()

    # -----------------------------------------------------------------

    @lazyproperty
    def data_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("data", create=True)

    # -----------------------------------------------------------------

    def path_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.create_directory_in(self.data_path, name)

    # -----------------------------------------------------------------

    def frame_filepath_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.path_for_name(name)
        return fs.join(path, "frame.fits")

    # -----------------------------------------------------------------

    def write_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the frames ...")
        
        # Loop over the frames
        for name in self.names:

            # Get the frame
            frame = self.frames[name]

            # Determine path
            path = self.frame_filepath_for_name(name)

            # Save the frame
            frame.saveto(path)

    # -----------------------------------------------------------------

    def regions_path_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.path_for_name(name)
        return fs.create_directory_in(path, "regions")

    # -----------------------------------------------------------------

    def regions_filepath_for_name(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        path = self.regions_path_for_name(name)
        return fs.join(path, label + ".reg")

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the regions ...")

        # Loop over the names
        for name in self.names:

            # Loop over the regions
            for label in self.regions[name]:

                # Determine path
                path = self.regions_filepath_for_name(name, label)

                # Save the regions
                self.regions[name][label].saveto(path)

    # -----------------------------------------------------------------

    def masks_path_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.path_for_name(name)
        return fs.create_directory_in(path, "masks")

    # -----------------------------------------------------------------

    def mask_filepath_for_name(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        path = self.masks_path_for_name(name)
        return fs.join(path, label + ".")

    # -----------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks ...")

        # Loop over the names
        for name in self.names:

            # Loop over the masks
            for label in self.masks[name]:

                # Determine path
                path = self.mask_filepath_for_name(name, label)

                # Save the mask
                self.masks[name][label].saveto(path)

    # -----------------------------------------------------------------

    @property
    def npanels(self):

        """
        This function ...
        :return:
        """

        return self.nrows * self.ncolumns

    # -----------------------------------------------------------------

    def _plot_reference_frame(self):

        """
        This function ...
        :return:
        """

        # Get the name
        name = self.config.scale_reference

        # Get the index of the frame
        index = self.index_of_frame(self.config.scale_reference)

        # Get row and column index
        row = int(index / self.ncolumns)
        col = index % self.ncolumns

        # Plot the frame
        vmin, vmax = self._plot_frame(name, row, col)

        # Return
        return vmin, vmax, name

    # -----------------------------------------------------------------

    def get_row_and_col(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get row and column index
        row = int(index / self.ncolumns)
        col = index % self.ncolumns

        # Return
        return row, col

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Initialize vmin and vmax
        vmin = vmax = None

        # Initialize list
        plotted = []

        image = None
        normalization = None

        # First plot the image of which we use the scale as reference
        if self.config.share_scale and self.config.scale_reference is not None:

            # Plot
            vmin, vmax, name = self._plot_reference_frame()

            # Add to plotted
            plotted.append(name)

        # Loop over the images
        for index in range(self.npanels):

            # Get row and col
            row, col = self.get_row_and_col(index)

            # Plot a frame
            if index < self.nframes:

                # Get name
                name = self.names[index]
                if name in plotted: continue

                # Does this frame have any finite data?
                has_data = not (self.frames[name].all_zeroes or self.frames[name].all_nans or self.frames[name].all_infs)
                if has_data:

                    # Plot the frame
                    vmin_image, vmax_image, image, normalization = self._plot_frame(name, row, col, vmin=vmin, vmax=vmax, return_image=True, return_normalization=True)

                    # Set vmin and vmax
                    if self.config.share_scale:
                        vmin = vmin_image
                        vmax = vmax_image

                # No data
                else: self._plot_empty(row, col)

            # Empty
            else: self._plot_empty(row, col)

        # Check if any is plotted
        if image is None: raise RuntimeError("No image is plotted")

        # Set colorbar
        if self.colorbar is not None:

            self.figure.figure.colorbar(image, cax=self.colorbar)
            #print(normalization)
            #colorbar = self.colorbar.colorbar(image, norm=normalization) # doesn't set the scale

            #self.colorbar.solids.set_rasterized(True)
            #self.colorbar.outline.set_visible(False)
            #self.colorbar.set_ticks([])
            self.colorbar.get_xaxis().set_ticks([])
            self.colorbar.get_yaxis().set_ticks([])

            self.colorbar.set_axis_off()

            #print(vars(self.colorbar))

            # DOESN'T WORK
            #self.colorbar.yaxis._visible = False

            # DOESN'T WORK
            #from matplotlib.spines import Spine
            #for child in self.colorbar.get_children():
            #    if isinstance(child, Spine):
            #        child.set_color("white")

            # DOESN'T WORK
            # Color spines
            #self.colorbar.spines['bottom'].set_color("white")
            #self.colorbar.spines['top'].set_color("white")
            #self.colorbar.spines['left'].set_color("white")
            #self.colorbar.spines['right'].set_color("white")

            # Color ticks
            #self.colorbar.xaxis.label.set_color("white")
            #self.colorbar.yaxis.label.set_color("white")
            #self.colorbar.tick_params(axis='x', colors="white")
            #self.colorbar.tick_params(axis='y', colors="white")

        # Finish the plot
        self.finish_plot()

# -----------------------------------------------------------------

observation_index = 0
model_index = 1
residuals_index = 2
distribution_index = 3

# -----------------------------------------------------------------

class ResidualImageGridPlotter(ImageGridPlotter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ResidualImageGridPlotter, self).__init__(*args, **kwargs)

        # The rows of the grid
        self.rows = OrderedDict()

        # Masks and regions to be overlayed on the images
        self.masks = defaultdict(dict)
        self.regions = defaultdict(dict)

        # Error maps
        self.errors = dict()

        # The residual maps and distributions
        self.residuals = OrderedDict()
        self.distributions = OrderedDict()

        # The colorbars
        self.colorbars = None

        # Initialize list to contain the names of rows that have been plotted
        self._plotted_rows = []

        # Initialize vmin and vmax for the images and for the residuals
        self._vmin = None
        self._vmax = None
        self._vmin_res = None
        self._vmax_res = None

        # Images
        self._observation_images = OrderedDict()
        self._model_images = OrderedDict()
        self._residual_images = OrderedDict()
        self._distribution_images = OrderedDict()

        # Normalizations
        self._observation_normalizations = OrderedDict()
        self._model_normalizations = OrderedDict()
        self._residual_normalizations = OrderedDict()

    # -----------------------------------------------------------------

    def has_row(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return name in self.names

    # -----------------------------------------------------------------

    def _process_frames(self, name, observation, model, copy=True):

        """
        This function ...
        :param name:
        :param observation:
        :param model:
        :param copy:
        :return:
        """

        # Process observation
        if observation is not None:
            observation, normalization_sum = self._process_frame(name + " observation image", observation, copy=copy, return_normalization_sum=True)
            observation_unit = observation.unit
        else: observation_unit = normalization_sum = None

        # Process
        if model is not None: model = self._process_frame(name + " model image", model, copy=copy, unit=observation_unit, normalization_sum=normalization_sum)

        # Return
        return observation, model

    # -----------------------------------------------------------------

    def _check_masks(self, mask, name):

        """
        This function ...
        :param mask:
        :param name:
        :return:
        """

        if types.is_dictionary(mask): masks = mask
        elif isinstance(mask, MaskBase): masks = {name: mask}
        #elif mask is None: masks = None
        else: raise ValueError("Invalid input")
        return masks

    # -----------------------------------------------------------------

    def _check_regions(self, regions, name):

        """
        This function ...
        :param regions:
        :param name:
        :return:
        """

        if types.is_dictionary(regions): pass
        elif isinstance(regions, RegionList): regions = {name: regions}
        elif isinstance(regions, Region): regions = {name: region_to_region_list(regions)}
        #elif regions is None: pass
        else: raise ValueError("Invalid input")
        return regions

    # -----------------------------------------------------------------

    def add_row(self, observation, model, name, with_residuals=None, errors=None, masks=None, regions=None, copy=True,
                process=True, residuals=None):

        """
        This function ...
        :param observation:
        :param model:
        :param name:
        :param with_residuals: True or False (or None: automatic)
        :param errors:
        :param masks:
        :param regions:
        :param copy:
        :param process:
        :param residuals:
        :return:
        """

        # Maximum number of rows?
        if self.nrows_per_grid == self.config.max_nrows:
            log.warning("Cannot add the row '" + name + "': maximum number of rows has been reached")
            return False

        # Check name
        if self.has_row(name): raise ValueError("Name '" + name + "' is already in use")

        # Process the frames if necessary
        if process: observation, model = self._process_frames(name, observation, model, copy=copy)
        elif copy: observation, model = observation.copy(), model.copy()

        # Check masks
        if masks is not None: masks = self._check_masks(masks, name)

        # Check regions
        if regions is not None: regions = self._check_regions(regions, name)

        # Set residuals flag
        if with_residuals is None:

            # Residuals are passed
            if residuals is not None: with_residuals = True
            # If either observation or model is None, don't create residuals
            elif observation is None or model is None: with_residuals = False
            else: with_residuals = True

        # Make entry
        entry = Map()
        entry.observation = observation
        entry.model = model
        entry.residuals = with_residuals
        entry.errors = errors
        entry.masks = masks
        entry.regions = regions

        # Add the row
        self.rows[name] = entry

        # Add residuals
        if residuals is not None: self.add_residuals(residuals, name, existing_row=True)

        # If weighed residuals have to be plotted, we need error map
        elif self.config.weighed and errors is None: raise ValueError("Errors have to be specified to create weighed residuals")

        # Succesfully added the row
        return True

    # -----------------------------------------------------------------

    def add_observation(self, observation, name, replace=False, existing_row=False, copy=True, process=True, with_residuals=True):

        """
        This function ...
        :param observation:
        :param name:
        :param replace:
        :param existing_row:
        :param copy:
        :param process:
        :param with_residuals: assume residuals are going to be added or should be calculated
        :return:
        """

        # Already a row added
        if name in self.names:

            # Already an observation for this row
            if self.rows[name].observation is not None:

                # Replace?
                if replace: self._set_observation(observation, name, copy=copy, process=process)
                else: raise ValueError("Already an observation image for the row '" + name + "'")

            # Observation not yet defined for this row
            else: self._set_observation(observation, name, copy=copy, process=process)

            # Has been added succesfully
            return True

        # New row
        else:
            if existing_row: raise ValueError("There is no '" + name + "' row yet")
            return self.add_row(observation, None, name, copy=copy, process=process, with_residuals=with_residuals)

    # -----------------------------------------------------------------

    def _set_observation(self, observation, name, copy=True, process=True):

        """
        This function ...
        :param observation:
        :param name:
        :param copy:
        :param process:
        :return:
        """

        # Process the observation frame
        if process: observation = self._process_frame(name + " observation image", observation, copy=copy)

        # Set the observation frame
        self.rows[name].observation = observation

    # -----------------------------------------------------------------

    def add_observation_from_file(self, path, name=None, plane_index=None, plane=None, replace=False, existing_row=False,
                                  process=True, only_if_row=False):

        """
        This function ....
        :param path:
        :param name:
        :param plane_index:
        :param plane:
        :param replace:
        :param existing_row:
        :param process:
        :param only_if_row:
        :return:
        """

        # Determine name
        if name is None: name = fs.strip_extension(fs.name(path))
        if only_if_row and name not in self.names: return False

        # Load the frame
        observation = Frame.from_file(path, index=plane_index, plane=plane)

        # Add
        return self.add_observation(observation, name, replace=replace, existing_row=existing_row, copy=False, process=process)

    # -----------------------------------------------------------------

    def add_model(self, model, name, replace=False, existing_row=False, copy=True, process=True, with_residuals=True):

        """
        This function ...
        :param model:
        :param name:
        :param replace:
        :param existing_row:
        :param copy:
        :param process:
        :param with_residuals: assume residuals are going to be added or have to be calculated
        :return:
        """

        # Already a row added
        if name in self.names:

            # Already a model for this row
            if self.rows[name].model is not None:

                # Replace?
                if replace: self._set_model(model, name, copy=copy, process=process)
                else: raise ValueError("Already a model image for the '" + name + "' row")

            # Model not yet defined for this row
            else: self._set_model(model, name, copy=copy, process=process)

            # Has been added succesfully
            return True

        # New row
        else:
            if existing_row: raise ValueError("There is no '" + name + "' row yet")
            return self.add_row(None, model, name, copy=copy, process=process, with_residuals=with_residuals)

    # -----------------------------------------------------------------

    def _set_model(self, model, name, copy=True, process=True):

        """
        This function ...
        :param model:
        :param name:
        :param copy:
        :param process:
        :return:
        """

        # Get the observation unit
        observation_unit = self.get_observation_unit(name)

        # Process the model frame
        if process: model = self._process_frame(name + " model image", model, copy=copy, unit=observation_unit)

        # Set the model frame
        self.rows[name].model = model

    # -----------------------------------------------------------------

    def add_model_from_file(self, path, name=None, plane_index=None, plane=None, replace=False, existing_row=False,
                            process=True, only_if_row=False):

        """
        This function ...
        :param path:
        :param name:
        :param plane_index:
        :param plane:
        :param replace:
        :param existing_row:
        :param process:
        :param only_if_row:
        :return:
        """

        # Determine name
        if name is None: name = fs.strip_extension(fs.name(path))
        if only_if_row and name not in self.names: return False

        # Load the frame
        model = Frame.from_file(path, index=plane_index, plane=plane)

        # Add
        return self.add_model(model, name, replace=replace, existing_row=existing_row, copy=False, process=process)

    # -----------------------------------------------------------------

    def add_residuals(self, residuals, name, replace=False, existing_row=False, copy=True, process=True):

        """
        This function ...
        :param residuals:
        :param name:
        :param replace:
        :param existing_row:
        :param copy:
        :param process:
        :return:
        """

        # Existing row
        if name in self.names:

            # Already a residual map for this row
            if name in self.residuals:

                # Replace
                if replace: self._set_residuals(residuals, name, copy=copy, process=process)
                else: raise ValueError("Already a residual map for the '" + name + "' row")

            # Residual map not yet defined for this row
            else: self._set_residuals(residuals, name, copy=copy, process=process)

            # Has been added succesfully
            return True

        # New row
        else:
            if existing_row: raise ValueError("There is no '" + name + "' row yet")
            self._set_residuals(residuals, name, copy=copy, process=process)
            return True

    # -----------------------------------------------------------------

    def _process_residuals(self, name, residuals, copy=False):

        """
        This function ...
        :param name:
        :param residuals:
        :param copy:
        :return:
        """

        # Initialize flags
        _rebinned = False
        _converted = False

        # Get the row coordinate system
        row_wcs = self.get_wcs(name)

        # Rebin the frame
        if residuals.wcs != row_wcs is not None:

            # Debugging
            log.debug("Rebinning the '" + name + "' residuals frame ...")

            # Rebin
            if copy: residuals = residuals.rebinned(row_wcs)
            else: residuals.rebin(row_wcs)

            # Set flag
            _rebinned = True

        # Get the row unit
        row_unit = self.get_unit(name)

        # Convert the unit
        if residuals.unit != row_unit:

            # Debugging
            log.debug("Converting the unit of the '" + name + "' residuals frame to " + tostr(row_unit) + " ...")

            # Convert the unit
            if _rebinned or not copy: residuals.convert_to(row_unit)
            else: residuals = residuals.converted_to(row_unit)

            # Set flag
            _converted = True

        # If no processing is performed, make a copy if necessary
        processed = _rebinned or _converted
        if not processed and copy: residuals = residuals.copy()

        # Return the new residuals frame
        return residuals

    # -----------------------------------------------------------------

    def _set_residuals(self, residuals, name, copy=True, process=True):

        """
        This function ...
        :param residuals:
        :param name:
        :return:
        """

        # Process the residuals map
        if process: residuals = self._process_residuals(name, residuals, copy=copy)

        # Set the residuals frame
        self.residuals[name] = residuals

    # -----------------------------------------------------------------

    def add_residuals_from_file(self, path, name=None, replace=False, existing_row=False, process=True, only_if_row=False):

        """
        This function ....
        :param path:
        :param name:
        :param replace:
        :param existing_row:
        :param process:
        :param only_if_row:
        :return:
        """

        # Determine the name
        if name is None: name = fs.strip_extension(fs.name(path))
        if only_if_row and name not in self.names: return False

        # Load
        residuals = Frame.from_file(path)

        # Add
        return self.add_residuals(residuals, name, replace=replace, existing_row=existing_row, copy=False, process=process)

    # -----------------------------------------------------------------

    def add_distribution(self, distribution, name, replace=False, existing_row=False):

        """
        This function ...
        :param distribution:
        :param name:
        :param replace:
        :param existing_row:
        :return:
        """

        # Existing row
        if name in self.names:

            # Already a residual distribution for this row
            if name in self.distributions:

                # Replace
                if replace: self._set_distribution(distribution, name)
                else: raise ValueError("Already a distribution for the '" + name + "' row")

            # Distribution not yet defined for this row
            else: self._set_distribution(distribution, name)

            # Has been added succesfully
            return True

        # New row
        else:
            if existing_row: raise ValueError("There is no '" + name + "' row yet")
            self._set_distribution(distribution, name)
            return True

    # -----------------------------------------------------------------

    def _set_distribution(self, distribution, name):

        """
        This function ...
        :param distribution:
        :param name:
        :return:
        """

        self.distributions[name] = distribution

    # -----------------------------------------------------------------

    def add_distribution_from_file(self, path, name=None, replace=False, existing_row=False, only_if_row=False):

        """
        This function ...
        :param path:
        :param name:
        :param replace:
        :param existing_row:
        :param only_if_row:
        :return:
        """

        # Determine the name
        if name is None: name = fs.strip_extension(fs.name(path))
        if only_if_row and name not in self.names: return False

        # Load
        distribution = Distribution.from_file(path)

        # Add
        return self.add_distribution(distribution, name, replace=replace, existing_row=existing_row)

    # -----------------------------------------------------------------

    @property
    def has_all_residuals(self):

        """
        This function ...
        :return:
        """

        return sequences.same_contents(self.with_residuals_row_names, self.residuals_names)

    # -----------------------------------------------------------------

    @property
    def needs_residuals(self):

        """
        This function ...
        :return:
        """

        return not self.has_all_residuals

    # -----------------------------------------------------------------

    @property
    def has_all_distributions(self):

        """
        This function ...
        :return:
        """

        return sequences.same_contents(self.with_residuals_row_names, self.distribution_names)

    # -----------------------------------------------------------------

    @property
    def needs_distributions(self):

        """
        This function ...
        :return:
        """

        return self.config.distributions and not self.has_all_distributions

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Uniformize the frames
        # if self.needs_uniform: self.uniformize()

        # 2. Create the residuals
        if self.needs_residuals: self.create_residuals()

        # 3. Create the distributions
        if self.needs_distributions: self.create_distributions()

        # 4. Write
        if self.config.write: self.write()

        # 5. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    @property
    def single_grid(self):

        """
        This function ...
        :return:
        """

        return self.config.ngrids == 1

    # -----------------------------------------------------------------

    @property
    def add_residuals_colorbars(self):

        """
        This function ...
        :return:
        """

        return not self.config.distributions # if distributions are plotted (with the colors), colormap is not necessary

    # -----------------------------------------------------------------

    @property
    def share_residuals_colorbars(self):

        """
        This function ...
        :return:
        """

        return self.config.share_scale_residuals

    # -----------------------------------------------------------------

    @lazyproperty
    def use_coordinates(self):

        """
        This function ...
        :return:
        """

        if self.config.coordinates is not None: return self.config.coordinates
        else: return self.all_have_wcs

    # -----------------------------------------------------------------

    @lazyproperty
    def projections(self):

        """
        This function ...
        :return:
        """

        # Don't plot based on the coordinates
        if not self.use_coordinates: return None

        # Initialize list structure
        projections = sequences.create_nested_3d(self.config.ngrids, self.max_nrows_per_grid, self.ncolumns)

        # Loop over the images
        for index, name in enumerate(self.names):

            # Get the coordinate system
            wcs = self.get_wcs(name)

            # Get grid index and relative row index
            grid, rel_row_index = self.get_grid_and_relative_row(index)

            # Loop over the column indices
            for col in self.column_indices_not_distributions: projections[grid][rel_row_index][col] = wcs.to_astropy()

        # Return the projections
        return projections

    # -----------------------------------------------------------------

    @property
    def share_x(self):

        """
        This function ...
        :return:
        """

        #return False
        return self.use_coordinates

    # -----------------------------------------------------------------

    @property
    def share_y(self):

        """
        This function ...
        :return:
        """

        #return False
        return self.use_coordinates

    # -----------------------------------------------------------------

    @property
    def share_per_row(self):

        """
        This function ...
        :return:
        """

        #return True
        return False

    # -----------------------------------------------------------------

    @property
    def share_per_column(self):

        """
        This function ...
        :return:
        """

        #return True
        return False

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ResidualImageGridPlotter, self).setup(**kwargs)

        # Check config
        #if self.config.weighed and not self.config.relative: raise ValueError("Cannot enable 'weighed' and disable 'relative' residual maps")
        if not self.config.relative and self.config.normalize: raise ValueError("Cannot plot the actual (not relative) values of the residuals when normalizing the frames is enabled")

        # Other check
        #if self.config.share_scale_residuals and (self.config.same_residuals_scale and not self.config.share_scale): raise ValueError("Cannot share scales between residual maps and images, while also sharing the scale between the different residual maps, while the scales of the different images are not shared")

        # Load the images
        if self.no_rows:
            if self.config.from_data: self.load_data()
            else: self.load_images()

        # Still no rows?
        if self.no_rows: raise RuntimeError("No rows are loaded")

        # Initialize the figure
        self.initialize_figure()

        # Setup the plots
        #if self.single_grid:
            #self.plots = self.figure.create_grid(self.nrows, self.ncolumns)
            #print(self.nrows, self.ncolumns)
            #if self.config.coordinates: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white", projection=self.projection)
            #else: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white")
        #else: self.plots, self.colorbars = self.figure.create_row_of_image_grids(self.max_nrows_per_grid, self.ncolumns, self.config.ngrids, return_colorbars=True)

        #self.plots, self.colorbars = self.figure.create_row_of_image_grids(self.max_nrows_per_grid, self.ncolumns,
        #                                                                   self.config.ngrids, return_colorbars=True,
        #                                                                   share_colorbars=self.share_residuals_colorbars,
        #                                                                   adjust_grid=self.config.adjust_grid)

        #self.plots = self.figure.create_grid(self.)

        #print(self.projections)

        # Debugging
        log.debug("Creating the grid of plots ...")

        # Create plots
        self.plots = self.figure.create_row_of_grids(self.max_nrows_per_grid, self.ncolumns, self.config.ngrids,
                                                     last_column_not_shared_y=self.last_column_not_shared_y,
                                                     sharex=self.share_x, sharey=self.share_y,
                                                     projections=self.projections, share_per_row=self.share_per_row,
                                                     share_per_column=self.share_per_column)

        # Sort the frames on filter
        if self.config.sort_filters: self.sort()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading images from file ...")

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading data from file ...")

        # Observations
        self.load_observations()

        # Models
        self.load_models()

        # Residuals
        self.load_residuals()

        # Distributions
        self.load_distributions()

    # -----------------------------------------------------------------

    def load_observations(self):

        """
        This function ....
        :return:
        """

        # Debugging
        log.debug("Loading the observed frames ...")

        # Determine the path
        observations_path = self.input_path_directory("observations", check=True)

        # Load observations
        for path, name in fs.files_in_path(observations_path, extension="fits", returns=["path", "name"]):

            # Check the name
            if self.config.contains is not None and not strings.contains_any(name, self.config.contains): continue
            if self.config.not_contains is not None and strings.contains_any(name, self.config.not_contains): continue
            if self.config.exact_name is not None and name not in self.config.exact_name: continue
            if self.config.exact_not_name is not None and name in self.config.exact_not_name: continue
            if self.config.startswith is not None and not strings.startswith_any(name, self.config.startswith): continue
            if self.config.endswith is not None and not strings.endswith_any(name, self.config.endswith): continue

            # Add the observation
            self.add_observation_from_file(path)

    # -----------------------------------------------------------------

    def load_models(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the model frames ...")

        # Determine the path
        models_path = self.input_path_directory("models", check=True)

        # Load models
        for path in fs.files_in_path(models_path, extension="fits"): self.add_model_from_file(path, existing_row=True, only_if_row=True)

    # -----------------------------------------------------------------

    def load_residuals(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the residual frames ...")

        # Determine the path
        residuals_path = self.input_path_directory("residuals")

        # Load residuals
        if fs.is_directory(residuals_path):
            for path in fs.files_in_path(residuals_path, extension="fits"): self.add_residuals_from_file(path, existing_row=True, only_if_row=True)

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the residual distributions ...")

        # Determine the path
        distributions_path = self.input_path_directory("distributions")

        # Load distributions
        if fs.is_directory(distributions_path):
            for path in fs.files_in_path(distributions_path, extension="dat"): self.add_distribution_from_file(path, existing_row=True, only_if_row=True)

    # -----------------------------------------------------------------

    def sort(self):

        """
        This function ...
        :return:
        """

        # Sort on filter
        self.rows = ordered_by_value(self.rows, key=lambda entry: entry.observation.wavelength if entry.observation is not None else entry.model.wavelength)

    # -----------------------------------------------------------------

    @property
    def width(self):

        """
        This function ...
        :return:
        """

        #print(self.config.plot.xsize, self.ncolumns, self.mean_width_to_height)
        return self.config.plot.xsize * self.ncolumns * self.mean_width_to_height * self.config.ngrids

    # -----------------------------------------------------------------

    @property
    def height(self):

        """
        This function ...
        :return:
        """

        #print(self.config.plot.ysize, self.nrows)
        return self.config.plot.ysize * self.max_nrows_per_grid

    # -----------------------------------------------------------------

    @property
    def mean_width_to_height(self):

        """
        This function ...
        :return:
        """

        values = []

        # Loop over the frames
        for name in self.names: values.append(self.get_width_to_height(name))

        # Return the mean
        return numbers.arithmetic_mean(*values)

    # -----------------------------------------------------------------

    def get_width_to_height(self, name):

        """
        This function ...
        :return:
        """

        width = self.get_width(name)
        height = self.get_height(name)
        return float(width) / float(height)

    # -----------------------------------------------------------------

    def get_width(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_xsize(name)

    # -----------------------------------------------------------------

    def get_height(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_ysize(name)

    # -----------------------------------------------------------------

    def get_unit(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data
        observation, model = self.get_data(name)

        if observation is None and model is None: raise ValueError("No data for the '" + name + "' row")
        elif observation is None: return model.unit
        elif model is None: return observation.unit
        else:
            if model.unit != observation.unit: raise ValueError("Inconsistent units for '" + name + "' row")
            return observation.unit

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data
        observation, model = self.get_data(name)

        if observation is None and model is None: raise ValueError("No data for the '" + name + "' row")
        elif observation is None: return model.wcs
        elif model is None: return observation.wcs
        else:
            if model.wcs != observation.wcs: raise ValueError("Inconsistent coordinate systems for '" + name + "' row")
            return observation.wcs

    # -----------------------------------------------------------------

    def has_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_wcs(name) is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def all_have_wcs(self):

        """
        This function ...
        :return:
        """

        # Loop over the rows
        for name in self.names:
            if not self.has_wcs(name): return False
        return True

    # -----------------------------------------------------------------

    def get_shape(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data
        observation, model = self.get_data(name)

        if observation is None and model is None: raise ValueError("No data for the '" + name + "' row")
        elif observation is None: return model.shape
        elif model is None: return observation.shape
        else:
            if model.shape != observation.shape: raise ValueError("Inconsistent shape for '" + name + "' row")
            return observation.shape

    # -----------------------------------------------------------------

    def get_xsize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data
        observation, model = self.get_data(name)

        if observation is None and model is None: raise ValueError("No data for the '" + name + "' row")
        elif observation is None: return model.xsize
        elif model is None: return observation.xsize
        else:
            if observation.xsize != model.xsize: raise ValueError("Inconsistent xsize for '" + name + "' row")
            return observation.xsize

    # -----------------------------------------------------------------

    def get_ysize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data
        observation, model = self.get_data(name)

        if observation is None and model is None: raise ValueError("No data for the '" + name + "' row")
        elif observation is None: return model.ysize
        elif model is None: return observation.ysize
        else:
            if observation.ysize != model.ysize: raise ValueError("Inconsistent ysize for '" + name + "' row")
            return observation.ysize

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        return len(self.rows)

    # -----------------------------------------------------------------

    @property
    def nrows_per_grid(self):

        """
        This function ...
        :return:
        """

        return float(self.nrows) / self.config.ngrids

    # -----------------------------------------------------------------

    @property
    def max_nrows_per_grid(self):

        """
        This function ...
        :return:
        """

        return int(math.ceil(self.nrows_per_grid))

    # -----------------------------------------------------------------

    @property
    def min_nrows_per_grid(self):

        """
        This function ...
        :return:
        """

        return self.nrows // self.config.ngrids

    # -----------------------------------------------------------------

    @property
    def nremaining_rows(self):

        """
        This function ...
        :return:
        """

        return self.nrows % self.config.ngrids

    # -----------------------------------------------------------------

    @lazyproperty
    def first_row_grids(self):

        """
        This function ...
        :return:
        """

        rows = []
        for grid in range(self.config.ngrids):
            if grid < self.nremaining_rows: start = grid * (self.min_nrows_per_grid + 1);
            else: start = self.nremaining_rows * (self.min_nrows_per_grid + 1) + (grid - self.nremaining_rows) * self.min_nrows_per_grid
            rows.append(start)
        return rows

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_assignment(self):

        """
        This function ...
        :return:
        """

        assignment = []

        # Loop over the rows
        for row in range(self.nrows):

            if row < self.nremaining_rows * (self.min_nrows_per_grid + 1): grid = row // (self.min_nrows_per_grid + 1)
            else:
                tprime = row - self.nremaining_rows * (self.min_nrows_per_grid + 1)
                grid = self.nremaining_rows + tprime // self.min_nrows_per_grid
            assignment.append(grid)
        return assignment

    # -----------------------------------------------------------------

    def get_relative_row(self, row):

        """
        This function ...
        :param row:
        :return:
        """

        # Get grid for row
        grid = self.grid_assignment[row]

        # Get relative row in grid
        return row - self.first_row_grids[grid]

    # -----------------------------------------------------------------

    def get_grid_and_relative_row(self, row):

        """
        This function ...
        :param row:
        :return:
        """

        # Get grid for row
        grid = self.grid_assignment[row]

        # Get relative row in grid
        relative_row = row - self.first_row_grids[grid]

        # Return
        return grid, relative_row

    # -----------------------------------------------------------------

    @property
    def has_rows(self):

        """
        This function ...
        :return:
        """

        return self.nrows > 0

    # -----------------------------------------------------------------

    @property
    def no_rows(self):

        """
        This function ...
        :return:
        """

        return self.nrows == 0

    # -----------------------------------------------------------------

    @property
    def nmasks(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.masks)

    # -----------------------------------------------------------------

    @property
    def has_masks(self):

        """
        This function ...
        :return:
        """

        return self.nmasks > 0

    # -----------------------------------------------------------------

    @property
    def no_masks(self):

        """
        This function ...
        :return:
        """

        return self.nmasks == 0

    # -----------------------------------------------------------------

    @property
    def nregions(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.regions)

    # -----------------------------------------------------------------

    @property
    def has_regions(self):

        """
        This function ...
        :return:
        """

        return self.nregions > 0

    # -----------------------------------------------------------------

    @property
    def no_regions(self):

        """
        This function ...
        :return:
        """

        return self.nregions == 0

    # -----------------------------------------------------------------

    def nregions_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.names: raise ValueError("Invalid row '" + name + "'")

        if self.rows[name].regions is None: return 0
        else: return len(self.rows[name].regions)

    # -----------------------------------------------------------------

    def has_regions_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nregions_for_row(name) > 0

    # -----------------------------------------------------------------

    def no_regions_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nregions_for_row(name) == 0

    # -----------------------------------------------------------------

    def nmasks_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.names: raise ValueError("Invalid row '" + name + "'")

        if self.rows[name].masks is None: return 0
        else: return len(self.rows[name].masks)

    # -----------------------------------------------------------------

    def has_masks_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nmasks_for_row(name) > 0

    # -----------------------------------------------------------------

    def no_masks_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.nmasks_for_row(name) == 0

    # -----------------------------------------------------------------

    @lazyproperty
    def with_residuals_row_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in self.rows:
            entry = self.rows[name]
            if not entry.residuals: continue
            names.append(name)
        return names

    # -----------------------------------------------------------------

    def do_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.with_residuals_row_names

    # -----------------------------------------------------------------

    @property
    def nwith_residuals_rows(self):

        """
        This function ...
        :return:
        """

        return len(self.with_residuals_row_names)

    # -----------------------------------------------------------------

    @property
    def has_residual_rows(self):

        """
        This function ...
        :return:
        """

        return self.nwith_residuals_rows > 0

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        if not self.has_residual_rows: return 2
        else:
            if self.config.distributions: return 4
            else: return 3

    # -----------------------------------------------------------------

    @property
    def last_column_not_shared_y(self):

        """
        This function ...
        :return:
        """

        return self.config.distributions

    # -----------------------------------------------------------------

    @property
    def first_name(self):

        """
        This function ...
        :return:
        """

        return self.names[0]

    # -----------------------------------------------------------------

    @property
    def last_name(self):

        """
        This function ...
        :return:
        """

        return self.names[-1]

    # -----------------------------------------------------------------

    @property
    def first_wcs(self):

        """
        This function ...
        :return:
        """

        return self.get_wcs(self.first_name)

    # -----------------------------------------------------------------

    @property
    def projection(self):

        """
        This function ...
        :return:
        """

        return self.first_wcs

    # -----------------------------------------------------------------

    def create_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the residual frames ...")

        # Loop over the rows for which residuals have to be calculated
        for name in self.with_residuals_row_names:

            # Check
            if not self.has_data(name):
                log.warning("No data for the '" + name + "' row to create residuals: skipping ...")
                continue
                #raise ValueError("No data for the '" + name + "' row to create residuals")

            # Get the observation and model
            observation = self.get_observation(name)
            model = self.get_model(name)
            errors = self.get_errors(name)

            # Error-weighed residuals
            if self.config.weighed: residual = (model - observation) / errors

            # Absolute residuals
            elif self.config.relative: residual = (model - observation) / model

            # Relative residuals
            else: residual = model - observation

            # ABSOLUTE VALUES?
            if self.config.absolute: residual = residual.absolute

            # Set the coordinate system
            residual.wcs = self.get_wcs(name)

            # Add the residual frame
            self.residuals[name] = residual

    # -----------------------------------------------------------------

    @property
    def residuals_names(self):

        """
        Thisn function ...
        :return:
        """

        return self.residuals.keys()

    # -----------------------------------------------------------------

    def get_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.residuals_names: raise ValueError("No residuals map for the '" + name + "' row")
        return self.residuals[name]

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the residual distributions ...")

        # Loop over the residual maps
        for name in self.residuals_names:

            # Debugging
            log.debug("Creating distribution for the '" + name + "' residuals ...")

            # Get the residual map
            residuals = self.get_residuals(name)

            # Create the distribution
            distribution = Distribution.from_data("Residual", residuals, sigma_clip=self.config.sigma_clip_distributions, sigma_level=self.config.sigma_clip_level)

            # Add the distribution
            self.distributions[name] = distribution

    # -----------------------------------------------------------------

    @property
    def distribution_names(self):

        """
        This function ...
        :return:
        """

        return self.distributions.keys()

    # -----------------------------------------------------------------

    def get_distribution(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.distribution_names: raise ValueError("No residuals distribution for the '" + name + "' row")
        return self.distributions[name]

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the observation
        self.write_observations()

        # Write the models
        self.write_models()

        # Write the regions
        if self.has_regions: self.write_regions()

        # Write the masks
        if self.has_masks: self.write_masks()

        # Write the residual frames
        self.write_residuals()

        # Write the distributions
        if self.config.distributions: self.write_distributions()

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in self.names:
            if not self.has_observation_data(name): continue
            names.append(name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def observations_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("observations", create=True)

    # -----------------------------------------------------------------

    def write_observations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the observation frames ...")

        # Loop over the frames
        for name in self.observation_names:

            # Get the observation
            observation = self.get_observation(name)

            # Determine path
            path = fs.join(self.observations_path, name + ".fits")

            # Save the observation image
            observation.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in self.names:
            if not self.has_model_data(name): continue
            names.append(name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def models_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("models", create=True)

    # -----------------------------------------------------------------

    def write_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model frames ...")

        # Loop over the frames
        for name in self.model_names:

            # Get the model
            model = self.get_model(name)

            # Determine the path
            path = fs.join(self.models_path, name + ".fits")

            # Save the model image
            model.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def regions_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("regions", create=True)

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the regions ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def masks_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("masks", create=True)

    # -----------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing masks ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("residuals", create=True)

    # -----------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual frames ...")

        # Loop over the residuals
        for name in self.residuals_names:

            # Get the residuals frame
            residuals = self.get_residuals(name)

            # Determine path
            path = fs.join(self.residuals_path, name + ".fits")

            # Write the residuals map
            residuals.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def distributions_path(self):
        
        """
        This function ...
        :return: 
        """

        return self.output_path_directory("distributions", create=True)
    
    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual distributions ...")

        # Loop over the distributions
        for name in self.distribution_names:

            # Get the distribution
            distribution = self.get_distribution(name)

            # Determine the path
            path = fs.join(self.distributions_path, name + ".dat")

            # Write the distribution
            distribution.saveto(path)

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.rows.keys()

    # -----------------------------------------------------------------

    def index_for_row(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.names.index(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_row_name(self):

        """
        This function ...
        :return:
        """

        if self.config.scale_reference is None: return None
        elif self.config.scale_reference not in self.names: raise ValueError("Invalid reference row name")
        else: return self.config.scale_reference

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_reference_row_name(self):

        """
        This function ...
        :return:
        """

        if self.config.scale_residuals_reference is None: return None
        elif self.config.scale_residuals_reference not in self.names: raise ValueError("Invalid residuals reference row name")
        else: return self.config.scale_residuals_reference

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_row_index(self):

        """
        This function ...
        :return:
        """

        return self.index_for_row(self.reference_row_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_reference_row_index(self):

        """
        This function ...
        :return:
        """

        return self.index_for_row(self.residuals_reference_row_name)

    # -----------------------------------------------------------------

    def plot_reference(self):

        """
        This function ...
        :return:
        """

        # Plot the row
        return self.plot_row(self.reference_row_name, self.reference_row_index)

    # -----------------------------------------------------------------

    def plot_residuals_reference(self):

        """
        This function ...
        :return:
        """

        # Plot the row
        return self.plot_row(self.residuals_reference_row_name, self.residuals_reference_row_index)

    # -----------------------------------------------------------------

    @property
    def adjust_grid(self):

        """
        This function ...
        :return:
        """

        if self.config.adjust_grid is not None: return self.config.adjust_grid
        elif self.same_shapes: return True
        elif self.config.uniformize: return True
        else: return False  # if grid is adjusted when plotting images with different shapes, it looks bad

    # -----------------------------------------------------------------

    def plot_row(self, name, index):

        """
        This function ...
        :param name:
        :param index:
        :return:
        """

        # Check
        if self.is_plotted(name): log.warning("The '" + name + "' row has already been plotted")

        # Debugging
        log.debug("Adding the '" + name + "' row to the plot ...")

        # Get properties
        unit = self.get_unit(name)
        xsize = self.get_xsize(name)
        ysize = self.get_ysize(name)

        # Debugging
        log.debug("Dimensions of row images: " + str(xsize) + " x " + str(ysize))
        log.debug("Unit of row: " + str(unit))

        # Determine image scale
        if self.config.share_scale: vmin, vmax = self._vmin, self._vmax
        else: vmin = vmax = None

        # Plot the observation
        vmin_observation, vmax_observation, obs_image, obs_normalization = self.plot_observation(name, index, vmin=vmin, vmax=vmax)

        # Plot the model, on the same scale
        self.plot_model(name, index, vmin=vmin_observation, vmax=vmax_observation, set_minmax=False)

        # Determine residuals scale # NO, WHY WOULD WE SHOW THE RESIDUALS ON THE SCALE OF THE IMAGES?
        #if self.config.same_residuals_scale: vmin_res, vmax_res = vmin_observation, vmax_observation
        #else: vmin_res = vmax_res = None

        # Plot the residuals
        if self.do_residuals(name): vmin_res, vmax_res, res_image, res_normalization = self.plot_residuals(name, index) # vmin=vmin_res, vmax=vmax_res)
        else: vmin_res = vmax_res = res_image = res_normalization = None

        # Plot the distribution
        if self.config.distributions: self.plot_distribution(name, index, vmin=vmin_res, vmax=vmax_res)

        # Plot the colorbar of the residual
        # if self.share_residuals_colorbars: pass
        # elif name == "0":
        #     colorbars = self.colorbars[index]
        #     for i, colorbar in enumerate(colorbars):
        #         #if i == 0: continue # FIRST ROW
        #         #if i == 1: continue # SECOND ROW
        #         self.figure.figure.colorbar(res_image, cax=colorbar)
        #     #colorbar = self.colorbars[index]
        #     #self.figure.figure.colorbar(res_image, cax=colorbar)
        #     #colorbar.get_xaxis().set_ticks([])
        #     #colorbar.get_yaxis().set_ticks([])
        #     #colorbar.set_axis_off()

        # Add to plotted row names
        self._plotted_rows.append(name)

    # -----------------------------------------------------------------

    @property
    def column_indices(self):

        """
        This function ...
        :return:
        """

        return list(range(self.ncolumns)) # [0, 1, 2] or [0, 1, 2, 3]

    # -----------------------------------------------------------------

    @property
    def column_indices_not_distributions(self):

        """
        This function ...
        :return:
        """

        return [0, 1, 2]

    # -----------------------------------------------------------------

    def get_plot(self, abs_row_index, col_index):

        """
        This function ...
        :param abs_row_index:
        :param col_index:
        :return:
        """

        grid, rel_row_index = self.get_grid_and_relative_row(abs_row_index)
        #print(abs_row_index, col_index, grid, rel_row_index)
        return self.plots[grid][rel_row_index][col_index]

    # -----------------------------------------------------------------

    def get_colorbar(self, abs_row_index, col_index):

        """
        This function ...
        :param abs_row_index:
        :param col_index:
        :return:
        """

        grid, rel_row_index = self.get_grid_and_relative_row(abs_row_index)
        return self.colorbars[grid][rel_row_index][col_index]

    # -----------------------------------------------------------------

    def get_observation_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._observation_images[name]

    # -----------------------------------------------------------------

    def get_observation_normalization(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._observation_normalizations[name]

    # -----------------------------------------------------------------

    def get_model_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._model_images[name]

    # -----------------------------------------------------------------

    def get_model_normalization(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._model_normalizations[name]

    # -----------------------------------------------------------------

    def get_residual_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._residual_images[name]

    # -----------------------------------------------------------------

    def get_residual_normalization(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._residual_normalizations[name]

    # -----------------------------------------------------------------

    def get_distribution_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._distribution_images[name]

    # -----------------------------------------------------------------

    def plot_empty_row(self, name, index):

        """
        This function ...
        :param name:
        :param index:
        :return:
        """

        # Debugging
        log.debug("Formatting empty row '" + name + "' ...")

        # Loop over the columns
        for col in self.column_indices:

            # Get the plot
            plot = self.get_plot(index, col)

            # Color spines
            plot.axes.spines['bottom'].set_color("white")
            plot.axes.spines['top'].set_color("white")
            plot.axes.spines['left'].set_color("white")
            plot.axes.spines['right'].set_color("white")

            # Color ticks
            # plot.axes.xaxis.label.set_color("white")
            # plot.axes.yaxis.label.set_color("white")
            plot.axes.tick_params(axis='x', colors="white", direction="inout")
            plot.axes.tick_params(axis='y', colors="white", direction="inout")

            # Set background color: otherwise NaNs are not plotted (-> white/transparent)
            if self.config.background: plot.axes.set_facecolor(self.background_color)

            # Add the label
            if col == 0: plot.axes.text(0.95, 0.95, name, color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")

    # -----------------------------------------------------------------

    def get_mask(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        return self.rows[name].masks[label]

    # -----------------------------------------------------------------

    def get_regions(self, name, label):

        """
        This function ...
        :param name:
        :param label:
        :return:
        """

        return self.rows[name].regions[label]

    # -----------------------------------------------------------------

    @property
    def colormap_name(self):

        """
        This function ...
        :return:
        """

        return self.config.colormap

    # -----------------------------------------------------------------

    @lazyproperty
    def colormap(self):

        """
        This function ...
        :return:
        """

        return cm.get_cmap(self.colormap_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def background_color(self):

        """
        This function ...
        :return:
        """

        return self.colormap(0.0)

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_colormap_name(self):

        """
        This function ...
        :return:
        """

        # Colormap for residuals is specified
        if self.config.residuals_colormap is not None: return cm.get_cmap(self.config.residuals_colormap)

        # Not specified, plotting absolute values of residuals
        elif self.config.absolute: return "BuPu" # sequential colormaps

        # Plotting relative or weighed residuals of residuals
        else: return "seismic" # diverging colormaps

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_colormap(self):

        """
        This function ...
        :return:
        """

        return cm.get_cmap(self.residuals_colormap_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_background_color(self):

        """
        This function ...
        :return:
        """

        return self.residuals_colormap(0.0)

    # -----------------------------------------------------------------

    def _plot_frame(self, name, frame, row, col, vmin=None, vmax=None, add_label=True, colormap=None, scale=None):

        """
        This function ...
        :param name:
        :param frame:
        :param row:
        :param col:
        :param vmin:
        :param vmax:
        :param add_label:
        :param colormap:
        :param scale:
        :return:
        """

        # Get the plot
        plot = self.get_plot(row, col)

        # Color spines
        plot.axes.spines['bottom'].set_color("white")
        plot.axes.spines['top'].set_color("white")
        plot.axes.spines['left'].set_color("white")
        plot.axes.spines['right'].set_color("white")

        # Color ticks
        # plot.axes.xaxis.label.set_color("white")
        # plot.axes.yaxis.label.set_color("white")
        plot.axes.tick_params(axis='x', colors="white", direction="inout")
        plot.axes.tick_params(axis='y', colors="white", direction="inout")

        # Set background color: otherwise NaNs are not plotted (-> white/transparent)
        if self.config.background: plot.axes.set_facecolor(self.background_color)

        # Add mask if present
        if self.has_masks_for_row(name):
            for label in self.rows[name].masks:
                mask = self.get_mask(name, label)
                frame[mask] = nan

        # Set interval
        if vmin is not None and vmax is not None: interval = [vmin, vmax]
        elif (vmin is not None) or (vmax is not None): raise ValueError("Must pass vmin and vmax together")
        else: interval = self.config.interval

        #else: aspect = "equal"
        #aspect = (1., 1.)
        #aspect = "equal"
        #aspect = "auto"
        if self.use_coordinates: aspect = "auto"
        else: aspect = "equal"

        # Determine the scale
        if scale is None: scale = self.config.scale # take the configured default when none is passed

        # Determine the colormap
        if colormap is None: colormap = self.colormap # take the configured default when none is passed

        # Plot
        #plot.axes.set_adjustable('box-forced')
        vmin_image, vmax_image, image, normalization = plotting.plot_box(frame.data, axes=plot.axes, interval=interval,
                                                                         scale=scale, cmap=colormap,
                                                                         alpha=self.config.alpha, return_image=True,
                                                                         return_normalization=True, aspect=aspect)

        # Add region if present
        if self.has_regions_for_row(name):
            for label in self.rows[name].regions:
                regions = self.get_regions(name, label)
                for patch in regions.to_mpl_patches(): plot.axes.add_patch(patch)

        # Add the label
        if add_label: plot.axes.text(0.95, 0.95, name, color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")

        # Return
        return vmin_image, vmax_image, image, normalization

    # -----------------------------------------------------------------

    def _plot_residuals(self, name, residuals, row, col, colormap=None, scale=None, vmin=None, vmax=None):

        """
        This function ...
        :param name:
        :param residuals:
        :param row:
        :param col:
        :param colormap:
        :param scale:
        :param vmin:
        :param vmax:
        :return:
        """

        # Get the plot
        plot = self.get_plot(row, col)

        # Color spines
        plot.axes.spines['bottom'].set_color("white")
        plot.axes.spines['top'].set_color("white")
        plot.axes.spines['left'].set_color("white")
        plot.axes.spines['right'].set_color("white")
        plot.axes.tick_params(axis='x', colors="white", direction="inout")
        plot.axes.tick_params(axis='y', colors="white", direction="inout")

        # Set background color: otherwise NaNs are not plotted (-> white/transparent)
        #if self.config.background: plot.axes.set_facecolor(self.residuals_background_color)

        # Add mask if present
        if self.has_masks_for_row(name):
            for label in self.rows[name].masks:
                mask = self.get_mask(name, label)
                residuals[mask] = nan

        #aspect = "equal"
        if self.use_coordinates: aspect = "auto"
        else: aspect = "equal"

        # Set interval
        if vmin is not None and vmax is not None: interval = [vmin, vmax]
        elif (vmin is not None) or (vmax is not None): raise ValueError("Must pass vmin and vmax together")
        else: interval = self.config.residuals_interval

        # Determine the scale
        if scale is None: scale = "linear"

        # Determine the colormap
        if colormap is None: colormap = self.residuals_colormap

        # Set flags
        if self.config.absolute: around_zero = symmetric = False
        else: around_zero = symmetric = True

        #print("interval", interval)
        #print("scale", scale)
        #print("colormap", colormap)
        #print("around_zero", around_zero)
        #print("symmetric", symmetric)

        # Plot
        vmin_image, vmax_image, image, normalization = plotting.plot_box(residuals.data, axes=plot.axes, interval=interval,
                                                                         scale=scale, cmap=colormap,
                                                                         alpha=self.config.alpha, return_image=True,
                                                                         return_normalization=True, aspect=aspect,
                                                                         around_zero=around_zero, symmetric=symmetric,
                                                                         check_around_zero=False)

        #print("limits res", vmin_image, vmax_image)

        # Add region if present
        if self.config.regions_on_residuals:
            if self.has_regions_for_row(name):
                for label in self.rows[name].regions:
                    regions = self.get_regions(name, label)
                    for patch in regions.to_mpl_patches(): plot.axes.add_patch(patch)

        # Return
        return vmin_image, vmax_image, image, normalization

    # -----------------------------------------------------------------

    def plot_observation(self, name, index, vmin=None, vmax=None, set_minmax=None):

        """
        This function ...
        :param name:
        :param index:
        :param vmin:
        :param vmax:
        :param set_minmax:
        :return:
        """

        # Debugging
        log.debug("Plotting the '" + name + "' observation image ...")

        # Get the observation
        observation = self.get_observation(name)

        # Set the normalization
        if vmin is None and vmax is None and self.config.share_scale: vmin, vmax = self._vmin, self._vmax

        # Plot
        vmin_image, vmax_image, image, normalization = self._plot_frame(name, observation, index, observation_index, vmin=vmin, vmax=vmax, add_label=True)

        # Set the image and normalization
        self._observation_images[name] = image
        self._observation_normalizations[name] = normalization

        # Set vmin and vmax
        if set_minmax is None: set_minmax = self.config.share_scale
        if set_minmax: self._vmin, self._vmax = vmin_image, vmax_image

        # Return
        return vmin_image, vmax_image, image, normalization

    # -----------------------------------------------------------------

    def plot_model(self, name, index, vmin=None, vmax=None, set_minmax=None):

        """
        This function ...
        :param name:
        :param index:
        :param vmin:
        :param vmax:
        :param set_minmax:
        :return:
        """

        # Debugging
        log.debug("Plotting the '" + name + "' model image ...")

        # Get the model
        model = self.get_model(name)

        # Set the normalization if necessary
        if vmin is None and vmax is None and self.config.share_scale: vmin, vmax = self._vmin, self._vmax

        # Plot
        vmin_image, vmax_image, image, normalization = self._plot_frame(name, model, index, model_index, vmin=vmin, vmax=vmax, add_label=False)

        # Set the image and the normalization
        self._model_images[name] = image
        self._model_normalizations[name] = normalization

        # Set vmin and vmax
        if set_minmax is None: set_minmax = self.config.share_scale
        if set_minmax: self._vmin, self._vmax = vmin_image, vmax_image

        # Return
        return vmin_image, vmax_image, image, normalization

    # -----------------------------------------------------------------

    def plot_residuals(self, name, index, vmin=None, vmax=None, set_minmax=None):

        """
        This function ...
        :param name:
        :param index:
        :param vmin:
        :param vmax:
        :param set_minmax:
        :return:
        """

        # Debugging
        log.debug("Plotting the '" + name + "' residual map ...")

        # Get the residual map
        residuals = self.get_residuals(name)
        #distribution = Distribution.from_data("residuals", residuals)
        #plot_distribution(distribution)

        # Set the vmin and vmax if necessary
        if vmin is None and vmax is None and self.config.share_scale_residuals: vmin, vmax = self._vmin_res, self._vmax_res

        #print(vmin, vmax)

        # Plot
        vmin_image, vmax_image, image, normalization = self._plot_residuals(name, residuals, index, residuals_index, vmin=vmin, vmax=vmax)

        # Set the residual image and the normalization
        self._residual_images[name] = image
        self._residual_normalizations[name] = normalization

        # Set min and max
        if set_minmax is None: set_minmax = self.config.share_scale_residuals
        if set_minmax: self._vmin_res, self._vmax_res = vmin_image, vmax_image

        # Plot colorbar
        # Get the colorbar axes
        #colorbar = self.get_colorbar(index, residuals_index)
        #self.figure.figure.colorbar(image, cax=colorbar)

        # Return
        return vmin_image, vmax_image, image, normalization

    # -----------------------------------------------------------------

    def plot_distribution(self, name, index, vmin=None, vmax=None):

        """
        This function ...
        :param name:
        :param index:
        :param vmin:
        :param vmax:
        :return:
        """

        # Debugging
        log.debug("Plotting the '" + name + "' residuals distribution ...")

        # Get the distribution
        distribution = self.get_distribution(name)

        # Get the plot
        plot = self.get_plot(index, distribution_index)

        # Color spines
        # print(plot.axes.spines.keys())
        plot.axes.spines['bottom'].set_color("white")
        plot.axes.spines['top'].set_color("white")
        plot.axes.spines['left'].set_color("white")
        plot.axes.spines['right'].set_color("white")

        # Color ticks
        # plot.axes.xaxis.label.set_color("white")
        # plot.axes.yaxis.label.set_color("white")
        plot.axes.tick_params(axis='x', colors="white", direction="inout")
        plot.axes.tick_params(axis='y', colors="white", direction="inout")

        # Define the colors
        normalization = self.get_residual_normalization(name)
        colors = []
        for value in distribution.values:
            #print("value", value)
            normalized = normalization(value)
            #print("normalized", normalized)
            rgb = self.residuals_colormap(normalized)[0]
            #print("color", rgb)
            hex = rgb2hex(rgb)
            #print(hex)
            colors.append(hex)

        # Set distribution x limits
        if vmin is not None and vmax is not None: x_limits = (vmin, vmax)
        elif (vmin is not None) or (vmax is not None): raise ValueError("Cannot specify only vmin or vmax")
        else: x_limits = None

        print("x_limits", x_limits)

        # Plot the distribution
        image = plot_distribution(distribution, axes=plot.axes, colors=colors, x_limits=x_limits)

        # Set the image
        self._distribution_images[index] = image

    # -----------------------------------------------------------------

    def get_data(self, name):

        """
        This function ....
        :param name:
        :return:
        """

        return self.get_observation(name), self.get_model(name)

    # -----------------------------------------------------------------

    def get_observation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.rows[name].observation

    # -----------------------------------------------------------------

    def get_observation_unit(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation(name): return self.get_observation(name).unit
        else: return None

    # -----------------------------------------------------------------

    def get_observation_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation(name): return self.get_observation(name).unit
        else: return None

    # -----------------------------------------------------------------

    def get_observation_shape(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation_data(name): return self.get_observation(name).shape
        else: return None

    # -----------------------------------------------------------------

    def get_observation_xsize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation(name): return self.get_observation(name).xsize
        else: return None

    # -----------------------------------------------------------------

    def get_observation_ysize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation(name): return self.get_observation(name).ysize
        else: return None

    # -----------------------------------------------------------------

    def get_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.rows[name].model

    # -----------------------------------------------------------------

    def get_model_unit(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_model(name): return self.get_model(name).unit
        else: return None

    # -----------------------------------------------------------------

    def get_model_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_model(name): return self.get_model(name).wcs
        else: return None

    # -----------------------------------------------------------------

    def get_model_shape(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_model(name): return self.get_model(name).shape
        else: return None

    # -----------------------------------------------------------------

    def get_model_xsize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_model(name): return self.get_model(name).xsize
        else: return None

    # -----------------------------------------------------------------

    def get_model_ysize(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_model(name): return self.get_model(name).ysize
        else: return None

    # -----------------------------------------------------------------

    def get_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.rows[name].errors

    # -----------------------------------------------------------------

    def has_observation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_row(name) and self.rows[name].observation is not None

    # -----------------------------------------------------------------

    def has_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_row(name) and self.rows[name].model is not None

    # -----------------------------------------------------------------

    def has_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.residuals and self.residuals[name] is not None

    # -----------------------------------------------------------------

    def has_observation_data(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        frame = self.get_observation(name)
        return not (frame is None or frame.all_zeroes or frame.all_nans or frame.all_infs)

    # -----------------------------------------------------------------

    def has_model_data(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        frame = self.get_model(name)
        return not (frame is None or frame.all_zeroes or frame.all_nans or frame.all_infs)

    # -----------------------------------------------------------------

    def has_data(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_observation_data(name) and self.has_model_data(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def scale_reference(self):

        """
        This function ...
        :return:
        """

        if self.config.share_scale: return self.config.scale_reference
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_scale_reference(self):

        """
        Thisfunction ...
        :return:
        """

        return self.scale_reference is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def scale_residuals_reference(self):

        """
        This function ...
        :return:
        """

        if self.config.share_scale_residuals: return self.config.scale_residuals_reference
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_scale_residuals_reference(self):

        """
        This function ...
        :return:
        """

        return self.scale_residuals_reference is not None

    # -----------------------------------------------------------------

    @property
    def one_scale_reference(self):

        """
        This function ...
        :return:
        """

        if self.has_scale_reference and self.has_scale_residuals_reference: return self.scale_reference == self.scale_residuals_reference
        elif not self.has_scale_reference and not self.has_scale_residuals_reference: return False
        else: return True

    # -----------------------------------------------------------------

    def plot_references(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting the reference row(s) ...")

        #print(self.has_scale_reference, self.has_scale_residuals_reference)

        # Image and residual map reference rows are the same
        if self.has_scale_reference and self.has_scale_residuals_reference and self.one_scale_reference: self.plot_reference()

        # Image and residual map reference rows are different
        elif self.has_scale_reference and self.has_scale_residuals_reference:

            self.plot_reference()
            self.plot_residuals_reference()

        # Reference row for the images
        elif self.has_scale_reference: self.plot_reference()

        # Reference row for the residual maps
        elif self.has_scale_residuals_reference: self.plot_residuals_reference()

    # -----------------------------------------------------------------

    def is_plotted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self._plotted_rows

    # -----------------------------------------------------------------

    @lazyproperty
    def min_xsize(self):

        """
        This function ...
        :return:
        """

        value = None
        for name in self.names:
            xsize = self.get_xsize(name)
            if value is None or xsize < value: value = xsize
        return value

    # -----------------------------------------------------------------

    @lazyproperty
    def max_xsize(self):

        """
        This function ...
        :return:
        """

        value = None
        for name in self.names:
            xsize = self.get_xsize(name)
            if value is None or xsize > value: value = xsize
        return value

    # -----------------------------------------------------------------

    @lazyproperty
    def min_ysize(self):

        """
        This function ...
        :return:
        """

        value = None
        for name in self.names:
            ysize = self.get_ysize(name)
            if value is None or ysize < value: value = ysize
        return value

    # -----------------------------------------------------------------

    @lazyproperty
    def max_ysize(self):

        """
        This function ...
        :return:
        """

        value = None
        for name in self.names:
            ysize = self.get_ysize(name)
            if value is None or ysize > value: value = ysize
        return value

    # -----------------------------------------------------------------

    @lazyproperty
    def same_shapes(self):

        """
        Thisn function ...
        :return:
        """

        shapes = []
        for name in self.names:
            shape = self.get_shape(name)
            shapes.append(shape)
        return sequences.all_equal(shapes)

    # -----------------------------------------------------------------

    @property
    def needs_uniform(self):

        """
        This function ...
        :return:
        """

        return self.config.uniformize and not self.same_shapes

    # -----------------------------------------------------------------

    def uniformize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uniformizing the rows ...")

        #print("min xsize:", self.min_xsize)
        #print("max xsize:", self.max_xsize)
        #print("min ysize:", self.min_ysize)
        #print("max ysize:", self.max_ysize)

        # Loop over the rows
        for name in self.names:

            # Get the xsize and ysize
            xsize = self.get_xsize(name)
            ysize = self.get_ysize(name)

            # Check size
            if xsize < self.max_xsize and ysize < self.max_ysize:

                # Get x and y factor
                x_factor = float(self.max_xsize) / xsize
                y_factor = float(self.max_ysize) / ysize

                # Determine factor
                if x_factor > y_factor: factor = int(math.ceil(x_factor))
                else: factor = int(math.ceil(y_factor))

                # Upsample
                self.upsample_row(name, factor)

            # Check xsize
            elif xsize < self.max_xsize:

                # Determine the factor
                factor = float(self.max_xsize) / xsize
                factor = int(math.ceil(factor))

                # Usample
                self.upsample_row(name, factor)

            # Check ysize
            elif ysize < self.max_ysize:

                # Determine the factor
                factor = float(self.max_ysize) / ysize
                factor = int(math.ceil(factor))

                # Usample
                self.upsample_row(name, factor)

    # -----------------------------------------------------------------

    def upsample_row(self, name, factor):

        """
        This function ...
        :param name:
        :param factor:
        :return:
        """

        # Debugging
        log.debug("Upsampling the '" + name + "' row with a factor of " + str(factor) + " ...")

        # Observation
        if self.has_observation(name):

            frame = self.get_observation(name)
            frame.upsample(factor, order=0)

        # Model
        if self.has_model(name):

            frame = self.get_model(name)
            frame.upsample(factor, order=0)

        # Residual map
        if self.has_residuals(name):

            frame = self.get_residuals(name)
            frame.upsample(factor, order=0)

    # -----------------------------------------------------------------

    def plot_aplpy(self):

        """
        This function ...
        :return:
        """

        import aplpy

        left = [0.1, 0.3, 0.5, 0.7]
        width = [0.2, 0.2, 0.2, 0.2]

        rectLS = []
        for x in left:
            for y in left:
                rectLS.append([x, y, 0.2, 0.2])

        axLS = []

        #fig = plt.figure()

        # First row
        #axLS.append(fig.add_axes(rectLS[0]))
        #for i in [1, 2, 3]: axLS.append(fig.add_axes(rectLS[i], sharey=axLS[-1]))

        # Second row
        #axLS.append(fig.add_axes(rectLS[4]))
        #for i in [1, 2, 3]: axLS.append(fig.add_axes(rectLS[i + 4], sharex=axLS[i], sharey=axLS[-1]))

        # Third row
        #axLS.append(fig.add_axes(rectLS[8]))
        #for i in [5, 6, 7]: axLS.append(fig.add_axes(rectLS[i + 4], sharex=axLS[i], sharey=axLS[-1]))

        # Fourth row
        #axLS.append(fig.add_axes(rectLS[12]))
        #for i in [9, 10, 11]: axLS.append(fig.add_axes(rectLS[i + 4], sharex=axLS[i], sharey=axLS[-1]))

        # Loop over the images
        for index, name in enumerate(self.names):

            # Check
            if self.is_plotted(name): continue

            # Has observation and model data
            #if self.has_data(name): self.plot_row(name, index)

            # Only observation data
            #elif self.has_observation_data(name): self.plot_observation(name, index)

            # Only model data
            #elif self.has_model_data(name): self.plot_model(name, index)

            frame = self.get_observation(name)

            #f1 = aplpy.FITSFigure(frame.wcs, figure=fig, subplot=plotloc[0])
            f1 = aplpy.FITSFigure(frame.hdu, figure=self.figure.figure, subplot=rectLS[index])
            #standard_setup(f1)
            ## f1.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
            f1.show_colorscale(vmin=0, vmax=0.00025, cmap='hot')
            f1.add_beam(major=0.01, minor=0.01, angle=0, fill=True, color='white')

            # f1.axis_labels.show_y()
            #f1.tick_labels.set_xposition('top')
            #f1.tick_labels.show()

        plt.show()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot using APLpy
        if self.config.aplpy: self.plot_aplpy()

        # Plot using imagegrid
        else: self.plot_imagegrid()

    # -----------------------------------------------------------------

    def plot_imagegrid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting using ImageGrid ...")

        # Plot the references
        self.plot_references()

        # Loop over the images
        for index, name in enumerate(self.names):

            # Check
            if self.is_plotted(name): continue

            # Debugging
            log.debug("Plotting the '" + name + "' row ...")

            # Has observation and model data
            if self.has_data(name): self.plot_row(name, index)

            # Only observation data
            elif self.has_observation_data(name): self.plot_observation(name, index)

            # Only model data
            elif self.has_model_data(name): self.plot_model(name, index)

            # No data
            else: self.plot_empty_row(name, index)

        # Plot empty rows
        self.plot_empty_rows()

        # Plot colorbars
        self.plot_colorbars()

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    def plot_colorbars(self):

        """
        This function ...
        :return:
        """

        # Plot colorbar
        # Get the colorbar axes
        # colorbar = self.get_colorbar(index, residuals_index)
        # self.figure.figure.colorbar(image, cax=colorbar)

        pass

    # -----------------------------------------------------------------

    def plot_empty_rows(self):

        """
        This function ...
        :return:
        """

        plotted = [[] for _ in range(self.config.ngrids)]

        # Loop over the images
        for index, name in enumerate(self.names):
            grid, rel_row_index = self.get_grid_and_relative_row(index)
            plotted[grid].append(rel_row_index)

        # Loop over the empty rows
        for grid_index, plotted_rows in enumerate(plotted):

            for relative_row in range(self.max_nrows_per_grid):
                if relative_row in plotted_rows: continue

                # Loop over the columns
                for col in self.column_indices:

                    # Get the plot
                    plot = self.plots[grid_index][relative_row][col]

                    # Color spines
                    plot.axes.spines['bottom'].set_color("white")
                    plot.axes.spines['top'].set_color("white")
                    plot.axes.spines['left'].set_color("white")
                    plot.axes.spines['right'].set_color("white")

                    # Color ticks
                    # plot.axes.xaxis.label.set_color("white")
                    # plot.axes.yaxis.label.set_color("white")
                    plot.axes.tick_params(axis='x', colors="white", direction="inout")
                    plot.axes.tick_params(axis='y', colors="white", direction="inout")

                    # Set background color: otherwise NaNs are not plotted (-> white/transparent)
                    #if self.config.background: plot.axes.set_facecolor(self.background_color)

    # -----------------------------------------------------------------

    def _plot_colorbar(self, image, grid_index):

        """
        This function ...
        :param image:
        :param grid_index:
        :return:
        """

        # Get the colorbar
        colorbar = self.colorbars[grid_index]

        # Plot
        self.figure.figure.colorbar(image, cax=colorbar)

        # Set properties
        colorbar.get_xaxis().set_ticks([])
        colorbar.get_yaxis().set_ticks([])
        colorbar.set_axis_off()

# -----------------------------------------------------------------

def sort_numbs(arr):
  numbers = []
  for k in range(len(arr)): 
    numb = str(arr[k].split('/')[-1].split('_')[-1].split('.fits'))
    #print numb
    numbers.append(numb)
  a = sorted(numbers)
  new_arr = []
  for k in range(len(a)):
    ind = numbers.index(a[k])
    new_arr.append(arr[ind])
  return new_arr

# -----------------------------------------------------------------

def line_reg(header1):
  ima_pix2sec = float(header1['PIXSCALE_NEW'])
  nx = int(header1['NAXIS1'])
  ny = int(header1['NAXIS2'])
  scale = int(round(nx/8.*ima_pix2sec,-1))
  x2 = nx*9.8/10.
  x1 = x2 - scale/ima_pix2sec
  y1 = ny/7.
  y2 = y1
  
  return x1,y1,x2,y2,scale

# -----------------------------------------------------------------

def create_discrete_colormap(N=8):

    """
    Define new colormap for residuals
    :param N:
    :return:
    """

    # Define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000','#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = colors.ListedColormap(cpool[:N], 'i8')
    cm.register_cmap(cmap=cmap_i8)
    return cmap_i8

# -----------------------------------------------------------------

def define_scale_bar_length(x_extent,pix2sec):
    scale_bar = round((x_extent * pix2sec) / 6.,0)
    return int(5. * round(float(scale_bar)/5.)) # Length of the bar in arcsec
    
# -----------------------------------------------------------------

def crea_scale_bar(ax, x0, x1, y0, y1, pix2sec):

  offset_x_factor = 0.98
  offset_y_factor = 0.1 
  x_extent = x1 - x0

  scale_bar_length = define_scale_bar_length(x_extent, pix2sec) / 2. #### divide by 2 !!!

  xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
  yc = fabs(y0) + (y1-y0)* offset_y_factor
  ax.errorbar(xc, yc, xerr=scale_bar_length/pix2sec,color='black',capsize=1,c='black')
  ax.text(xc, yc, str(int(scale_bar_length*2.))+'\"', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')

# -----------------------------------------------------------------

def standard_setup(sp):

    """
    This function ...
    :param sp:
    :return:
    """

    sp.set_frame_color('black')
    sp.set_tick_labels_font(size='10')
    sp.set_axis_labels_font(size='12')
    #sp.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')
    sp.set_xaxis_coord_type('scalar')
    sp.set_yaxis_coord_type('scalar')
    sp.set_tick_color('black')
    sp.recenter(x=0.0, y=0.0,width=3.,height=0.6)
    sp.set_tick_xspacing(0.4)
    sp.set_tick_yspacing(0.25)
    sp.set_system_latex(True)
    sp.tick_labels.hide()
    sp.axis_labels.hide()

# -----------------------------------------------------------------
