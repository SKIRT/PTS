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
from scipy import ndimage
import copy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import colors
from matplotlib import cm
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty
from collections import defaultdict

# Import astronomical modules
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from ...core.tools import types
from ..basics.mask import MaskBase
from ..region.list import RegionList
from ..region.region import Region
from ..region.list import region_to_region_list
from ..core.image import Image
from ..core.frame import Frame
from ...core.tools import filesystem as fs
from ...core.tools.numbers import nan
from ...core.basics.containers import ordered_by_value
from ...core.basics.plot import MPLFigure, BokehFigure, mpl, bokeh, dark_pretty_colors, pretty_colors, filled_markers
from ..tools import plotting
from ...core.tools.utils import lazyproperty
from ..region.list import load_region_list
from ..core.mask import Mask
from ...core.tools import numbers

# -----------------------------------------------------------------

# EXISTS:

# astropy.visualization.wcsaxes.WCSAxesSubplot(fig, *args, **kwargs)[source] [edit on github]¶
# For making subplots with WCS

# A subclass class for WCSAxes
# fig is a matplotlib.figure.Figure instance.
# args is the tuple (numRows, numCols, plotNum), where the array of subplots in the figure has dimensions numRows, numCols, and where plotNum is the number of the subplot being created. plotNum starts at 1 in the upper left corner and increases to the right.
# If numRows <= numCols <= plotNum < 10, args can be the decimal integer numRows * 100 + numCols * 10 + plotNum.

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

        # Crop to
        self.crop_to = None
        self.cropping_factor = 1.

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

        # Get 'crop to'
        if "crop_to" in kwargs: self.crop_to = kwargs.pop("crop_to")
        if "cropping_factor" in kwargs: self.cropping_factor = kwargs.pop("cropping_factor")

    # -----------------------------------------------------------------

    def initialize_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Initializing the figure ...")

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

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
        #plt.tight_layout()

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

        # 2. Write
        if self.config.write: self.write()

        # 3. Make the plot
        self.plot()

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

        # Initialize the figure
        self.initialize_figure()

        # Setup the plots
        #self.plots = self.figure.create_grid(self.nrows, self.ncolumns)
        if self.config.coordinates: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white", projection=self.projection)
        else: self.plots, self.colorbar = self.figure.create_image_grid(self.nrows, self.ncolumns, return_colorbar=True, edgecolor="white")

        # Sort the frames on filter
        if self.config.sort_filters: self.sort()

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.frames[name].wcs

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

        # Downsampling required?
        copied = False
        if process and self.config.downsample and image.xsize > self.config.max_npixels or image.ysize > self.config.max_npixels:

            # Determine the downsample factor
            downsample_factor = max(image.xsize, image.ysize) / float(self.config.max_npixels)

            # Downsample
            if copy:
                image = image.copy()
                copied = True
            image.downsample(downsample_factor, convert=False, dilate_nans=False, dilate_infs=False)

        # Is copy still required?
        copy = copy and not copied

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

    def _process_frame(self, name, frame, copy=True):

        """
        This function ...
        :param name:
        :param frame:
        :param copy:
        :return:
        """

        # Initialize flags
        _cropped = False
        _downsampled = False
        _normalized = False

        # Crop the frame
        if self.crop_to is not None:

            # Debugging
            log.debug("Cropping the '" + name + "' frame ...")

            # Create cropped frame
            if copy: frame = frame.cropped_to(self.crop_to, factor=self.cropping_factor, out_of_bounds="expand")
            else: frame.crop_to(self.crop_to, factor=self.cropping_factor, out_of_bounds="expand")

            # Set flag
            _cropped = True

        # Downsampling required?
        if self.config.downsample and frame.xsize > self.config.max_npixels or frame.ysize > self.config.max_npixels:

            # Debugging
            log.debug("Downsampling the '" + name + "' frame ...")

            # Determine the downsample factor
            downsample_factor = max(frame.xsize, frame.ysize) / float(self.config.max_npixels)

            #print(frame.name, downsample_factor)

            # Downsample
            if _cropped or not copy: frame.downsample(downsample_factor, convert=False, dilate_nans=False, dilate_infs=False)
            else: frame = frame.downsampled(downsample_factor, convert=False, dilate_nans=False, dilate_infs=False)

            # Set flag
            _downsampled = True

        # Normalize the frame
        if self.config.normalize:

            # Debugging
            log.debug("Normalizing the '" + name + "' frame ...")

            # Create normalized frame
            if (_cropped or _downsampled) or not copy: frame.normalize()
            else: frame = frame.normalized()

            # Set flag
            _normalized = True

        # If no processing is performed, make a copy if necessary
        processed = _normalized or _cropped or _downsampled
        if not processed and copy: frame = frame.copy()

        # Return the new frame
        return frame

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
        if types.is_dictionary(mask): self.add_region_lists(name, regions)
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

    @property
    def colormap(self):

        """
        This function ...
        :return:
        """

        return cm.get_cmap(self.config.colormap)

    # -----------------------------------------------------------------

    @property
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

    def _plot_frame(self, name, row, col, vmin=None, vmax=None, return_image=False, return_normalization=False):

        """
        This function ...
        :param name:
        :param row:
        :param col:
        :param vmin:
        :param vmax:
        :param return_image:
        :param return_normalization:
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
        if self.config.background: plot.axes.set_axis_bgcolor(self.background_color)

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
        plot.axes.text(0.95, 0.95, name, color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")  # fontweight='bold'

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
        if self.config.background: plot.axes.set_axis_bgcolor(self.background_color)

        # Create NaNs image
        #nans = np.full((100,100), np.nan)

        #vmin_image, vmax_image = plotting.plot_box(nans, axes=plot.axes, interval=interval, scale=self.config.scale, cmap=self.colormap, alpha=self.config.alpha)

        #extent = None
        #aspect = "equal"
        #plot.axes.imshow(nans, origin="lower", vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, aspect=aspect, extent=extent)

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

            # Get row and column index
            row = int(index / self.ncolumns)
            col = index % self.ncolumns

            # Plot a frame
            if index < self.nframes:

                # Get name
                name = self.names[index]
                if name in plotted: continue

                # Plot the frame
                vmin_image, vmax_image, image, normalization = self._plot_frame(name, row, col, vmin=vmin, vmax=vmax, return_image=True, return_normalization=True)

                # Set vmin and vmax
                if self.config.share_scale:
                    vmin = vmin_image
                    vmax = vmax_image

            # Empty
            else: self._plot_empty(row, col)

        # Set colorbar
        if image is None: raise RuntimeError("No image is plotted")
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

        # Set the weighed flag
        self.weighed = None

        # Set the histogram flag
        self.histogram = None

        # The rows of the grid
        self.rows = OrderedDict()
        self.plot_residuals = []

        # Error maps (for weighed residuals)
        self.errors = dict()

        # The names of the columns
        if self.weighed: self.column_names = ["Observation", "Model", "Weighed residuals"]
        else: self.column_names = ["Observation", "Model", "Residuals"]

        # Box (SkyRectangle) where to cut off the maps
        self.box = None

        self._plotted_rows = 0

        self.absolute = False

    # -----------------------------------------------------------------

    def set_bounding_box(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        self.box = box

    # -----------------------------------------------------------------

    def add_row(self, frame_a, frame_b, name, residuals=True, errors=None):

        """
        This function ...
        :param frame_a:
        :param frame_b:
        :param name:
        :param residuals:
        :param errors:
        :return:
        """

        # Add the row
        self.rows[name] = (frame_a, frame_b)

        # If residuals have to be plotted
        if residuals:

            self.plot_residuals.append(name)
            if self.weighed:
                if errors is None: raise ValueError("Errors have to be specified to create weighed residuals")
                self.errors[name] = errors

    # -----------------------------------------------------------------

    def set_column_names(self, name_a, name_b, name_residual="Residual"):

        """
        This function ...
        :param name_a:
        :param name_b:
        :param name_residual:
        :return:
        """

        # Set names
        self.column_names = [name_a, name_b, name_residual]

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Write
        if self.config.write: self.write()

        # 3. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ResidualImageGridPlotter, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set default values for all attributes
        self.title = None
        self.rows = OrderedDict()
        self.errors = dict()
        self.plot_residuals = []
        self.column_names = ["Observation", "Model", "Residual"]
        self._plotted_rows = 0

    # -----------------------------------------------------------------

    @property
    def nresidual_rows(self):

        """
        This function ...
        :return:
        """

        return len(self.plot_residuals)

    # -----------------------------------------------------------------

    @property
    def has_residual_rows(self):

        """
        This function ...
        :return:
        """

        return self.nresidual_rows > 0

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        if not self.has_residual_rows: return 2
        else:
            if self.histogram: return 4
            else: return 3

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

    # -----------------------------------------------------------------

    def write_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the frames ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Determine the wcs with the smallest pixelscale
        reference_wcs = None
        for label in self.rows:
            if reference_wcs is None or reference_wcs.average_pixelscale > self.rows[label][0].average_pixelscale: reference_wcs = copy.deepcopy(self.rows[label][0].wcs)

        number_of_rows = len(self.rows)
        axisratio = float(self.rows[self.rows.keys()[0]][0].xsize) / float(self.rows[self.rows.keys()[0]][0].ysize)
        #print("axisratio", axisratio)

        one_frame_x_size = 3.
        fig_x_size = 3. * one_frame_x_size
        #fig_y_size = number_of_rows * one_frame_x_size / axisratio
        fig_y_size = one_frame_x_size * number_of_rows * 0.7

        # Create a figure
        #self._figure = plt.figure(figsize=(fig_x_size, fig_y_size))
        #self._figure.subplots_adjust(left=0.05, right=0.95)

        # Create grid
        self._grid = AxesGrid(self._figure, 111,
                                nrows_ncols=(len(self.rows), self.ncolumns),
                                axes_pad=0.02,
                                label_mode="L",
                                share_all=True,
                                cbar_location="right",
                                cbar_mode="single",
                                cbar_size="0.5%",
                                cbar_pad="0.5%",
                                )  # cbar_mode="single"

        for cax in self._grid.cbar_axes:
            cax.toggle_label(False)

        #rectangle_reference_wcs = self.box.to_pixel(reference_wcs)

        data = OrderedDict()

        greatest_shape = None

        if self.box is not None:

            for label in self.rows:
                wcs = self.rows[label][0].wcs

                rectangle = self.box.to_pixel(wcs)

                print(rectangle)
                print(rectangle.lower_left)
                print(rectangle.upper_right)

                y_min = rectangle.lower_left.y
                y_max = rectangle.upper_right.y
                x_min = rectangle.lower_left.x
                x_max = rectangle.upper_right.x

                reference = self.rows[label][0][y_min:y_max, x_min:x_max]
                model = self.rows[label][1][y_min:y_max, x_min:x_max]
                data[label] = (reference, model)

                #print(label, "box height/width ratio:", float(reference.shape[0])/float(reference.shape[1]))

                if greatest_shape is None or greatest_shape[0] < reference.shape[0]: greatest_shape = reference.shape

        else:

            for label in self.rows:
                reference = self.rows[label][0]
                model = self.rows[label][1]
                data[label] = (reference, model)
                if greatest_shape is None or greatest_shape[0] < reference.shape[0]: greatest_shape = reference.shape

        # Loop over the rows
        for label in self.rows:

            #wcs = self.rows[label][0].wcs

            if data[label][0].shape == greatest_shape:
                reference = data[label][0]
                model = data[label][1]
            else:
                factor = float(greatest_shape[0]) / float(data[label][0].shape[0])
                order = 0
                reference = ndimage.zoom(data[label][0], factor, order=order)
                model = ndimage.zoom(data[label][1], factor, order=order)

            # CREATE THE RESIDUAL
            if self.weighed:
                errors = self.errors[label]
                residual = (model - reference) / errors
            else:
                if self.absolute: residual = model - reference
                else: residual = (model - reference)/model

            # Plot the reference image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(reference, label, 0)

            # Plot the model image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(model, label, 1, vlimits=(vmin,vmax))

            # Plot the residual image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(residual, label, 2, vlimits=(vmin,vmax))

            self._plotted_rows += 3

        #self._grid.axes_llc.set_xlim(x0, x1)
        #self._grid.axes_llc.set_ylim(y0, y1)

        self._grid.axes_llc.set_xticklabels([])
        self._grid.axes_llc.set_yticklabels([])
        self._grid.axes_llc.get_xaxis().set_ticks([])  # To remove ticks
        self._grid.axes_llc.get_yaxis().set_ticks([])  # To remove ticks

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    def _plot_frame(self, name, row, col, vmin=None, vmax=None):

        """
        This function ...
        :param name:
        :param row:
        :param col:
        :param vmin:
        :param vmax:
        :return:
        """

        # Debugging
        log.debug("Adding the '" + name + "' frame to the plot ...")

    # -----------------------------------------------------------------

    def plot_frame(self, frame, row_label, column_index, borders=(0,0,0,0), vlimits=None):

        """
        This function ...
        :param frame:
        :param column_index:
        :param row_label:
        :param borders:
        :param vlimits:
        :return:
        """

        grid_index = self._plotted_rows + column_index

        x0 = borders[0]
        y0 = borders[1]

        #x1 = frame.xsize
        #y1 = frame.ysize
        x1 = frame.shape[1]
        y1 = frame.shape[0]

        #vmax = np.max(frame)  # np.mean([np.max(data_ski),np.max(data_ref)])
        #vmin = np.min(frame)  # np.mean([np.min(data_ski),np.min(data_ref)])
        #if min_int == 0.: min_int = vmin
        #else: vmin = min_int
        #if max_int == 0.: max_int = vmax
        #else: vmax = max_int

        if vlimits is None:
            min_value = self.vmin if self.vmin is not None else np.nanmin(frame)
            max_value = 0.5 * (np.nanmax(frame) + min_value)
        else:
            min_value = vlimits[0]
            max_value = vlimits[1]

        aspect = "equal"
        if column_index != 2:

            # Get the color map
            cmap = cm.get_cmap(self.colormap)

            # Set background color
            background_color = cmap(0.0)
            self._grid[grid_index].set_axis_bgcolor(background_color)

            # Plot
            frame[np.isnan(frame)] = 0.0
            norm = ImageNormalize(stretch=LogStretch())
            im = self._grid[grid_index].imshow(frame, cmap=cmap, vmin=min_value, vmax=max_value, interpolation="nearest", origin="lower", aspect=aspect, norm=norm) # 'nipy_spectral_r', 'gist_ncar_r'

        else:

            if self.absolute:
                # Get the color map
                cmap = cm.get_cmap(self.colormap)
                norm = ImageNormalize(stretch=LogStretch())
            else:
                cmap = discrete_cmap()
                min_value = 0.001
                max_value = 1.
                norm = None

            print(min_value, max_value)

            im = self._grid[grid_index].imshow(frame, cmap=cmap, vmin=min_value, vmax=max_value, interpolation="nearest", origin="lower", aspect=aspect, norm=norm)
            cb = self._grid[grid_index].cax.colorbar(im)

            # cb.set_xticklabels(labelsize=1)
            # grid[number+numb_of_grid].cax.toggle_label(True)
            for cax in self._grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                # cax.axis[cax.orientation].set_fontsize(3)
                cax.tick_params(labelsize=3)
                cax.set_ylim(min_value, max_value)
                # cax.set_yticklabels([0, 0.5, 1])

        if column_index == 0:
            self._grid[grid_index].text(0.03, 0.95, row_label, color='black', transform=self._grid[grid_index].transAxes, fontsize=fsize + 2, fontweight='bold', va='top')

        # if numb_of_grid==0:
        #    crea_scale_bar(grid[number+numb_of_grid],x0,x1,y0,y1,pix2sec)
        #    crea_scale_bar(grid[number+numb_of_grid],x0,x1,y0,y1,pix2sec)

        return x0, x1, y0, y1, min_value, max_value

# -----------------------------------------------------------------

fsize = 2

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

# Define new colormap for residuals
def discrete_cmap(N=8):
    # define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000','#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = colors.ListedColormap(cpool[0:N], 'i8')
    cm.register_cmap(cmap=cmap_i8)
    return cmap_i8


def define_scale_bar_length(x_extent,pix2sec):
    scale_bar = round((x_extent * pix2sec) / 6.,0)
    return int(5. * round(float(scale_bar)/5.)) # Length of the bar in arcsec
    

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
