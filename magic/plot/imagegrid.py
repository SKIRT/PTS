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
import matplotlib.gridspec as gridspec
import glob
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import LogNorm
import pyfits
from collections import OrderedDict
from textwrap import wrap

from astropy.io import fits
from pyfits import PrimaryHDU, Header
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import aplpy
import wcsaxes
import matplotlib.colors as mpl_colors
import matplotlib.colorbar as mpl_colorbar

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

# -----------------------------------------------------------------

class ImageGridPlotter(object):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        :param title:
        """

        # Set the title
        self.title = title

        # Figure and grid
        self._figure = None
        self._grid = None

        # Properties
        self.style = "dark" # "dark" or "light"
        self.transparent = True
        self.format = None

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

# -----------------------------------------------------------------

class StandardImageGridPlotter(ImageGridPlotter):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        :param title:
        """

        # Call the constructor of the base class
        super(StandardImageGridPlotter, self).__init__(title)

        # -- Attributes --

        # The images to be plotted
        self.images = OrderedDict()

        # Properties
        self.ncols = 7
        self.width = 16

        self.vmin = None

        self.colormap = "viridis"

    # -----------------------------------------------------------------

    def run(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def add_image(self, image, label):

        """
        This function ...
        :param image:
        :param label:
        :return:
        """

        self.images[label] = image

    # -----------------------------------------------------------------

    @property
    def nimages(self):

        """
        This function ...
        :return:
        """

        return len(self.images)

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine the necessary number of rows
        nrows = int(math.ceil(self.nimages / self.ncols))

        ratio = float(nrows) / float(self.ncols)
        height = ratio * self.width

        # Create the figure
        self._figure = plt.figure(figsize=(self.width, height))

        self._figure.subplots_adjust(hspace=0.0, wspace=0.0)

        #self._figure.text(0.385, 0.97, "Offset from centre (degrees)", color='black', size='16', weight='bold')
        #self._figure.text(0.02, 0.615, "Offset from centre (degrees)", color='black', size='16', weight='bold', rotation='vertical')

        def standard_setup(sp):
            sp.set_frame_color('black')
            sp.set_tick_labels_font(size='10')
            sp.set_axis_labels_font(size='12')
            # sp.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')
            sp.set_xaxis_coord_type('scalar')
            sp.set_yaxis_coord_type('scalar')
            sp.set_tick_color('black')
            sp.recenter(x=0.0, y=0.0, width=3., height=0.6)
            sp.set_tick_xspacing(0.4)
            sp.set_tick_yspacing(0.25)
            sp.set_system_latex(True)
            sp.tick_labels.hide()
            sp.axis_labels.hide()

        # Create grid
        #self._grid = AxesGrid(self._figure, 111,
        #                      nrows_ncols=(nrows, self.ncols),
        #                      axes_pad=0.0,
        #                      label_mode="L",
        #                      #share_all=True,
        #                      share_all=False,
        #                      cbar_location="right",
        #                      cbar_mode="single",
        #                      cbar_size="0.5%",
        #                      cbar_pad="0.5%")  # cbar_mode="single"

        gs = gridspec.GridSpec(nrows, self.ncols, wspace=0.0, hspace=0.0)

        # Loop over the images
        counter = 0
        ax = None
        for label in self.images:

            row = int(counter / self.ncols)
            col = counter % self.ncols

            frame = self.images[label]

            #ax = self._grid[counter]

            subplotspec = gs[row, col]

            #points = subplotspec.get_position(self._figure).get_points()
            #print(points)
            #x_min = points[0, 0]
            #x_max = points[1, 0]
            #y_min = points[0, 1]
            #y_max = points[1, 1]
            # width = x_max - x_min
            # height = y_max - y_min
            # ax = self._figure.add_axes([x_min, y_min, width, height])

            #ax = plt.subplot(subplotspec)
            #shareax = ax if ax is not None else None
            #ax = plt.subplot(subplotspec, projection=frame.wcs.to_astropy(), sharex=shareax, sharey=shareax)
            ax = plt.subplot(subplotspec, projection=frame.wcs.to_astropy())

            #lon = ax.coords[0]
            #lat = ax.coords[1]

            #overlay = ax.get_coords_overlay('fk5')
            #overlay.grid(color='white', linestyle='solid', alpha=0.5)

            # Determine the maximum value in the box and the mimimum value for plotting
            norm = ImageNormalize(stretch=LogStretch())
            #min_value = np.nanmin(frame)
            min_value = self.vmin if self.vmin is not None else np.nanmin(frame)
            max_value = 0.5 * (np.nanmax(frame) + min_value)

            #f1.show_colorscale(vmin=min_value, vmax=max_value, cmap="viridis")
            #f1.show_beam(major=0.01, minor=0.01, angle=0, fill=True, color='white')
            ## f1.axis_labels.show_y()
            #f1.tick_labels.set_xposition('top')
            #f1.tick_labels.show()

            ax.set_xticks([])
            ax.set_yticks([])
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])

            #ax.spines['bottom'].set_color("white")
            #ax.spines['top'].set_color("white")
            #ax.spines['left'].set_color("white")
            #ax.spines['right'].set_color("white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.tick_params(axis='x', colors="white")
            ax.tick_params(axis='y', colors="white")

            # Get the color map
            cmap = cm.get_cmap(self.colormap)

            # Set background color
            background_color = cmap(0.0)
            ax.set_axis_bgcolor(background_color)

            # Plot
            frame[np.isnan(frame)] = 0.0
            ax.imshow(frame, vmin=min_value, vmax=max_value, cmap=cmap, origin='lower', norm=norm, interpolation="nearest", aspect=1)

            # Add the label
            ax.text(0.95, 0.95, label, color='white', transform=ax.transAxes, fontsize=10, va="top", ha="right") # fontweight='bold'

            #ax.coords.grid(color='white')

            counter += 1

        all_axes = self._figure.get_axes()
        # show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
            #if ax.is_first_row():
            #    ax.spines['top'].set_visible(True)
            #if ax.is_last_row():
            #    ax.spines['bottom'].set_visible(True)
            #if ax.is_first_col():
            #    ax.spines['left'].set_visible(True)
            #if ax.is_last_col():
            #    ax.spines['right'].set_visible(True)

        # Add a colourbar

        #axisf3 = self._figure.add_axes(gs[row, col+1:])

        subplotspec = gs[row, col+1:]
        points = subplotspec.get_position(self._figure).get_points()
        #print("colorbar points:", points)

        x_min = points[0,0]
        x_max = points[1,0]
        y_min = points[0,1]
        y_max = points[1,1]

        #print((x_min, x_max), (y_min, y_max))

        #points_flattened = points.flatten()
        #print("colorbar:", points_flattened)

        x_center = 0.5 * (x_min + x_max)
        y_center = 0.5 * (y_min + y_max)

        width = 0.9* (x_max - x_min)
        height = 0.2 * (y_max - y_min)

        x_min = x_center - 0.5 * width
        x_max = x_center + 0.5 * width
        y_min = y_center - 0.5 * height
        y_max = y_center + 0.5 * height

        #ax_cm = plt.subplot(points)

        #ax_cm = plt.axes(points_flattened)

        ax_cm = self._figure.add_axes([x_min, y_min, width, height])

        cm_cm = cm.get_cmap(self.colormap)
        norm_cm = mpl_colors.Normalize(vmin=0, vmax=1)
        cb = mpl_colorbar.ColorbarBase(ax_cm, cmap=cm_cm, norm=norm_cm, orientation='horizontal')
        cb.set_label('Flux (arbitrary units)')



        # Set the title
        if self.title is not None: self._figure.suptitle("\n".join(wrap(self.title, 60)))

        #plt.tight_layout()

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the SED plot to a buffer ...")
        elif path is None: log.debug("Showing the SED plot ...")
        else: log.debug("Saving the SED plot to " + str(path) + " ...")

        if path is not None:
            # Save the figure
            plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------

class ResidualImageGridPlotter(ImageGridPlotter):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ResidualImageGridPlotter, self).__init__(title)

        # -- Attributes --

        # Set the title
        self.title = title

        # The rows of the grid
        self.rows = OrderedDict()
        self.plot_residuals = []

        # The names of the columns
        self.column_names = ["Observation", "Model", "Residual"]

        # Box (SkyRectangle) where to cut off the maps
        self.box = None

        self._plotted_rows = 0

    # -----------------------------------------------------------------

    def set_bounding_box(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        self.box = box

    # -----------------------------------------------------------------

    def add_row(self, image_a, image_b, label, residuals=True):

        """
        This function ...
        :param image_a:
        :param image_b:
        :param label:
        :param residuals:
        :return:
        """

        self.rows[label] = (image_a, image_b)
        if residuals: self.plot_residuals.append(label)

    # -----------------------------------------------------------------

    def set_column_names(self, name_a, name_b, name_residual="Residual"):

        """
        This function ...
        :param name_a:
        :param name_b:
        :param name_residual:
        :return:
        """

        self.column_names = [name_a, name_b, name_residual]

    # -----------------------------------------------------------------

    def run(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set default values for all attributes
        self.title = None
        self.rows = OrderedDict()
        self.plot_residuals = []
        self.column_names = ["Observation", "Model", "Residual"]
        self._figure = None
        self._grid = None
        self._plotted_rows = 0

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine the wcs with the smallest pixelscale
        reference_wcs = None
        for label in self.rows:
            if reference_wcs is None or reference_wcs.xy_average_pixelscale > self.rows[label][0].xy_average_pixelscale: reference_wcs = copy.deepcopy(self.rows[label][0].wcs)

        number_of_rows = len(self.rows)
        axisratio = float(self.rows[self.rows.keys()[0]][0].xsize) / float(self.rows[self.rows.keys()[0]][0].ysize)
        #print("axisratio", axisratio)

        one_frame_x_size = 3.
        fig_x_size = 3. * one_frame_x_size
        #fig_y_size = number_of_rows * one_frame_x_size / axisratio
        fig_y_size = one_frame_x_size * number_of_rows * 0.7

        # Create a figure
        self._figure = plt.figure(figsize=(fig_x_size, fig_y_size))
        self._figure.subplots_adjust(left=0.05, right=0.95)

        # Create grid
        self._grid = AxesGrid(self._figure, 111,
                                nrows_ncols=(len(self.rows), 3),
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
        for label in self.rows:
            wcs = self.rows[label][0].wcs

            rectangle = self.box.to_pixel(wcs)

            y_min = rectangle.lower_left.y
            y_max = rectangle.upper_right.y
            x_min = rectangle.lower_left.x
            x_max = rectangle.upper_right.x

            reference = self.rows[label][0][y_min:y_max, x_min:x_max]
            model = self.rows[label][1][y_min:y_max, x_min:x_max]
            data[label] = (reference, model)

            print(label, "box height/width ratio:", float(reference.shape[0])/float(reference.shape[1]))

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

            residual = (model - reference)/model

            # Plot the reference image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(reference, label, 0)

            # Plot the model image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(model, label, 1, min_int=vmin, max_int=vmax)

            # Plot the residual image
            x0, x1, y0, y1, vmin, vmax = self.plot_frame(residual, label, 2)

            self._plotted_rows += 3

        #self._grid.axes_llc.set_xlim(x0, x1)
        #self._grid.axes_llc.set_ylim(y0, y1)

        self._grid.axes_llc.set_xticklabels([])
        self._grid.axes_llc.set_yticklabels([])
        self._grid.axes_llc.get_xaxis().set_ticks([])  # To remove ticks
        self._grid.axes_llc.get_yaxis().set_ticks([])  # To remove ticks

        # Add title if requested
        #if self.title is not None: self._figure.suptitle(self.title, fontsize=12, fontweight='bold')

        plt.tight_layout()

        # Debugging
        log.debug("Saving the SED plot to " + path + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25, format=self.format, transparent=self.transparent)
        plt.close()

    # -----------------------------------------------------------------

    def plot_frame(self, frame, row_label, column_index, borders=(0,0,0,0), min_int=0., max_int=0.):

        """
        This function ...
        :param frame:
        :param column_index:
        :param row_label:
        :param borders:
        :return:
        """

        grid_index = self._plotted_rows + column_index

        x0 = borders[0]
        y0 = borders[1]

        #x1 = frame.xsize
        #y1 = frame.ysize
        x1 = frame.shape[1]
        y1 = frame.shape[0]

        #if column_index == 0:

        #    mod_data = fabs(data)
        #    min_data = np.min(mod_data[np.nonzero(mod_data)])
        #    for k in range(yaxis):
        #        for i in range(xaxis):
        #            if data[k, i] <= 0.: data[k, i] = min_data

        #if column_index == 2:

        #    for k in range(yaxis):
        #        for i in range(xaxis):
        #            if np.isnan(data[k, i]) == True:
        #                data[k, i] = 0.

        vmax = np.max(frame)  # np.mean([np.max(data_ski),np.max(data_ref)])
        vmin = np.min(frame)  # np.mean([np.min(data_ski),np.min(data_ref)])

        if min_int == 0.: min_int = vmin
        else: vmin = min_int

        if max_int == 0.: max_int = vmax
        else: vmax = max_int

        #aspect = "auto"
        aspect = "equal"
        if column_index != 2:
            im = self._grid[grid_index].imshow(frame, cmap='nipy_spectral_r', vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect=aspect)  # 'gist_ncar_r'
        else:
            im = self._grid[grid_index].imshow(frame, cmap=discrete_cmap(), vmin=0.001, vmax=1, interpolation='none', origin="lower", aspect=aspect)
            cb = self._grid[grid_index].cax.colorbar(im)

            # cb.set_xticklabels(labelsize=1)
            # grid[number+numb_of_grid].cax.toggle_label(True)
            for cax in self._grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                # cax.axis[cax.orientation].set_fontsize(3)
                cax.tick_params(labelsize=3)
                cax.set_ylim(0, 1)
                # cax.set_yticklabels([0, 0.5, 1])

        if column_index == 0:
            self._grid[grid_index].text(0.03, 0.95, row_label, color='black', transform=self._grid[grid_index].transAxes, fontsize=fsize + 2, fontweight='bold', va='top')

        # if numb_of_grid==0:
        #    crea_scale_bar(grid[number+numb_of_grid],x0,x1,y0,y1,pix2sec)
        #    crea_scale_bar(grid[number+numb_of_grid],x0,x1,y0,y1,pix2sec)

        return x0, x1, y0, y1, vmin, vmax

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
