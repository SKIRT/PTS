#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.colours Contains the ColourAnalyser class

# -----------------------------------------------------------------

# Import standard modules
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import glob
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import LogNorm
import pyfits
from collections import OrderedDict

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
        """

        # Set the title
        self.title = title

        # The rows of the grid
        self.rows = OrderedDict()
        self.plot_residuals = []

        # The names of the columns
        self.column_names = ["Observation", "Model", "Residual"]

        self._figure = None
        self._grid = None
        self._plotted_rows = 0

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

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

        pass

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        #hdulist = fits.open(ref_arr[0])
        #inframe1 = hdulist[0].data
        #header1 = hdulist[0].header
        #ySize, xSize = inframe1.shape

        #parser.add_argument("Scale", nargs='?', const=1., help="Input scale in [arcsec/pix]", type=float, default=1.)
        #parser.add_argument("borders", nargs='?', const='0,0,0,0', help="Input the borders of the frame to be plotted: x1,y1,x2,y2", type=str, default='0,0,0,0')

        # Create a figure
        self._figure = plt.figure(0, (6, 20))
        self._figure.subplots_adjust(left=0.05, right=0.95)

        # Creat grid
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

        #number = 0
        #min_int = min_Int

        # Loop over the rows
        for label in self.rows:

            # Plot the reference image
            x0, x1, y0, y1, vmin, vmax = self.plot_indiv_frame(0, label)

            # Plot the model image
            x0, x1, y0, y1, vmin, vmax = self.plot_indiv_frame(1, label, min_int=vmin, max_int=vmax)

            # Plot the residual image
            x0, x1, y0, y1, vmin, vmax = self.plot_indiv_frame(2, label)

            self._plotted_rows += 3

            #number = number + 3

        #min_int = min_Int

        self._grid.axes_llc.set_xlim(x0, x1)
        self._grid.axes_llc.set_ylim(y0, y1)
        self._grid.axes_llc.set_xticklabels([])
        self._grid.axes_llc.set_yticklabels([])
        self._grid.axes_llc.get_xaxis().set_ticks([])  # To remove ticks
        self._grid.axes_llc.get_yaxis().set_ticks([])  # To remove ticks


        # Add title if requested
        if self.title is not None: self._figure.suptitle(self.title, fontsize=14, fontweight='bold')

        # Debugging
        log.debug("Saving the SED plot to " + path + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    def plot_indiv_frame(self, column_index, row_label, borders=(0,0,0,0), min_int=0., max_int=0.):

        """
        This function ...
        :param column_index:
        :param row_label:
        :param borders:
        :return:
        """

        grid_index = self._plotted_rows + column_index

        if column_index == 2:
            frame = self.rows[row_label][1] - self.rows[row_label][0]
        else:
            # Get the image frame to be plotted
            frame = self.rows[row_label][column_index]

        #borders = borders.split(',')
        x0 = borders[0]
        y0 = borders[1]
        #x1 = borders[2]
        #y1 = borders[3]

        x1 = frame.xsize
        y1 = frame.ysize

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

        if column_index != 2:
            im = self._grid[grid_index].imshow(frame, cmap='nipy_spectral_r', vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')  # 'gist_ncar_r'
        else:
            im = self._grid[grid_index].imshow(frame, cmap=discrete_cmap(), vmin=0.001, vmax=1, interpolation='none', origin="lower", aspect='equal')
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

# EXAMPLE: python ~/CurrentWork/ImaPrep/IMAN/SKIRT/plot_fitskirt_images.py imf_output.fski /
# /home/amosenko/CurrentWork/HEROES/models/IC2531/reference /home/amosenko/CurrentWork/HEROES/models/IC2531/RESULTS/fit3


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
    
    
    
def sort_bands(arr):
  new_arr = arr[:]
  out_arr = arr[:]
  bands = arr[:]
  bands_out = arr[:]
  #print arr
  for k in range(len(arr)): 
    band = str(arr[k].split('/')[-1].split('_')[1])
    telescope = str(arr[k].split('/')[-1].split('_')[0])
    #print telescope
    if band=='U':
      numb = 1
    elif band=='u':
      numb = 2
    elif band=='B':
      numb = 3    
    elif band=='g':
      numb = 4
    elif band=='V':
      numb = 5
    elif band=='R':
      numb = 6
    elif band=='r':
      numb = 7    
    elif band=='I':
      numb = 8
    elif band=='i':
      numb = 9
    elif band=='z':
      numb = 10
    elif band=='J':
      #print telescope,arr[k]
      numb = 11
      #print numb
    elif band=='H':
      numb = 12
    elif band=='K':
      if telescope=='2MASS':
	numb = 13
      else:
	numb = 14
    elif band=='w1' or band=='1':
      if telescope=='WISE':
	numb = 15
      else:
	numb = 16
    elif band=='w2' or band=='2':
      if telescope=='WISE':
	numb = 18
      else:
	numb = 17
    elif band=='w3':
      numb = 19
    elif band=='w4':
      numb = 20
    new_arr[k]=numb
    bands[k]=band
  a = sorted(new_arr)
  #print a
  for k in range(len(a)):
    number = new_arr.index(a[k])
    out_arr[k] = arr[number]
    bands_out[k] = out_arr[k].split('/')[-1].split('_')[1]
    if (bands_out[k]=='w1' or bands_out[k]=='w2') and out_arr[k].split('/')[-1].split('_')[0]=='Spitzer':
      if bands_out[k]=='w1':
	bands_out[k]='3.6'
      if bands_out[k]=='w2':
	bands_out[k]='I2'
    if bands_out[k]=='1':
      bands_out[k]='I1'
  return out_arr,bands_out

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
    

def crea_scale_bar(ax,x0,x1,y0,y1,pix2sec):
  offset_x_factor = 0.98
  offset_y_factor = 0.1 
  x_extent = x1 - x0

  scale_bar_length = define_scale_bar_length(x_extent,pix2sec) / 2. #### divide by 2 !!!

  xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
  yc = fabs(y0) + (y1-y0)* offset_y_factor
  ax.errorbar(xc, yc, xerr=scale_bar_length/pix2sec,color='black',capsize=1,c='black')
  ax.text(xc, yc, str(int(scale_bar_length*2.))+'\"', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
    




def main(fski_file,ref_dir,res_dir,min_Int,pix2sec,borders='0,0,0,0',format_out='pdf'):

    fski_file = fski_file.split('.fski')[0]

    #1. Original images should be in ref_dir
    old_arr = glob.glob(ref_dir+'/*_norm.fits')  

    ref_arr,bands = sort_bands(old_arr)  


    #2. Results should be placed in res_dir
    old_arr = glob.glob(res_dir+'/' + str(fski_file) + '_Best_*.fits')
    mod_arr = sort_numbs(old_arr) 

    
    old_arr = glob.glob(res_dir+'/' + str(fski_file) + '_Residual_*.fits')
    res_arr = sort_numbs(old_arr)


    hdulist = pyfits.open(ref_arr[0])
    inframe1 = hdulist[0].data
    header1 = hdulist[0].header
    ySize, xSize = inframe1.shape

    fig = plt.figure(0, (6, 20))
    fig.subplots_adjust(left=0.05, right=0.95) 

    grid = AxesGrid(fig, 111,
			nrows_ncols=(len(res_arr), 3),
			axes_pad=0.02,
			label_mode="L",
			share_all=True,
			cbar_location="right",
			cbar_mode="single",
			cbar_size="0.5%",
			cbar_pad="0.5%",
			) # cbar_mode="single"

    for cax in grid.cbar_axes:
	cax.toggle_label(False)

    number = 0
    min_int = min_Int
    for i in range(len(res_arr)): 
        # Reference image:
        x0,x1,y0,y1,vmin,vmax = plot_indiv_frame(grid,number,ref_arr[i],'Data',pix2sec=pix2sec,borders=borders,min_int=min_int,max_int = 0., name = str(bands[i]))

        # Model image:
        x0,x1,y0,y1,vmin,vmax = plot_indiv_frame(grid,number,mod_arr[i],'Model',pix2sec=pix2sec,borders=borders,min_int=vmin,max_int = vmax,name = str(bands[i]))        

        # Residual image:
        x0,x1,y0,y1,vmin,vmax = plot_indiv_frame(grid,number,res_arr[i],'Residual',pix2sec=pix2sec,borders=borders,min_int=min_int,max_int = 0.,name = str(bands[i]))
	number = number + 3
	min_int = min_Int
    grid.axes_llc.set_xlim(x0,x1)
    grid.axes_llc.set_ylim(y0,y1)
    grid.axes_llc.set_xticklabels([])
    grid.axes_llc.set_yticklabels([])
    grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
    grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks


    plt.draw()  
    plt.savefig('fitskirt_models.%s' % (format_out), bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()   
    #plt.show()
    return '2d_decom_res.%s' % (format_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plotting FitSKIRT model images, residuals and reference images")
    parser.add_argument("fski_file", help="Input .fski file")
    parser.add_argument("ref_dir", nargs='?', const=1, help="Input directory with the reference images",type=str,default='./reference')
    parser.add_argument("res_dir", nargs='?', const=1, help="Input directory with the FitSKIRT models",type=str,default='./models')
    parser.add_argument("min_int", nargs='?', const=0., help="Input minimum DN level to highlight the structure in [DN]",type=float,default=0.0000000001) 
    parser.add_argument("Scale", nargs='?', const=1., help="Input scale in [arcsec/pix]",type=float,default=1.) 
    parser.add_argument("borders", nargs='?', const='0,0,0,0', help="Input the borders of the frame to be plotted: x1,y1,x2,y2",type=str,default='0,0,0,0')    
    parser.add_argument("format_out", nargs='?', const=0., help="Input the format of the output file",type=str,default='pdf')     
    args = parser.parse_args()
  
    fski_file = args.fski_file
    ref_dir = args.ref_dir
    res_dir = args.res_dir
    min_int = float(args.min_int)
    pix2sec = float(args.Scale)
    borders = str(args.borders)
    format_out = str(args.format_out)
    
    main(fski_file,ref_dir,res_dir,min_int,pix2sec,borders='0,0,0,0',format_out='pdf')