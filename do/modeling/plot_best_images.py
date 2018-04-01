# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package

# Import standard modules

import aplpy
import scipy.constants as cst
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.tools.plotting import get_vmin_vmax
from pts.core.filter.filter import parse_filter
from pts.core.tools import filesystem as fs
from pts.magic.tools.headers import get_header, get_filter
from pts.magic.core.frame import Frame
from pts.core.basics.log import log

# ------------------------------------------------------------------------------

# Set fundamental physiscal constasnts
c = cst.c # m/s
h = cst.h # J s
k = cst.k # J/K
Mpc = 3.08567758e22 # m
Msun = 1.9891e30 # kg
Lsun = 3.828e26  # W

# ------------------------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# Get the galaxy center
center = environment.galaxy_center
ra = center.ra.to("deg").value
dec = center.dec.to("deg").value

# ------------------------------------------------------------------------------

light_theme = "light"
dark_theme = "dark"
themes = [light_theme, dark_theme]

# ------------------------------------------------------------------------------

default_filters = [parse_filter("FUV"), parse_filter("SDSS r"), parse_filter("I1"), parse_filter("MIPS 24mu"), parse_filter("Pacs red")]

# ------------------------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation name
definition.add_required("generation", "string", "generation name")

# Filters
definition.add_positional_optional("filters", "filter_list", "filters for which to plot the images", default_filters)

# Colormap
definition.add_optional("cmap", "string", "colormap", "viridis")

# Options
definition.add_optional("axes_label_size", "positive_integer", "axes label size", 14)
definition.add_optional("ticks_label_size", "positive_integer", "ticks label size", 8)
definition.add_optional("legend_fontsize", "positive_integer", "legend fontsize", 14)
definition.add_optional("legend_markers_cale", "positive_integer", "legend marker scale", 0)
definition.add_optional("lines_marker_size", "positive_real", "lines marker size", 2.5)

# Dark or light theme?
definition.add_optional("theme", "string", "theme for the plot", light_theme, choices=themes)

# Pixelscale
definition.add_optional("pixelscale", "angle", "pixelscale for the images", "5 arcsec", convert_default=True)
definition.add_optional("radius", "angle", "radius of the field of view", "0.1 deg", convert_default=True)

# Interval
definition.add_optional("interval", "string", "interval", "pts")

definition.add_optional("filter_fontsize", "positive_real", "fontsize for filter titles", 18)

# ------------------------------------------------------------------------------

# Get configuration
config = parse_arguments("plot_best_images", definition, "plot images of best simulation")

# ------------------------------------------------------------------------------

# Load the generation
run = runs.load(config.run)
generation = run.get_generation(config.generation)

# Get the chi squared table
chi_squared = generation.chi_squared_table

# Get the name of the best simulation
simulation_name = chi_squared.best_simulation_name

# ------------------------------------------------------------------------------

#plt.rcParams.update({'font.size':20})
plt.rcParams["axes.labelsize"]  = config.axes_label_size #16 #default 20
plt.rcParams["xtick.labelsize"] = config.ticks_label_size #10 #default 16
plt.rcParams["ytick.labelsize"] = config.ticks_label_size #10 #default 16
plt.rcParams["legend.fontsize"] = config.legend_fontsize #10 #default 14
plt.rcParams["legend.markerscale"] = config.legend_markers_cale
plt.rcParams["lines.markersize"] = config.lines_marker_size #4 #default 4
plt.rcParams["axes.linewidth"]  = config.linewidth

# Set light theme
if config.theme == light_theme:

    text_color = "black"
    frame_color = "black"
    background_color = "white"

    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['xtick.color'] = 'black'
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams["axes.labelcolor"] = 'black'
    plt.rcParams["text.color"] = 'black'
    # plt.rcParams["axes.titlecolor"]='black'

# Dark theme
elif config.theme == dark_theme:

    text_color = "white"
    frame_color = "white"
    background_color = "black"

    plt.rcParams['axes.facecolor']='black'
    plt.rcParams['savefig.facecolor']='black'
    plt.rcParams['axes.edgecolor']='white'
    plt.rcParams['xtick.color']='white'
    plt.rcParams['ytick.color']='white'
    plt.rcParams["axes.labelcolor"]='white'
    plt.rcParams["text.color"]='white'
    #plt.rcParams["axes.titlecolor"]='white'

else: raise ValueError("Invalid value for 'theme': " + config.theme )

#plt.rcParams['xtick.major.size'] = 5
#plt.rcParams['xtick.major.width'] = 2
#plt.rcParams['ytick.major.size'] = 5
#plt.rcParams['ytick.major.width'] = 2

# ------------------------------------------------------------------------------

# Get the mock images
image_paths = generation.get_image_paths_for_simulation(simulation_name)

# Create a dictionary for the images
#mock_images = OrderedDict()

nfilters = len(config.filters)

# Create a list for the images
mock_images = [None] * nfilters

# Loop over the image paths
for filepath in image_paths:

    # Get image name
    name = fs.strip_extension(fs.name(filepath))

    # Get header
    header = get_header(filepath)

    # Get the filter
    fltr = get_filter(name, header=header)

    # Check whether the filter is in the list of filters to be plotted
    if fltr not in config.filters: continue

    # Get the index for this filter
    index = config.filters.index(fltr)

    # Load the image
    frame = Frame.from_file(filepath)

    # Replace zeroes and negatives
    frame.replace_zeroes_by_nans()
    frame.replace_negatives_by_nans()

    # Set the image
    mock_images[index] = frame

# ------------------------------------------------------------------------------

# Get the observed images
observed_images = []
for fltr in config.filters:

    # Check
    if not environment.has_photometry_image_for_filter(fltr):
        observed_images.append(None)
        continue

    # Get the photometry image
    filepath = environment.photometry_image_paths_for_filters[fltr]

    # Load the image
    frame = Frame.from_file(filepath)

    # Replace zeroes and negatives
    frame.replace_zeroes_by_nans()
    frame.replace_negatives_by_nans()

    # Set the image
    observed_images.append(frame)

# ------------------------------------------------------------------------------

# Create residuals
residuals = [None] * nfilters

for index in range(nfilters):

    fltr = config.filters[index]

    mock_image = mock_images[index]
    observed_image = observed_images[index]

    if mock_image is None and observed_image is None:
        raise ValueError("No data for the " + str(fltr) + " filter")

    elif mock_image is None:
        log.warning("No mock image for the " + str(fltr) + " filter")
        continue

    elif observed_image is None:
        log.warning("No observed image for the " + str(fltr) + " filter")
        continue

    # Create the residual frame
    res = Frame((observed_image - mock_image) / observed_image, wcs=observed_image.wcs)

    # Add the residuals frame
    residuals[index] = res

# ------------------------------------------------------------------------------

# Create the figure
fig = plt.figure(figsize=(15, 10))

# ------------------------------------------------------------------------------

pixelscale_deg = config.pixelscale.to("deg").value
radius_deg = config.radius.to("deg").value

# ------------------------------------------------------------------------------

#f1 = aplpy.FITSFigure(path_to_nan+'GALEX_FUV_nan.fits', slices=[0], figure=fig, subplot=[0.1,0.65,0.168,0.3])
f1 = aplpy.FITSFigure(observed_images[0].to_hdu(), figure=fig, subplot=[0.1,0.65,0.168,0.3])

f1_vmin, f1_vmax = get_vmin_vmax(observed_images[0].data, interval=config.interval)
f1.show_colorscale(vmin=f1_vmin, vmax=f1_vmax, cmap=config.cmap)

f1.recenter(ra, dec, radius=radius_deg)
f1.ticks.set_xspacing(pixelscale_deg)

f1.frame.set_color(frame_color)

f1._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f1._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)
f1._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f1._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f1.set_nan_color('#440154')

f1._ax1.scatter(ra,dec, marker='.', label='Observation')

legend1 = f1._ax1.legend(loc='upper right', fontsize=12, fancybox=True, framealpha=0, numpoints=None)
plt.setp(legend1.get_texts(), color=text_color)

# Set title
f1._ax1.set_title(str(config.filters[0]), fontsize=config.filter_fontsize)

# ------------------------------------------------------------------------------

#f2 = aplpy.FITSFigure(path_to_nan+'R_nan.fits', slices=[0], figure=fig, subplot=[0.268,0.65,0.168,0.3])
f2 = aplpy.FITSFigure(observed_images[1].to_hdu(), figure=fig, subplot=[0.268,0.65,0.168,0.3])

# Set colorscale
f2_vmin, f2_vmax = get_vmin_vmax(observed_images[1].data, interval=config.interval)
f2.show_colorscale(vmin=f2_vmin, vmax=f2_vmax, cmap=config.cmap)

f2.recenter(ra, dec, radius=radius_deg)
f2.ticks.set_xspacing(pixelscale_deg)

f2.frame.set_color(frame_color)

f2._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f2._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f2.set_nan_color('#440154')

# Set title
f2._ax1.set_title(str(config.filters[1]), fontsize=config.filter_fontsize)

# ------------------------------------------------------------------------------

#f3 = aplpy.FITSFigure(path_to_nan+'IRAC_I1_nan.fits', slices=[0], figure=fig, subplot=[0.436,0.65,0.168,0.3])
f3 = aplpy.FITSFigure(observed_images[2].to_hdu(), figure=fig, subplot=[0.436,0.65,0.168,0.3])

# Set colorscale
f3_vmin, f3_vmax = get_vmin_vmax(observed_images[2].data, interval=config.interval)
f3.show_colorscale(vmin=f3_vmin, vmax=f3_vmax, cmap=config.cmap)

f3.recenter(ra, dec, radius=radius_deg)
f3.ticks.set_xspacing(pixelscale_deg)

f3.frame.set_color(frame_color)

f3._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f3._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f3.set_nan_color('#440154')

# Set title
f3._ax1.set_title(str(config.filters[2]), fontsize=config.filter_fontsize)

# ------------------------------------------------------------------------------

#f4 = aplpy.FITSFigure(path_to_nan+'MIPS_24mu_nan.fits', slices=[0], figure=fig, subplot=[0.604,0.65,0.168,0.3])
f4 = aplpy.FITSFigure(observed_images[3].to_hdu(), figure=fig, subplot=[0.604,0.65,0.168,0.3])

# Set colorscale
f4_vmin, f4_vmax = get_vmin_vmax(observed_images[3].data, interval=config.interval)
f4.show_colorscale(vmin=f4_vmin, vmax=f4_vmax, cmap=config.cmap)

f4.recenter(ra, dec, radius=radius_deg)
f4.ticks.set_xspacing(pixelscale_deg)

f4.frame.set_color(frame_color)

f4._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f4._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f4.set_nan_color('#440154')

# Set title
f4._ax1.set_title(str(config.filters[3]), fontsize=config.filter_fontsize)

# ------------------------------------------------------------------------------

#f5 = aplpy.FITSFigure(path_to_nan+'PACS_red_nan.fits', slices=[0], figure=fig, subplot=[0.772,0.65,0.168,0.3])
f5 = aplpy.FITSFigure(observed_images[4].to_hdu(), figure=fig, subplot=[0.772,0.65,0.168,0.3])

# Set colorscale
f5_vmin, f5_vmax = get_vmin_vmax(observed_images[4].data, interval=config.interval)
f5.show_colorscale(vmin=f5_vmin, vmax=f5_vmax, cmap=config.cmap)

f5.recenter(ra, dec, radius=radius_deg)
f5.ticks.set_xspacing(pixelscale_deg)

f5.frame.set_color(frame_color)

f5._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f5._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f5.set_nan_color('#440154')

# Set title
f5._ax1.set_title(str(config.filters[4]), fontsize=config.filter_fontsize)

# ------------------------------------------------------------------------------

#f6 = aplpy.FITSFigure(path_to_nan+'GALEX FUV_nan.fits', figure=fig, subplot=[0.1,0.35,0.168,0.3])
f6 = aplpy.FITSFigure(mock_images[0].to_hdu(), figure=fig, subplot=[0.1,0.35,0.168,0.3])

# Set colorscale
f6_vmin, f6_vmax = get_vmin_vmax(mock_images[0].data, interval=config.interval)
f6.show_colorscale(vmin=f6_vmin, vmax=f6_vmax, cmap=config.cmap)

f6.tick_labels.set_font(size='small')
f6.recenter(ra, dec, radius=radius_deg)
f6.ticks.set_xspacing(pixelscale_deg)

f6.frame.set_color(frame_color)

f6._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f6._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)
f6._ax1.scatter(ra,dec, marker='.', label='Model')

legend6 = f6._ax1.legend(loc='upper right', fontsize=12, fancybox=False, framealpha=0, numpoints=None)
plt.setp(legend6.get_texts(), color='w')

f6.set_nan_color('#440154')

# ------------------------------------------------------------------------------

#f7 = aplpy.FITSFigure(path_to_nan+'SDSS r_nan.fits', figure=fig, subplot=[0.268,0.35,0.168,0.3])
f7 = aplpy.FITSFigure(mock_images[1].to_hdu(), figure=fig, subplot=[0.268,0.35,0.168,0.3])

# Get colorscale
f7_vmin, f7_vmax = get_vmin_vmax(mock_images[1].data, interval=config.interval)
f7.show_colorscale(vmin=f7_vmin, vmax=f7_vmax, cmap=config.cmap)

f7.recenter(ra, dec, radius=radius_deg)
f7.ticks.set_xspacing(pixelscale_deg)

f7.frame.set_color(frame_color)

f7._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f7._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f7.set_nan_color('#440154')

# ------------------------------------------------------------------------------

#f8 = aplpy.FITSFigure(path_to_nan+'IRAC I1_nan.fits', figure=fig, subplot=[0.436,0.35,0.168,0.3])
f8 = aplpy.FITSFigure(mock_images[2].to_hdu(), figure=fig, subplot=[0.436,0.35,0.168,0.3])

# Set colorscale
f8_vmin, f8_vmax = get_vmin_vmax(mock_images[2].data, interval=config.interval)
f8.show_colorscale(vmin=f8_vmin, vmax=f8_vmax, cmap=config.cmap)

f8.recenter(ra, dec, radius=radius_deg)
f8.ticks.set_xspacing(pixelscale_deg)

f8.frame.set_color(frame_color)

f8._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f8._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f8.set_nan_color('#440154')

# ------------------------------------------------------------------------------

#f9 = aplpy.FITSFigure(path_to_nan+'MIPS 24mu_nan.fits', figure=fig, subplot=[0.604,0.35,0.168,0.3])
f9 = aplpy.FITSFigure(mock_images[3].to_hdu(), figure=fig, subplot=[0.604,0.35,0.168,0.3])

# Set colorscale
f9_vmin, f9_vmax = get_vmin_vmax(mock_images[3].data, interval=config.interval)
f9.show_colorscale(vmin=f9_vmin, vmax=f9_vmax, cmap=config.cmap)

f9.recenter(ra, dec, radius=radius_deg)
f9.ticks.set_xspacing(pixelscale_deg)

f9.frame.set_color(frame_color)

f9._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f9._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f9.set_nan_color('#440154')

# ------------------------------------------------------------------------------

# Hidden image (this should go fisrt then f10)
#f16 = aplpy.FITSFigure(path_to_nan+'Pacs red_nan.fits', figure=fig, subplot=[0.772,0.35,0.168,0.3])
f16 = aplpy.FITSFigure(mock_images[4].to_hdu(), figure=fig, subplot=[0.772,0.35,0.168,0.3])

# Set colorscale
f16_vmin, f16_vmax = get_vmin_vmax(mock_images[4].data, interval=config.interval)
f16.show_colorscale(vmin=f16_vmin, vmax=f16_vmax, cmap=config.cmap)

f16.recenter(ra, dec, radius=radius_deg)

f16.ticks.hide_x()
f16.ticks.hide_y()

f16.set_nan_color('#440154')

f16.add_colorbar()

f16.colorbar.set_box([0.941, 0.35, 0.01, 0.6], box_orientation='vertical')
f16.colorbar.set_axis_label_text('Flux (Arbitrary Units)')

# ------------------------------------------------------------------------------

#f10 = aplpy.FITSFigure(path_to_nan+'Pacs red_nan.fits', figure=fig, subplot=[0.772,0.35,0.168,0.3])
f10 = aplpy.FITSFigure(mock_images[4].to_hdu(), figure=fig, subplot=[0.772,0.35,0.168,0.3])

# Set colorscale
f10_vmin, f10_vmax = get_vmin_vmax(mock_images[4].data, interval=config.interval)
f10.show_colorscale(vmin=f10_vmin, vmax=f10_vmax, cmap=config.cmap)

f10.recenter(ra, dec, radius=radius_deg)
f10.ticks.set_xspacing(pixelscale_deg)

f10.frame.set_color(frame_color)

f10._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f10._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f10.set_nan_color('#440154')

# ------------------------------------------------------------------------------

#f11 = aplpy.FITSFigure(path_to_res+'GALEX_FUV_res.fits', slices=[0], figure=fig, subplot=[0.1,0.05,0.168,0.3])
f11 = aplpy.FITSFigure(residuals[0].to_hdu(), figure=fig, subplot=[0.1,0.05,0.168,0.3])

# Set colorscale
f11.show_colorscale(vmin=-1.0, vmax=1.0, cmap='RdBu')

# Set ...
f11.recenter(ra, dec, radius=radius_deg)
f11.ticks.set_xspacing(pixelscale_deg)

f11.frame.set_color(frame_color)

f11._ax1.scatter(ra,dec, marker='.', label='Relative \nResidual')
f11._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f11._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

# Set legend
legend11 = f11._ax1.legend(loc='lower right', fontsize=12, fancybox=False, framealpha=0, numpoints=None)
plt.setp(legend11.get_texts(), color='w')

f11.set_nan_color(background_color)

# ------------------------------------------------------------------------------

#f12 = aplpy.FITSFigure(path_to_res+'R_res.fits', slices=[0], figure=fig, subplot=[0.268,0.05,0.168,0.3])
f12 = aplpy.FITSFigure(residuals[1].to_hdu(), figure=fig, subplot=[0.268,0.05,0.168,0.3])

# Set colorscale
f12.show_colorscale(vmin=-1.0, vmax=1.0,cmap='RdBu')

# Set ...
f12.recenter(ra, dec, radius=radius_deg)
f12.ticks.set_xspacing(pixelscale_deg)

f12.frame.set_color(frame_color)

f12._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f12._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f12.set_nan_color(background_color)

# ------------------------------------------------------------------------------

#f13 = aplpy.FITSFigure(path_to_res+'IRAC_I1_res.fits', slices=[0], figure=fig, subplot=[0.436,0.05,0.168,0.3])
f13 = aplpy.FITSFigure(residuals[2].to_hdu(), figure=fig, subplot=[0.436,0.05,0.168,0.3])

f13.tick_labels.set_font(size='small')

# Set colorscale
f13.show_colorscale(vmin=-1.0, vmax=1.0, cmap='RdBu')

# Set ...
f13.recenter(ra, dec, radius=radius_deg)
f13.ticks.set_xspacing(pixelscale_deg)

f13.frame.set_color(frame_color)

f13._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f13._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f13.set_nan_color(background_color)

# ------------------------------------------------------------------------------

#f14 = aplpy.FITSFigure(path_to_res+'MIPS_24mu_res.fits', slices=[0], figure=fig, subplot=[0.604,0.05,0.168,0.3])
f14 = aplpy.FITSFigure(residuals[3].to_hdu(), figure=fig, subplot=[0.604,0.05,0.168,0.3])

# Set colorscale
f14.show_colorscale(vmin=-1.0, vmax=1.0,cmap='RdBu')

f14.recenter(ra, dec, radius=radius_deg)
f14.ticks.set_xspacing(pixelscale_deg)

f14.frame.set_color(frame_color)

f14._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f14._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f14.set_nan_color(background_color)

# ------------------------------------------------------------------------------

# Hidden figure (this should go fisrt then f15)
#f17 = aplpy.FITSFigure(path_to_res+'PACS_red_res.fits', slices=[0], figure=fig, subplot=[0.772,0.05,0.168,0.1])
f17 = aplpy.FITSFigure(residuals[4].to_hdu(), figure=fig, subplot=[0.772,0.05,0.168,0.1])

# Set colorscale
f17.show_colorscale(vmin=-1, vmax=1, cmap='RdBu')

f17.recenter(ra, dec, radius=radius_deg)

f17.ticks.hide_x()
f17.ticks.hide_y()

f17.set_nan_color(background_color)

f17.add_colorbar()

f17.colorbar.set_box([0.941, 0.05, 0.01, 0.299],box_orientation='vertical')

# ------------------------------------------------------------------------------

#f15 = aplpy.FITSFigure(path_to_res+'PACS_red_res.fits', slices=[0], figure=fig, subplot=[0.772,0.05,0.168,0.3])
f15 = aplpy.FITSFigure(residuals[4].to_hdu(), slices=[0], figure=fig, subplot=[0.772,0.05,0.168,0.3])

# Set colorscale
f15.show_colorscale(vmin=-1.0, vmax=1.0, cmap='RdBu')

f15.recenter(ra, dec, radius=radius_deg)
f15.ticks.set_xspacing(pixelscale_deg)

f15.frame.set_color(frame_color)

f15._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f15._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)
f15._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
f15._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

f15.set_nan_color(background_color)

# ------------------------------------------------------------------------------

# Hide xaxis & xtick labels with aplpy for the appropriate subplots
# First row
f1.axis_labels.hide_x()
f1.tick_labels.hide_x()

f2.axis_labels.hide_x()
f2.tick_labels.hide_x()

f3.axis_labels.hide_x()
f3.tick_labels.hide_x()

f4.axis_labels.hide_x()
f4.tick_labels.hide_x()

f5.axis_labels.hide_x()
f5.tick_labels.hide_x()

# Second row
f6.axis_labels.hide_x()
f6.tick_labels.hide_x()

f7.axis_labels.hide_x()
f7.tick_labels.hide_x()

f8.axis_labels.hide_x()
f8.tick_labels.hide_x()

f9.axis_labels.hide_x()
f9.tick_labels.hide_x()

f10.axis_labels.hide_x()
f10.tick_labels.hide_x()

# Third row
f11.axis_labels.hide_x()

f12.axis_labels.hide_x()

f14.axis_labels.hide_x()

f15.axis_labels.hide_x()

# Hidden images :P

# Second row last column
f16.axis_labels.hide_x()
f16.tick_labels.hide_x()

# Third row last column
f17.axis_labels.hide_x()
f17.tick_labels.hide_x()

# Hide yaxis & ytick labels with aplpy for the appropriate subplots

# First row
f1.axis_labels.hide_y()

f2.axis_labels.hide_y()
f2.tick_labels.hide_y()

f3.axis_labels.hide_y()
f3.tick_labels.hide_y()

f4.axis_labels.hide_y()
f4.tick_labels.hide_y()

f5.axis_labels.hide_y()
f5.tick_labels.hide_y()

# Second row
f7.axis_labels.hide_y()
f7.tick_labels.hide_y()

f8.axis_labels.hide_y()
f8.tick_labels.hide_y()

f9.axis_labels.hide_y()
f9.tick_labels.hide_y()

f10.axis_labels.hide_y()
f10.tick_labels.hide_y()

# Third row
f11.axis_labels.hide_y()

f12.axis_labels.hide_y()
f12.tick_labels.hide_y()

f13.axis_labels.hide_y()
f13.tick_labels.hide_y()

f14.axis_labels.hide_y()
f14.tick_labels.hide_y()

f15.axis_labels.hide_y()
f15.tick_labels.hide_y()

# Hidden images :P

# Second row last column
f16.axis_labels.hide_y()
f16.tick_labels.hide_y()

# Third row last column
f17.axis_labels.hide_y()
f17.tick_labels.hide_y()

# ------------------------------------------------------------------------------

# DRAW

fig.canvas.draw()

#fig.savefig(galaxy_name+'_maps.png',dpi=800)

plt.show()

plt.close(fig)

# ------------------------------------------------------------------------------
