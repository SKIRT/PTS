#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_truncated Plot the truncated images for a certain factor.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import ipywidgets as widgets
from IPython.display import display
#from ipywidgets import embed_snippet, embed_minimal_html
import matplotlib.pyplot as plt
#from IPython.terminal.embed import embed_snippet
from ipyvolume.transferfunction import TransferFunctionWidgetJs3
from ipyvolume.style import dark, light

# Import astronomical modules
from astropy.io.fits import Header

import ipyvolume.pylab as p3
import ipyvolume
from ipyvolume.embed import embed_html

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools.logging import log
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter, parse_arguments
from pts.core.tools import introspection
from pts.modeling.basics.models import load_2d_model, load_3d_model
from pts.magic.core.frame import Frame
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.models import DeprojectionModel3D
from pts.modeling.basics.properties import GalaxyProperties

# -----------------------------------------------------------------

# Create log
definition = ConfigurationDefinition()
definition.add_optional("style", "string", "style: light or dark", "light", choices=["light", "dark"])
config = parse_arguments("test_m81_model_plot", definition)

# -----------------------------------------------------------------

def xyz(shape=128, limits=[-3, 3], spherical=False, sparse=True, centers=False):

    """
    This function ...
    :param shape:
    :param limits:
    :param spherical:
    :param sparse:
    :param centers:
    :return:
    """

    dim = 3

    try: shape[0]
    except: shape = [shape] * dim

    try: limits[0][0]
    except: limits = [limits] * dim

    if centers: v = [slice(vmin+(vmax-vmin)/float(N)/2, vmax-(vmax-vmin)/float(N)/4, (vmax-vmin)/float(N)) for (vmin, vmax), N in zip(limits, shape)]
    else: v = [slice(vmin, vmax+(vmax-vmin)/float(N)/2, (vmax-vmin)/float(N-1)) for (vmin, vmax), N in zip(limits, shape)]

    if sparse: x, y, z = np.ogrid.__getitem__(v)
    else: x, y, z = np.mgrid.__getitem__(v)

    # RETURN
    if spherical:

        r = np.linalg.norm([x, y, z])
        theta = np.arctan2(y, x)
        phi = np.arccos(z / r)
        return x, y, z, r, theta, phi

    else: return x, y, z

# -----------------------------------------------------------------

def plot_galaxy_components(components, draw=True, show=True, shape=128, unit="pc"):

    """
    This function ....
    :param components:
    :param draw:
    :param show:
    :param shape:
    :param unit:
    :param kwargs:
    :return:
    """

    #limits = [-10000., 10000.]

    x_min = 0.0
    x_max = 0.0
    y_min = 0.0
    y_max = 0.0
    z_min = 0.0
    z_max = 0.0

    #print("")
    # Determine limits, loop over components
    for name in components:

        component = components[name]
        x_min_scalar = component.xmin.to(unit).value
        x_max_scalar = component.xmax.to(unit).value
        y_min_scalar = component.ymin.to(unit).value
        y_max_scalar = component.ymax.to(unit).value
        z_min_scalar = component.zmin.to(unit).value
        z_max_scalar = component.zmax.to(unit).value

        if x_min_scalar < x_min: x_min = x_min_scalar
        if x_max_scalar > x_max: x_max = x_max_scalar
        if y_min_scalar < y_min: y_min = y_min_scalar
        if y_max_scalar > y_max: y_max = y_max_scalar
        if z_min_scalar < z_min: z_min = z_min_scalar
        if z_max_scalar > z_max: z_max = z_max_scalar

        #print(name + " limits: ")
        #print("")
        #print(" - x: " + str(component.xrange))
        #print(" - y: " + str(component.yrange))
        #print(" - z: " + str(component.zrange))
        #print("")

    minvalue = min(x_min, y_min, z_min)
    maxvalue = max(x_max, y_max, z_max)

    # Define limits
    #limits = [[x_min, x_max], [y_min, y_max], [z_min, z_max]]

    limits = [minvalue, maxvalue]

    # Debugging
    log.debug("Plot limits: " + str(limits))

    # Create coordinate data
    x, y, z, r, theta, phi = xyz(shape=shape, limits=limits, spherical=True)
    data = r * 0

    # Loop over the components
    for name in components:

        # Debugging
        log.debug("Computing the density of the " + name + " component ...")

        component = components[name]
        density = component.density_function(normalize=True)(x, y, z)
        data += density

    # DRAW FIGURE
    if draw:

        # :param lighting: boolean, to use lighting or not, if set to false, lighting parameters will be overriden
        # :param data_min: minimum value to consider for data, if None, computed using np.nanmin
        # :param data_max: maximum value to consider for data, if None, computed using np.nanmax
        # :param tf: transfer function (see ipyvolume.transfer_function, or use the argument below)
        # :param stereo: stereo view for virtual reality (cardboard and similar VR head mount)
        # :param width: width of rendering surface
        # :param height: height of rendering surface
        # :param ambient_coefficient: lighting parameter
        # :param diffuse_coefficient: lighting parameter
        # :param specular_coefficient: lighting parameter
        # :param specular_exponent: lighting parameter
        # :param downscale: downscale the rendering for better performance, for instance when set to 2, a 512x512 canvas will show a 256x256 rendering upscaled, but it will render twice as fast.
        # :param level: level(s) for the where the opacity in the volume peaks, maximum sequence of length 3
        # :param opacity: opacity(ies) for each level, scalar or sequence of max length 3
        # :param level_width: width of the (gaussian) bumps where the opacity peaks, scalar or sequence of max length 3
        # :param kwargs: extra argument passed to Volume and default transfer function

        # DEFAULT:

        # lighting=False, data_min=None, data_max=None, tf=None, stereo=False,
        # width=400, height=500,
        # ambient_coefficient=0.5, diffuse_coefficient=0.8,
        # specular_coefficient=0.5, specular_exponent=5,
        # downscale=1,
        # level=[0.1, 0.5, 0.9], opacity=[0.01, 0.05, 0.1], level_width=0.1,

        level = [0.2]
        opacity = [0.05, 0.0, 0.0]
        level_width = 0.2
        level_width = [level_width] * 3

        kwargs = dict()
        kwargs["width"] = 700
        kwargs["height"] = 800
        kwargs["stereo"] = False
        kwargs["level"] = level
        kwargs["opacity"] = opacity
        kwargs["level_width"] = level_width
        kwargs["downscale"] = 1

        # Create transfer function arguments
        tf_kwargs = {}

        # Clip off lists
        min_length = min(len(level), len(level_width), len(opacity))
        level = list(level[:min_length])
        opacity = list(opacity[:min_length])
        level_width = list(level_width[:min_length])
        # append with zeros
        while len(level) < 3:
            level.append(0)
        while len(opacity) < 3:
            opacity.append(0)
        while len(level_width) < 3:
            level_width.append(0)
        for i in range(1,4):
            tf_kwargs["level"+str(i)] = level[i-1]
            tf_kwargs["opacity"+str(i)] = opacity[i-1]
            tf_kwargs["width"+str(i)] = level_width[i-1]
        tf = TransferFunctionWidgetJs3(**tf_kwargs)

        # Set the transfer function
        kwargs["tf"] = tf

        # Set style
        if config.style == "dark": kwargs["style"] = dark
        elif config.style == "light": kwargs["style"] = light
        else: raise ValueError("Invalid style: " + config.style)

        # Create the volume plot
        vol = ipyvolume.quickvolshow(data=data.T, **kwargs)
        #vol = p3.volshow(data=data, **kwargs)

        # SHOW?
        if show:
            #p3.volshow()
            vol.show()
        return vol

    # ONLY RETURN THE DATA
    else: return data

# -----------------------------------------------------------------

# Create the configuration definition
#definition = ConfigurationDefinition()

# Get configuration
#setter = ArgumentConfigurationSetter("test_m81_model_plot")
#config = setter.run(definition)

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M81
m81_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "M81")

# -----------------------------------------------------------------

# Determine the path
path = fs.join(m81_data_path, "properties.dat")

# Load
properties = GalaxyProperties.from_file(path)

# -----------------------------------------------------------------

models_path = fs.join(m81_data_path, "models")
disk2d_path = fs.join(models_path, "disk.mod")
bulge2d_path = fs.join(models_path, "bulge.mod")

# -----------------------------------------------------------------

# Load the models
disk2d_model = load_2d_model(disk2d_path)
bulge2d_model = load_2d_model(bulge2d_path)

# Get the scale heights
old_scale_height = disk2d_model.scalelength / 8.26  # De Geyter et al. 2014
young_scale_height = 0.5 * old_scale_height
ionizing_scale_height = 0.25 * old_scale_height
dust_scale_height = 0.25 * old_scale_height

# -----------------------------------------------------------------

old_filename = "old_stars.fits"
young_filename = "young_stars.fits"
ionizing_filename = "ionizing_stars.fits"
dust_filename = "dust.fits"

# -----------------------------------------------------------------

# Determine paths
path = fs.join(m81_data_path, "components")
bulge_path = fs.join(path, "bulge.mod")
disk_path = fs.join(path, "disk.mod")

# -----------------------------------------------------------------

# Load bulge model
bulge = load_3d_model(bulge_path)

# Load disk model
disk = load_3d_model(disk_path)

# No y flattening: this is a mistake in the file
bulge.y_flattening = 1.

print("bulge:")
print(bulge)

print("disk:")
print(disk)

# -----------------------------------------------------------------

# Determine path to maps directory
maps_path = fs.join(m81_data_path, "maps")

# Determine the path to the header file
header_path = fs.join(maps_path, "header.txt")
header = Header.fromtextfile(header_path)
wcs = CoordinateSystem(header=header)

# Old stars
old_map_path = fs.join(maps_path, old_filename)
old_map = Frame.from_file(old_map_path)
old_map.wcs = wcs

# young stars
young_map_path = fs.join(maps_path, young_filename)
young_map = Frame.from_file(young_map_path)
young_map.wcs = wcs

# Ionizing stars
ionizing_map_path = fs.join(maps_path, ionizing_filename)
ionizing_map = Frame.from_file(ionizing_map_path)
ionizing_map.wcs = wcs

# Dust
dust_map_path = fs.join(maps_path, dust_filename)
dust_map = Frame.from_file(dust_map_path)
dust_map.wcs = wcs

# -----------------------------------------------------------------

# CREATE DEPROEJCTIONS (also sets the distance)
old_deprojection = DeprojectionModel3D.from_wcs(wcs, properties.center, properties.distance, properties.position_angle, properties.inclination, old_map_path, old_scale_height)
young_deprojection = DeprojectionModel3D.from_wcs(wcs, properties.center, properties.distance, properties.position_angle, properties.inclination, young_map_path, young_scale_height)
ionizing_deprojection = DeprojectionModel3D.from_wcs(wcs, properties.center, properties.distance, properties.position_angle, properties.inclination, ionizing_map_path, ionizing_scale_height)
dust_deprojection = DeprojectionModel3D.from_wcs(wcs, properties.center, properties.distance, properties.position_angle, properties.inclination, dust_map_path, dust_scale_height)

# -----------------------------------------------------------------

# Path for rendering
filename = "render.html"
filepath = fs.join(introspection.pts_temp_dir, filename)

#components = {"disk": disk, "bulge": bulge}
#components = {"old": old_deprojection, "bulge": bulge}
#components = {"ionizing": ionizing_deprojection}
#components = {"young": young_deprojection}
components = {"dust": dust_deprojection}

box = plot_galaxy_components(components, draw=True, show=False)

embed_html(filepath, box)

# -----------------------------------------------------------------

# Open HTML
fs.open_file(filepath)

# -----------------------------------------------------------------

# Remove
#fs.remove_file(filepath)

# -----------------------------------------------------------------
