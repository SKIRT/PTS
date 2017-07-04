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

# Import astronomical modules
from astropy.io.fits import Header

import ipyvolume.pylab as p3
import ipyvolume
from ipyvolume.embed import embed_html

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import introspection
from pts.modeling.basics.models import load_2d_model, load_3d_model
from pts.magic.core.frame import Frame
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.models import DeprojectionModel3D
from pts.modeling.basics.properties import GalaxyProperties

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

def plot_galaxy_components(components, draw=True, show=True, shape=128, **kwargs):

    """
    This function ....
    :param components:
    :param draw:
    :param show:
    :param shape:
    :param kwargs:
    :return:
    """

    limits = [-10000., 10000.]

    x, y, z, r, theta, phi = xyz(shape=shape, limits=limits, spherical=True)
    data = r * 0

    # Loop over the components
    for name in components:

        component = components[name]
        density = component.density_function()(x, y, z)
        if name == "disk": density *= 40.
        data += density

        #print(name, np.mean(density))

    #kwargs["draw"] = True
    #kwargs["show"] = True

    # Return the plot
    return ipyvolume.quickvolshow(data=data.T, **kwargs)

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

components = {"disk": disk, "bulge": bulge}
box = plot_galaxy_components(components)
embed_html(filepath, box)

# -----------------------------------------------------------------

# Open HTML
fs.open_file(filepath)

# -----------------------------------------------------------------

# Remove
#fs.remove_file(filepath)

# -----------------------------------------------------------------
