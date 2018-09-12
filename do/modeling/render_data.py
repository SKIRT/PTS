#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.render_data Render (or plot) 3D data

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from matplotlib import pyplot as plt

# Import astronomical modules
import yt
from yt.units import parsec, Msun

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.modeling.core.data import Data3D
from pts.core.basics.log import log
from pts.core.simulation.tree import get_cell_coordinates

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "3D data file path")
definition.add_optional("tree", "file_path", "tree filepath")
config = parse_arguments("render_data", definition)

# -----------------------------------------------------------------

def make_testdata():

    n_particles = 5000000

    ppx, ppy, ppz = 1e6*np.random.normal(size=[3, n_particles])

    ppm = np.ones(n_particles)

    #print(ppx.shape)
    #print(ppy.shape)
    #print(ppz.shape)
    #print(ppm.shape)

    return ppx, ppy, ppz, ppm

# -----------------------------------------------------------------

def make_data(ppx, ppy, ppz, ppm):

    """
    This function ...
    :param ppx:
    :param ppy:
    :param ppz:
    :param ppm:
    :return:
    """

    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_mass': ppm}

    bbox = 1.1 * np.array([[np.min(ppx), np.max(ppx)], [np.min(ppy), np.max(ppy)], [np.min(ppz), np.max(ppz)]])
    return data, bbox

# -----------------------------------------------------------------

def plot(data, bbox):

    """
    This function ...
    :param data:
    :param bbox:
    :return:
    """

    #yt.load_amr_grids()

    ds = yt.load_particles(data, length_unit=parsec, mass_unit=1e8*Msun, n_ref=256, bbox=bbox)

    #xwidth = (data.x_span,"pc",)
    #ywidth = (data.y_span,"pc",)
    xwidth = (bbox[0][1]-bbox[0][0], "pc",)
    ywidth = (bbox[1][1]-bbox[1][0], "pc",)

    #print(ds.derived_field_list)

    what = ('deposit', 'io_density')
    slc = yt.SlicePlot(ds, "z", what)
    slc.set_width((xwidth,ywidth,))

    #slc.show() # only for ipython notebook? :(
    slc.save("test.pdf")

# -----------------------------------------------------------------

# Load the data
data = Data3D.from_file(config.filename)
#_x, _y, _z, values = data.valid_x, data.valid_y, data.valid_z, data.valid_values
_x, _y, _z, values = data.x, data.y, data.z, data.values

#x, y, z, values = make_testdata()
#values = np.ones(data.nx)
#print(x.shape, type(x), x.mask, np.sum(x.mask))
#print(y.shape, type(y), y.mask, np.sum(y.mask))
#print(z.shape, type(z), z.mask, np.sum(z.mask))
#x, y, z = np.random.normal(size=[3, data.nx])
#print(x.shape, type(x))
#print(y.shape, type(y))
#print(z.shape, type(z))
#data, bbox = make_data(x, y, z, values)
#plot(data, bbox)

# -----------------------------------------------------------------

# Define bounding box
#bbox = 1.1*np.array([[data.min_x, data.max_x], [data.min_y, data.max_y], [data.min_z, data.max_z]])

# -----------------------------------------------------------------

#print(data.valid_x.shape)
#print(data.valid_y.shape)
#print(data.valid_z.shape)
#print(data.valid_values.shape)
#x, y, z, values = data.valid_x, data.valid_y, data.valid_z, data.valid_values

# Create the particle data set
#dat = {"particle_position_x":data.x, "particle_position_y":data.y, "particle_position_z":data.z, "sfr":data.values}
#dat = {"particle_position_x":x, "particle_position_y":y, "particle_position_z":z, "particle_mass":values}
#dat = {"particle_position_x":data.valid_x, "particle_position_y":data.valid_y, "particle_position_z":data.valid_z, ("io", "StarFormationRate"):data.valid_values}
#ds = yt.load_particles(dat, length_unit=parsec, bbox=bbox, mass_unit=1e8*Msun, n_ref=256)

# -----------------------------------------------------------------

#print(ds.field_list)
#print(ds.derived_field_list)

# -----------------------------------------------------------------

#slc = yt.SlicePlot(ds, 2, ('deposit', 'all_cic'))
#slc.set_width((8, 'Mpc'))
#slc.show_axes()
#plt.show()
#slc.show()
#slc.save("test.pdf")

#yt.ParticlePlot(ds, "particle_position_x", "particle_position_y", z_fields=None)

# -----------------------------------------------------------------

# grid_data = [dict(left_edge = [0.0, 0.0, 0.0],
#               right_edge = [1.0, 1.0, 1.],
#               level = 0,
#               dimensions = [32, 32, 32],
#               number_of_particles = 0),
#          dict(left_edge = [0.25, 0.25, 0.25],
#               right_edge = [0.75, 0.75, 0.75],
#               level = 1,
#               dimensions = [32, 32, 32],
#               number_of_particles = 0)]

#
#for g in grid_data: g["density"] = (np.random.random(g["dimensions"])*2**g["level"], "g/cm**3")
#ds = load_amr_grids(grid_data, [32, 32, 32], length_unit=1.0)

#

domain_dimensions = [32, 32, 32]

# Get
xmin, xmax, ymin, ymax, zmin, zmax = get_cell_coordinates(config.tree, read_method="pandas")

ncells = len(xmin)
#grid_data = []
#for index in range(ncells):
#cell = dict(left_edge=[xmin[index],ymin[index],zmin[index]], right_edge=[xmax[index],ymax[index],zmax[index]],
#                dimensions=domain_dimensions, density=values[index], level=0)
#    grid_data.append(cell)

grid_data = dict(left_edge=[xmin,ymin,zmin], right_edge=[xmax,ymax,zmax], level=0, dimensions=domain_dimensions)

bbox = None
ds = yt.load_amr_grids(grid_data, domain_dimensions, bbox=bbox)

# load_amr_grids(grid_data, domain_dimensions,
#                    bbox=None, sim_time=0.0, length_unit=None,
#                    mass_unit=None, time_unit=None, velocity_unit=None,
#                    magnetic_unit=None, periodicity=(True, True, True),
#                    geometry="cartesian", refine_by=2, unit_system="cgs"):
#     r"""Load a set of grids of data into yt as a
#     :class:`~yt.frontends.stream.data_structures.StreamHandler`.
#     This should allow a sequence of grids of varying resolution of data to be
#     loaded directly into yt and analyzed as would any others.  This comes with
#     several caveats:
#
#     * Units will be incorrect unless the unit system is explicitly specified.
#     * Some functions may behave oddly, and parallelism will be
#       disappointing or non-existent in most cases.
#     * Particles may be difficult to integrate.
#     * No consistency checks are performed on the index
#
#     Parameters
#     ----------
#
#     grid_data : list of dicts
#         This is a list of dicts. Each dict must have entries "left_edge",
#         "right_edge", "dimensions", "level", and then any remaining entries are
#         assumed to be fields. Field entries must map to an NDArray. The grid_data
#         may also include a particle count. If no particle count is supplied, the
#         dataset is understood to contain no particles. The grid_data will be
#         modified in place and can't be assumed to be static.
