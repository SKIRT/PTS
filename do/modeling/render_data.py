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

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "3D data file path")
config = parse_arguments("render_data", definition)

# -----------------------------------------------------------------

def test():

    n_particles = 5000000

    ppx, ppy, ppz = 1e6*np.random.normal(size=[3, n_particles])

    ppm = np.ones(n_particles)

    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_mass': ppm}
    bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
    ds = yt.load_particles(data, length_unit=parsec, mass_unit=1e8*Msun, n_ref=256, bbox=bbox)

    slc = yt.SlicePlot(ds, 2, ('deposit', 'all_cic'))
    slc.set_width((8, 'Mpc'))
    slc.show_axes()
    plt.show()

test()
exit()

# -----------------------------------------------------------------

# Load the data
data = Data3D.from_file(config.filename)

# -----------------------------------------------------------------

# Define bounding box
bbox = 1.1*np.array([[data.min_x, data.max_x], [data.min_y, data.max_y], [data.min_z, data.max_z]])

# -----------------------------------------------------------------

# Create the particle data set
#dat = {"particle_position_x":data.x, "particle_position_y":data.y, "particle_position_z":data.z, "sfr":data.values}
#dat = {"particle_position_x":data.x, "particle_position_y":data.y, "particle_position_z":data.z, "particle_mass":data.values}
dat = {"particle_position_x":data.valid_x, "particle_position_y":data.valid_y, "particle_position_z":data.valid_z, ("io", "StarFormationRate"):data.valid_values}
ds = yt.load_particles(dat, length_unit=parsec, bbox=bbox)

# -----------------------------------------------------------------

print(ds.field_list)
print(ds.derived_field_list)

# -----------------------------------------------------------------

#slc = yt.SlicePlot(ds, 2, ('deposit', 'all_cic'))
yt.ParticlePlot(ds, "particle_position_x", "particle_position_y", z_fields=None)

# -----------------------------------------------------------------
