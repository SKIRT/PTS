#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.simulator Performing the SKIRT simulation for a given SKIRT-run record
#
# The function in this module performs the SKIRT simulation on an EAGLE galaxy according to the
# specifications defined in a particular SKIRT-runs database record.

# -----------------------------------------------------------------

import os
import os.path
import shutil
import numpy as np
from astropy import units
import warnings

from ..core.tools.geometry import Transform
from ..core.simulation.skifile import SkiFile
from .skirtrun import SkirtRun
from . import config as config

# -----------------------------------------------------------------

## This function performs the SKIRT simulation on an EAGLE galaxy as described by the specified
# SKIRT-runs database record. The function uses the following fields in the specified record:
# runid, eaglesim, snaptag, galaxyid, skitemplate, numpp, deltamax.
#
# It is assumed that the data for the EAGLE galaxy has already been extracted from the snapshot
# and stored as text files in the appropriate SKIRT results directory.
#
def simulate(record):
    # get access to the appropriate SKIRT result directories
    skirtrun = SkirtRun(record["runid"])
    filenameprefix = "{}_{}_".format(record["eaglesim"], record["galaxyid"])

    # get the extent of the dust distribution defined in the input gas particles files
    # limited within a range of [0.1 kpc, 30 kpc]
    dustextent = max(0.1, min(30., getdustextent(os.path.join(skirtrun.inpath(),filenameprefix+"gas.dat"))))

    # get the rotation unit vector of the original galaxy from the info file in the input folder
    for line in open(os.path.join(skirtrun.inpath(),filenameprefix+"info.txt")):
        if not line.startswith("#"):
            segments = line.split()
            if segments[0]=="original_rotation_axis_x": a = float(segments[2])
            if segments[0]=="original_rotation_axis_y": b = float(segments[2])
            if segments[0]=="original_rotation_axis_z": c = float(segments[2])

    # calculate the angles for the "random" instrument, in degrees
    inc, azi, pos = randomInstrumentAngles(a, b, c)

    # create an adjusted copy of the ski file for this run
    ski = SkiFile(os.path.join(config.templates_path, record["skitemplate"]+".ski"))
    ski.setstarfile(filenameprefix + "stars.dat")
    ski.setgasfile(filenameprefix + "gas.dat")
    ski.sethiifile(filenameprefix + "hii.dat")
    ski.setpackages(record["numpp"])
    ski.setmaxmassfraction(record["deltamax"])
    ski.setdustextent(dustextent*units.kpc)
    ski.set_instrument_inclination("rn", inc*units.deg)
    ski.set_instrument_azimuth("rn", azi*units.deg)
    ski.set_instrument_pa("rn", pos*units.deg)
    ski.saveto(os.path.join(skirtrun.runpath(), filenameprefix+record["skitemplate"]+".ski"))

    # copy the wavelength grid
    grid = ski.wavelengthsfile()
    shutil.copyfile(os.path.join(config.templates_path, grid), os.path.join(skirtrun.inpath(), grid))

    # run the SKIRT simulation
    # (specify more than one process to enable MPI; actual nr of processes is defined by SLURM batch script)
    simulation = skirtrun.execute(mpistyle='srun', processes=2, threads=8)
    if simulation.status() != 'Finished':
        raise ValueError("SKIRT simulation " + simulation.status())

# -----------------------------------------------------------------

## This private helper function returns the instrument angles (inclination, azimuth, position) in degrees
# given the unit vector of the galaxy's angular momentum before the coordinate axes were lined up
def randomInstrumentAngles(a, b, c):
    # reconstruct the rotation matrix used for lining up the galaxy
    transf = Transform()
    v = np.sqrt(b*b+c*c)
    if v > 0.3:
        transf.rotateX(c/v, -b/v)
        transf.rotateY(v, -a)
    else:
        v = np.sqrt(a*a+c*c)
        transf.rotateY(c/v, -a/v)
        transf.rotateX(v, -b)

    # rotate two of the coordinate axes like the galaxy was rotated
    tx = transf.transform(1,0,0, 1)[0:3]
    tz = transf.transform(0,0,1, 1)[0:3]

    # calculate angles (following James' vodoo)
    inc    = np.arccos(tz[2])
    azi    = np.arctan2(tz[1], tz[0])
    alpha  = np.array((-np.sin(azi), np.cos(azi), 0))
    cosrot = np.dot(alpha, tx)
    sinrot = (np.cross(alpha, tx)**2).sum()**0.5
    pos    = (np.sign(90 - np.abs(azi))) *  np.arctan2(sinrot, cosrot)

    # convert the angles to degrees
    degree = 180/np.pi
    return inc*degree, azi*degree, pos*degree

# -----------------------------------------------------------------

## This private helper function returns the extent of the dust distribution defined in the specified gas
# particle file. Specifically it returns the largest half-size along the coordinate axes in kpc, rounded up to 100 pc.
def getdustextent(gaspath):
    extent = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = np.loadtxt(gaspath, usecols=(0,1,2,3,6), unpack=True)
    if data.shape[0]!=0:
        x,y,z,h,T = data
        mask = (T <= 8000.)
        if np.any(mask):
            xmin = np.min(x[mask] - h[mask])
            ymin = np.min(y[mask] - h[mask])
            zmin = np.min(z[mask] - h[mask])
            xmax = np.max(x[mask] + h[mask])
            ymax = np.max(y[mask] + h[mask])
            zmax = np.max(z[mask] + h[mask])
            extent = max(abs(xmin), abs(ymin), abs(zmin), abs(xmax), abs(xmax), abs(zmax))
            extent = np.ceil(extent/100.) / 10.

    return extent

# -----------------------------------------------------------------
