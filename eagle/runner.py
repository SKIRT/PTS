#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.runner Functions to actually perform a SKIRT run from the SKIRT-run database
#
# This module offers functions to actually perform a SKIRT simulation on an EAGLE galaxy according to the
# specifications defined in a particular SKIRT-runs database record.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

import os
import os.path
import shutil

from . import config as config
from .database import Database
from .galaxy import Snapshot, Galaxy
from .skirtrun import SkirtRun

from ..core.simulation.skifile import SkiFile
from ..core.plot.seds import plotseds
from ..core.plot.rgbimages import makergbimages

# -----------------------------------------------------------------

## This function actually performs a SKIRT simulation on an EAGLE galaxy according to the
# specifications defined in the SKIRT-runs database record with the specified run-id.
#
# Upon invocation of the function, the run-status of the specified record must be 'scheduled'. An error is raised
# if this is not the case. The function sets the run-status to 'running' while it is running, and it finally
# updates the run-status to 'completed' or 'failed' before returning.
#
# Specifically, the script performs the following steps:
#  - get the relevant information from the database record
#  - verify that the run-status of the database record is 'scheduled'
#  - update the run-status of the database record to become 'running'
#  - create the appropriate SKIRT result directories
#  - create an adjusted copy of the ski file for this run
#  - extract the particle data for the galaxy from the EAGLE snapshot
#  - run the SKIRT simulation
#  - create RGB images and/or SED plots from the results
#  - update the run-status of the database record to become 'completed' or 'failed'
#
def run(runid):
    # get the relevant information from the database record
    db = Database()
    rows = db.select("runid = ?", (runid,))
    db.close()
    if len(rows) != 1: raise ValueError("The specified run-id does not match a database record: " + str(runid))
    record = rows[0]

    # verify the runstatus of the database record
    if record['runstatus'] != 'scheduled':
        raise ValueError("The database record has run-status '" + record['runstatus'] + "' rather than 'scheduled'")

    try:
        # set the runstatus of the database record to 'running'
        db = Database()
        with db.transaction():
            db.updatestatus((runid,), 'running')
        db.close()

        # create the appropriate SKIRT result directories
        skirtrun = SkirtRun(runid, create=True)

        # extract the particle data for the galaxy from the EAGLE snapshot
        galaxies = Snapshot(record['eaglesim'], redshift=record['redshift']).galaxies()
        galaxy = galaxies.galaxy(record['galaxyid'])
        galaxy.export(skirtrun.inpath())

        # create an adjusted copy of the ski file for this run
        ski = SkiFile(os.path.join(config.templates_path, record['skitemplate']+".ski"))
        prefix = galaxy.prefix()
        ski.setstarfile(prefix + "_stars.dat")
        ski.setgasfile(prefix + "_gas.dat")
        ski.sethiifile(prefix + "_hii.dat")
        ski.saveto(os.path.join(skirtrun.runpath(), prefix+"_"+record['skitemplate']+".ski"))

        # copy the wavelength grid
        grid = ski.wavelengthsfile()
        shutil.copyfile(os.path.join(config.templates_path, grid), os.path.join(skirtrun.inpath(), grid))

        # run the SKIRT simulation
        simulation = skirtrun.execute(mpistyle=config.mpistyle,
                                      processes=config.nodes_per_job*config.processes_per_node,
                                      threads=config.threads_per_process)
        if simulation.status() != 'Finished':
            raise ValueError("SKIRT simulation " + simulation.status())

        # create SED plot
        plotseds(simulation)

        # create basic RGB images at the SDSS gri wavelengths (not integrated over bands)
        makergbimages(simulation, wavelength_tuples=((0.753,0.617,0.470),) )

        # move any .png and .pdf files to the visualization directory
        for visfile in filter(lambda fn: fn.endswith((".png",".pdf")), os.listdir(skirtrun.outpath())):
            os.rename(os.path.join(skirtrun.outpath(), visfile),
                      os.path.join(skirtrun.vispath(), visfile))

        # set the runstatus of the database record to 'completed'
        db = Database()
        with db.transaction():
            db.updatestatus((runid,), 'completed')
        db.close()

    except:
        # set the runstatus of the database record to 'failed'
        db = Database()
        with db.transaction():
            db.updatestatus((runid,), 'failed')
        db.close()
        raise

# -----------------------------------------------------------------
