#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.runner Functions to actually perform a SKIRT run from the SKIRT-run database
#
# This module offers functions to actually perform a SKIRT simulation on an EAGLE galaxy according to the
# specifications defined in a particular SKIRT-runs database record.

# -----------------------------------------------------------------

import os.path
import shutil
import eagle.config as config
from eagle.database import Database
from eagle.galaxy import Snapshot, Galaxy
from eagle.skirtrun import SkirtRun
from pts.skifile import SkiFile
from pts.makergbimages import makergbimages
from pts.plotseds import plotseds

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

        # construct some names based on the galaxy group/subgroup numbers
        starsname = "galaxy_{0}_{1}_stars.dat".format(record['groupnr'], record['subgroupnr'])
        gasname = "galaxy_{0}_{1}_gas.dat".format(record['groupnr'], record['subgroupnr'])
        hiiname = "galaxy_{0}_{1}_hii.dat".format(record['groupnr'], record['subgroupnr'])
        skiname = "{0}_{1}_{2}.ski".format(record['skitemplate'], record['groupnr'], record['subgroupnr'])

        # create an adjusted copy of the ski file for this run
        ski = SkiFile(os.path.join(config.templates_path, record['skitemplate']+".ski"))
        ski.setstarfile(starsname)
        ski.setgasfile(gasname)
        ski.sethiifile(hiiname)
        ski.setpackages(100000 + 5*(record['starparticles']+record['gasparticles']))
        ski.saveto(os.path.join(skirtrun.runpath(), skiname))

        # extract the particle data for the galaxy from the EAGLE snapshot
        galaxies = Snapshot(record['eaglesim'], redshift=record['redshift']).galaxies()
        galaxy = galaxies.galaxy(record['groupnr'], record['subgroupnr'])
        galaxy.export(skirtrun.inpath())

        # run the SKIRT simulation
        simulation = skirtrun.execute()
        simulation.removetemporaryfiles()
        if simulation.status() != 'Finished':
            raise ValueError("SKIRT simulation " + simulation.status())

        # create RGB images and/or SED plots
        makergbimages(simulation)
        plotseds(simulation)

        # move any .png and .pdf files to the visualization directory
        for visfile in filter(lambda fn: fn.endswith((".png",".pdf")), os.listdir(skirtrun.outpath())):
            shutil.move(os.path.join(skirtrun.outpath(), visfile), skirtrun.vispath())

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
