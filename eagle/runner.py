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

import os
import os.path
import shutil
import eagle.config as config
from eagle.database import Database
from eagle.galaxy import Snapshot, Galaxy
from eagle.skirtrun import SkirtRun
from pts.skifile import SkiFile
from pts.plotseds import plotseds   # import this one first to avoid warnings about setting the matplotlib backend
from pts.makergbimages import makeintegratedrgbimages
from pts.filter import Filter

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
        ski.saveto(os.path.join(skirtrun.runpath(), skiname))

        # copy the wavelength grid
        grid = record['skitemplate']+"_wavelengths.dat"
        shutil.copyfile(os.path.join(config.templates_path, grid), os.path.join(skirtrun.inpath(), grid))

        # extract the particle data for the galaxy from the EAGLE snapshot
        galaxies = Snapshot(record['eaglesim'], redshift=record['redshift']).galaxies()
        galaxy = galaxies.galaxy(record['groupnr'], record['subgroupnr'])
        galaxy.export(skirtrun.inpath())

        # run the SKIRT simulation
        simulation = skirtrun.execute(mpistyle=config.mpistyle,
                                      processes=config.nodes_per_job*config.processes_per_node,
                                      threads=config.threads_per_process)
        if simulation.status() != 'Finished':
            raise ValueError("SKIRT simulation " + simulation.status())

        # create SED plot
        plotseds(simulation)

        # create RGB images integrated over the SDSS gri bands, and with enhanced IR and and FUV channels
        filterR = Filter('SDSS.i')
        filterG = Filter('SDSS.r')
        filterB = Filter('SDSS.g')
        filterIR = Filter('Pacs.green')
        filterUV = Filter('GALEX.FUV')
        fmin,fmax = makeintegratedrgbimages(simulation,
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1) ] )
        makeintegratedrgbimages(simulation,
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1), (filterIR, 0.02,0,0), (filterUV, 0,0,4) ],
            postfix="2", fmin=fmin, fmax=fmax)

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
