#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.run_queue Run a particular simulation queue (intended to be run on a remote system)

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools.logging import log
from pts.core.remote.queue import SimulationQueue
from pts.core.tools import filesystem as fs
from pts.core.simulation.execute import SkirtExec
from pts.core.remote.jobscript import JobScript

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("name", "string", "name of the simulation queue")
definition.add_required("walltime", "duration", "walltime for the jobs")
definition.add_positional_optional("runid", "integer", "start with this run ID", default=0)

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("run_queue", "Run particular simulation queue")
config = setter.run(definition)

# -----------------------------------------------------------------

time_left = None

# Open the simulation queue
queue_path = fs.join(pts_run_dir, name + ".queue")
queue = SimulationQueue.from_file(queue_path)

# -----------------------------------------------------------------

# Create SKIRT execution environment
skirt = SkirtExec()

# -----------------------------------------------------------------

index = 0

while time_left >= estimated_time:

    # update the records to indicate that a run has been scheduled and submit the job;
    # do this within a single transaction context to ensure that the scheduled job sees the updated records
    # print("Submitting job to queue", config.queue, "for run-ids", runids)

    # with db.transaction():
    #    db.updatestatus(runids, 'running')
    #    db.updatefield(runids, 'queue', config.queue)
    #    subprocess.call(("bsub",), stdin=open(jobscriptname))

    # Run the simulations

    # Get the next simulation

    #runids = [index]

    # get the records corresponding to the run-id sequence from the database
    #records = queue.select("runid in (" + ",".join(map(str, runids)) + ")")
    #assert len(records) == 1

    #record = records[0]

    #queue.

    #queue = Database()
    runid = index
    rows = queue.select("runid = ?", (runid,))
    #db.close()
    if len(rows) != 1: raise ValueError("The specified run-id does not match a database record: " + str(runid))
    record = rows[0]

    # verify the runstatus of the database record
    if record['runstatus'] != 'scheduled':
        raise ValueError("The database record has run-status '" + record['runstatus'] + "' rather than 'scheduled'")

    # set the runstatus of the database record to 'running'
    #db = Database()
    with queue.transaction():
        queue.updatestatus((runid,), 'running')
    #db.close()


    # RUN

    skirt.run(definition, logging_options, parallelization, emulate, wait, silent)

    index += 1



# Launch new job if there are simulations left
jobscript_path = fs.join()

#write("run_queue(...)") # start with runid = index

jobscript = JobScript()

#jobscript.save()

# Submit job script

# qsub jobscript_path

# -----------------------------------------------------------------
