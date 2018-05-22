#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.runner Contains the QueueRunner class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .queue import SimulationQueue
from ..tools import filesystem as fs
from ..simulation.execute import SkirtExec
from .jobscript import JobScript
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

class QueueRunner(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param config:
        :param interactive
        """

        # Call the constructor of the base class
        super(QueueRunner, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(QueueRunner, self).setup(**kwargs)

# -----------------------------------------------------------------

time_left = None

# Open the simulation queue
#queue_path = fs.join(pts_run_dir, name + ".queue")
#queue = SimulationQueue.from_file(queue_path)

# -----------------------------------------------------------------

# Create SKIRT execution environment
skirt = SkirtExec()

# -----------------------------------------------------------------

index = 0

#while time_left >= estimated_time:
while False:

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
#jobscript_path = fs.join()

#write("run_queue(...)") # start with runid = index

#jobscript = JobScript()

#jobscript.save()

# Submit job script

# qsub jobscript_path

# -----------------------------------------------------------------
