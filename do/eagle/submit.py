#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.submit Submit batch jobs for processing a given stage on a set of SKIRT-runs records
#
# This script prepares and submits batch jobs for processing a given stage on a set of SKIRT-runs records.
# Processing in the batch job occurs in one of two modes depending on the command line arguments:
#
#\verbatim
#pts eagle/submit <stage> loop <numjobs>
#pts eagle/submit <stage> force <runidspec>
#\endverbatim
#
# where:
#  - \<stage\> is one of "extract", "simulate", "observe", or "store";
#  - \<numjobs\> is the number of jobs to submit as a job array (or 1 for a single job);
#  - \<runidspec\> is a comma-seperated list of run-ids and/or run-id ranges (two run-ids with a dash in between).
#
# In loop mode, each batch job repeatedly performs the specified stage for a SKIRT-run record that is
# at the corresponding workflow stage and has the 'scheduled' status, until there are no more such records,
# or until the wall-time allocation of the job runs out.
#
# In force mode, the script performs the specified stage for each of the specified SKIRT-run records,
# regardless of the stage and status of those records.
#
# The current implementation submits all jobs to the Durham cosma6 queue using the SLURM queueing system.
#

# -----------------------------------------------------------------

# Import standard modules
import os.path
import subprocess
import sys

# Import the relevant PTS classes and modules
import pts.eagle.config as config
from pts.eagle.skirtrun import runids_in_range

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 4: raise ValueError("This script expects three command-line arguments: " \
                                        "(stage loop numjobs) or (stage force runidspec)")
stage = sys.argv[1]
if not stage in ("extract", "simulate", "observe", "store"):
    raise ValueError("First argument must be 'extract', 'simulate', 'observe', or 'store'")
mode = sys.argv[2]
if not mode in ("loop", "force"):
    raise ValueError("Second argument must be 'loop' or 'force'")
argum = sys.argv[3]

# parse and verify third argument according to requested mode
if mode == "loop":
    numjobs = int(argum)
    if numjobs<1 or numjobs>99: raise ValueError("Number of jobs must be in range [1,99]")
    print "Preparing {} {} jobs in loop mode".format(numjobs, stage)
if mode == "force":
    numjobs = 1
    runidspec = argum
    numruns = len(runids_in_range(runidspec))
    if numruns<1 or numruns>20: raise ValueError("Can't schedule more than 20 SKIRT-runs in single job")
    print "Preparing {} job with {} SKIRT-runs in force mode".format(stage, numruns)

# set some parameters according to requested stage
wallhours = 72
if stage=="simulate":
    tasks = 2
    cpuspertask = 8
    cpuspernode = 16
    memorypernode = 120  # can usually be lower for simulations without IFU
else:
    tasks = 1
    cpuspertask = 1
    cpuspernode = 1
    memorypernode = 30

# prepare an appropriate batch job script
jobscriptpath = os.path.join(config.jobs_path, "{}_jobscript.sh".format(stage))
jobscript = open(jobscriptpath,'w')
jobscript.write("#!/bin/bash -l\n")
if numjobs>1:
    jobscript.write("#SBATCH --array=1-{}\n".format(numjobs))
    jobscript.write("#SBATCH -o {}/{}_log_%A_%a.txt\n".format(config.jobs_path, stage))
    jobscript.write("#SBATCH -e {}/{}_log_%A_%a.txt\n".format(config.jobs_path, stage))
else:
    jobscript.write("#SBATCH -o {}/{}_log_%j.txt\n".format(config.jobs_path, stage))
    jobscript.write("#SBATCH -e {}/{}_log_%j.txt\n".format(config.jobs_path, stage))
jobscript.write("#SBATCH --ntasks={}\n".format(tasks))
jobscript.write("#SBATCH --cpus-per-task={}\n".format(cpuspertask))
jobscript.write("#SBATCH --time={}:00:00\n".format(wallhours))
jobscript.write("#SBATCH -J ES{}\n".format(stage))
jobscript.write("#SBATCH -D {}\n".format(config.jobs_path))
jobscript.write("#SBATCH -p cosma6\n")
jobscript.write("#SBATCH -A dp004\n")
jobscript.write("#SBATCH --exclusive\n")
jobscript.write("#SBATCH --mincpus={}\n".format(cpuspernode))
jobscript.write("#SBATCH --mem={}G\n".format(memorypernode))
jobscript.write("#SBATCH --hint=memory_bound\n")
jobscript.write("#SBATCH --mem_bind=local\n")
jobscript.write("echo Job started\n")
if mode=="loop":
    if numjobs>1:
        jobscript.write("echo sleep $SLURM_ARRAY_TASK_ID\n")  # spread out start times for job array members
        jobscript.write("sleep $SLURM_ARRAY_TASK_ID\n")  # spread out start times for job array members
    jobscript.write("python -u -m pts.do eagle/perform {} loop {}*3600\n".format(stage, wallhours-2))
if mode=="force":
    jobscript.write("python -u -m pts.do eagle/perform {} force {}\n".format(stage, runidspec))
jobscript.write("echo Job ended\n")
jobscript.write("sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode\n")
jobscript.write("exit\n")
jobscript.close()

# submit the batch job
subprocess.call(("sbatch",), stdin=open(jobscriptpath))

# -----------------------------------------------------------------010
