#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.prepareflybymovie Prepare a flyby movie by adjusting a particular ski file.
#
# This script adds the appropriate instruments to the specified ski file to create a flyby movie
# according to the specifications contained in one of the timelines hardcoded in this script.
# The first argument specifies the name of the timeline (currently one of "circle", "thru", "custom").
# The second argument specifies the file path of the ski file, or the keyword "interactive"
# to visualize the timeline interactively.

# -----------------------------------------------------------------

# Import standard modules
import sys

# Get the command line arguments
timelinename = sys.argv[1].lower() if len(sys.argv) > 1 else ""
skifilepath = sys.argv[2].lower() if len(sys.argv) > 2 else ""
interactive = (skifilepath=="interactive")

# -----------------------------------------------------------------

rate = 15
width = 21000
pixels = 1000

# -----------------------------------------------------------------

# construct the requested timeline
from pts.core.basics.timeline import Timeline
timeline = Timeline(rate=rate, lengthunit=(1 if interactive else width), shape=(pixels,pixels))

if timelinename == "circle":
    timeline.set(viewport=(1.5,1.5,1.5), crosshair=(0,0,0), upwards=(0.5,0.5,1), focal=2)
    timeline.advance((0,10), distance=1.5)
    timeline.circle((0,9), rotation=120, orientation=-15)

elif timelinename == "thru":
    timeline.set(viewport=(1.5,1.5,0), crosshair=(-1.5,-1.5,0), upwards=(0,0,1), focal=2)
    timeline.advance((0,10), distance=3)

elif timelinename == "custom":      # hardcode your own timeline
    timeline.set(viewport=(1,1,1), crosshair=(0,0,0), upwards=(0.5,0.5,1), focal=0.2)
    timeline.advance((0,9), distance=1.5)

else: raise ValueError("Unsupported timeline name")

# -----------------------------------------------------------------

# show interactively
if interactive:
    from pts.core.basics.viewangles import ViewAngles
    view = ViewAngles()
    view.flyby(timeline)
    view.close()

# or adjust ski file
else:
    from pts.core.plot.flybymovie import prepareflybymovie
    prepareflybymovie(skifilepath, timeline)
    print "adjusted " + skifilepath

# -----------------------------------------------------------------
