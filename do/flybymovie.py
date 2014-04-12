#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.flybymovie Create a flyby movie for a particular ski file.
#
# This script creates a flyby movie for a particular ski file.
# The in/out filenames and other parameters are hardcoded in the script.
# Thus the script mainly serves as an example of how to use the pts.flybymovie module.
#

# -----------------------------------------------------------------

interactive = False

skifile = "galaxy_1_14_movie.ski"
moviefile = "galaxy_1_14.mov"
rate = 15
width = 21000
pixels = 1000
contrast = True

from pts.timeline import Timeline
timeline = Timeline(rate=rate, lengthunit=(1 if interactive else width), shape=(pixels,pixels))
timeline.set(viewport=(1,1,1), crosshair=(0,0,0), upwards=(0.5,0.5,1), focal=0.2)
timeline.glide((0,5), position=(-1,-1,1))
timeline.advance((0,5), distance=1)
timeline.set(5, crosshair=(1,1,-1))
timeline.advance((5,10), distance=0.9)

# -----------------------------------------------------------------

if interactive:
    from pts.viewangles import ViewAngles
    view = ViewAngles()
    view.flyby(timeline)
    view.close()

else:
    from pts.flybymovie import *
    adjustskifile(skifile, timeline)
    simulation = execskirt("skirt", skifile)
    createmovie(simulation, moviefile, timeline.rate, contrast=contrast)
    #createmovie("flyby_out", moviefile, timeline.rate, contrast=contrast)

# -----------------------------------------------------------------
