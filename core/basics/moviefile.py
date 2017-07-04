#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.moviefile Creating a movie from a sequence of stills
#
# An instance of the MovieFile class in this module allows creating a movie file from a sequence of images.

# -----------------------------------------------------------------

# Import standard modules
import os.path
import subprocess

# -----------------------------------------------------------------
#  MovieFile class
# -----------------------------------------------------------------

# TODO: Avoid mencoder dependency (which is difficult to install on Mac) by using matplotlib's MovieWriter
# (http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files)
# or something else ?

## An instance of the MovieFile class allows creating a movie file from a sequence of images. Use the constructor
# to specify the filename and the movie format, insert the images one at a time with the add() function, and finally
# write the movie header information to the file with the close() function.
#
# This class requires MEncoder to be installed and in the PATH.
# See the Installation Guide for more information.
#
# This class can generate movies in two distinct formats:
#  - MPEG-1 codec in \c .mpg container: this format is compatible with most movie players, including Windows
#    Media Player, Apple's QuickTime Player, and the free opensource player VLC. However the files are typically
#    several times larger than those created with more modern compression technologies such as MPEG-4. Also,
#    in this format the aspect ratio of the movie frame is fixed at 4/3 and the frame rate is fixed at 24 fps.
#  - MPEG-4 codec in \c .mov container: this format is playable with Apple's QuickTime Player (standard on Mac OS X
#    and available as a free download on Windows) and the free opensource player VLC. This format uses a modern
#    compression technology and it supports any aspect ratio and frame rate (probably within certain limits).
#
class MovieFile:
    ## The constructor launches \c mencoder to create a movie output file (replacing any existing file with the same
    # name). When the constructor finishes, \c mencoder is waiting for movie frames to be sent to it.
    #
    # The constructor accepts the following arguments:
    # - filepath: the filepath of the movie output file, which \em must end with either the .mpg
    #   or the .mov filename extension, determining the output movie format.
    # - shape: the number of movie frame pixels in the x and y directions; the default value is (800,600).
    #   For the .mpg format the aspect ratio \em must be 4/3, i.e. <tt>shape[0]/shape[1]==4/3</tt>.
    # - rate: the number of frames per second in the output movie; the default value is 24 fps.
    #   For the .mpg format the frame rate \em must be 24 fps.
    #
    def __init__(self, filepath, shape=(800,600), rate=24):
        # verify restrictions on arguments
        assert filepath.endswith((".mov",".mpg"))
        if filepath.endswith(".mpg"):
            assert 3*shape[0] == 4*shape[1]
            assert rate == 24

        # remember total number of bytes in pixel buffer for a single frame encoded as RGBA
        self.buflen = 4*shape[0]*shape[1]

        # ensure that we have access rights to create the output file (since we ignore any messages from mencoder)
        filepath = os.path.expanduser(filepath)
        open(filepath,'w').close()

        # construct the first part of the command line for raw video input (identical for both output formats)
        # note: type '$ mencoder -vf format=fmt=help' for a list of valid pixel formats / byte orderings
        cmdline = [ 'mencoder', '/dev/stdin', '-demuxer', 'rawvideo',
                    '-rawvideo', 'w=%i:h=%i'%shape + ':fps=%i:format=rgba'%rate,
                    '-really-quiet', '-o', filepath ]

        # add the appropriate options for .mov output
        if filepath.endswith(".mov"):
            cmdline += [
                    '-ovc', 'lavc',
                    '-lavcopts', 'vcodec=mpeg4:vbitrate=10000:mbd=2:cmp=2:subcmp=2:trell=yes:v4mv=yes:aic=2:vglobal=1',
                    '-ffourcc', 'mp4v', '-of', 'lavf',
                    '-lavfopts', 'format=mp4:i_certify_that_my_video_stream_does_not_use_b_frames'
                    ]

        # set the appropriate options for .mpg output
        if filepath.endswith(".mpg"):
            cmdline += [
                    '-ovc', 'lavc',
                    '-lavcopts', 'vcodec=mpeg1video:vbitrate=1152:keyint=15:mbd=2:aspect=4/3',
                    '-of', 'mpeg',
                    '-mpegopts', 'format=mpeg1:tsaf:muxrate=2000'
                    ]

        # launch mencoder; pipe the input from this process and pipe any messages to the null device
        self.p = subprocess.Popen(cmdline, stdin=subprocess.PIPE,
                                  stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

    # This function writes a string or a buffer containing the pixel data for a single movie frame to \c mencoder.
    # Each pixel must be encoded as RGBA (one byte for each R,G,B channel and an additional byte that is ignored).
    # The length of the buffer must match <tt>4*shape[0]*shape[1]</tt>.
    def add(self, frame):
        assert len(frame) == self.buflen
        self.p.stdin.write(frame)
        #self.p.communicate(input=frame) better? (is recommended)

   # This function writes the proper header information for the movie, and closes the movie file.
   # The function \em must be called for the movie file to be playable.
    def close(self):
        self.p.stdin.close()

# -----------------------------------------------------------------
