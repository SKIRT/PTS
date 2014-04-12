#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory
#

# generate the html documentation in a temporary folder next to the svn folder
/Applications/Doxygen.app/Contents/Resources/doxygen doc/html.doxygen

# move Qt compressed help file
mkdir -p ../doc
mv -f ../html/PTS.qch ../doc/PTS.qch
