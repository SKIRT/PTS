#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory; use on Mac OS X only
#

# generate the html documentation in a folder next to the git folder
mkdir -p ../html
/Applications/Doxygen.app/Contents/Resources/doxygen doc/html.doxygen

# copy redirecting index.html file
cp doc/index_root.html ../html/index.html
