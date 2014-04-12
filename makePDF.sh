#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory
#

# generate the latex documentation in a temporary folder next to the svn folder
/Applications/Doxygen.app/Contents/Resources/doxygen doc/pdf.doxygen

# generate the pdf file
make -C ../latex

# move the result
mkdir -p ../doc
mv -f ../latex/refman.pdf ../doc/PTS.pdf

# remove the intermediate files
rm -r ../latex
