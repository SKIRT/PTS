#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory
#

# This function returns the Qt version associated with a certain qhelpgenerator executable
function qt_version {

    local VERSION_OUTPUT="$($1 -v | tr -d ' ')"
    local SPLITTED="$(sed s/'(Qt'/' '/g <<< $VERSION_OUTPUT)"
    local LIST=($SPLITTED)
    local SECOND_PART=${LIST[1]}
    local SPLITTED="$(sed s/')'/' '/g <<< $SECOND_PART)"
    local LIST=($SPLITTED)
    local VERSION=${LIST[0]}
    echo $VERSION
}


# On the Mac OS X platform
if [ "$(uname)" == "Darwin" ]
then

  # Generate html documentation in a temporary folder next to the git folder
  /Applications/Doxygen.app/Contents/Resources/doxygen doc/html.doxygen

# On the Linux platform
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]
then

  # Generate html documentation in a temporary folder next to the git folder
  doxygen doc/html_ubuntu.doxygen

# Unsupported platform
else
  echo "Platforms other than Mac OS X or Ubuntu are not supported!"
  exit
fi


# Copy the 'mouse over' SKIRT logo
cp doc/images/SkirtLogoSmall-home.png ../html/SkirtLogoSmall-home.png

# Add the MathJax script to the index.qhp file -> index_mathjax.qhp
python doc/enable_qch_mathjax.py

# Obtain the MathJax repository if it is not yet present
if [ ! -d ../html/mathjax ]; then

  # Clone the repository and checkout version 2.4
  git clone git://github.com/mathjax/MathJax.git ../html/mathjax
  git -C ../html/mathjax checkout -b v2.4-latest origin/v2.4-latest

  # Remove unnecessary files and folders
  xargs -I fname rm -r fname < doc/mathjax_delete.txt

fi

# Search for qhelpgenerator in the home directory
PATHLIST="$(find $HOME/Qt* -name qhelpgenerator -type f | tr '\n' ' ')"

# Search for qhelpgenerator in the $PATH
PATHLIST="$(which qhelpgenerator) $PATHLIST"

# Set the QHELPPATH to an empty string initially
QHELPPATH=""

# Loop over all the qhelpgenerator paths
for path in $PATHLIST; do

  # Get the associated Qt version
  VERSION="$(qt_version $path)"

  # Check whether the Qt version is supported
  if [[ $VERSION > '5.2.0' ]]
  then
      # If another supported qhelpgenerator was found, check whether this qhelpgenerator corresponds to
      # a more recent Qt installation
      if [[ ! "$QHELPPATH" == "" ]]
      then
        CURRENT_VERSION="$(qt_version $path)"
        if [[ $VERSION > $CURRENT_VERSION ]]
        then
            QHELPPATH=$path
        fi

      # Else, just use this qhelpgenerator
      else
        QHELPPATH=$path
      fi
  fi

done

# Exit with an error if we don't find a supported qhelpgenerator
if [ "$QHELPPATH" == "" ]
then
  echo error: Could not find a qhelpgenerator executable associated with a supported Qt version
  exit
else
  echo Using qhelpgenerator in $QHELPPATH
  # Generate the Qt compressed help file
  $QHELPPATH ../html/index_mathjax.qhp -o ../doc/PTS.qch
fi
