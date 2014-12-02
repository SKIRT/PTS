#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory
#

##############################
## On the Mac OS X platform ##
##############################

if [ "$(uname)" == "Darwin" ]
then

  # Generate html documentation in a temporary folder next to the git folder
  /Applications/Doxygen.app/Contents/Resources/doxygen doc/html.doxygen

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

  #####################################################
  ### Search for the most recent version of Qt first ##
  #####################################################

  if [ -f $HOME/Qt5.3.2/5.3/clang_64/bin/qhelpgenerator ]
  then

    # generate the Qt compressed help file
    $HOME/Qt5.3.2/5.3/clang_64/bin/qhelpgenerator ../html/index_mathjax.qhp -o ../doc/PTS.qch
    
  ############################################################
  ## Search for the older Qt version that is also supported ##
  ############################################################
    
  elif [ -f $HOME/Qt5.2.1/5.2.1/clang_64/bin/qhelpgenerator ]
  then
    
    # generate the Qt compressed help file
    $HOME/Qt5.2.1/5.2.1/clang_64/bin/qhelpgenerator ../html/index_mathjax.qhp -o ../doc/PTS.qch
    
  else
    echo "Error: could not find the Qt help file generator. If Qt is installed in a custom location, change this file accordingly."	
  fi
    
##########################
# On the Linux platform ##
##########################

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]
then
    
  # generate html documentation in a temporary folder next to the git folder
  doxygen doc/html_ubuntu.doxygen
    
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
    
  #####################################################
  ### Search for the most recent version of Qt first ##
  #####################################################
    
  if [ -f $HOME/Qt5.3.2/5.3/gcc_64/bin/qhelpgenerator ]
  then

    # generate the Qt compressed help file
    $HOME/Qt5.3.2/5.3/gcc_64/bin/qhelpgenerator ../html/index_mathjax.qhp -o ../doc/PTS.qch

  ############################################################
  ## Search for the older Qt version that is also supported ##
  ############################################################

  elif [ -f $HOME/Qt5.2.1/5.2.1/gcc_64/bin/qhelpgenerator ]
  then

    # generate the Qt compressed help file
    $HOME/Qt5.2.1/5.2.1/gcc_64/bin/qhelpgenerator ../html/index_mathjax.qhp -o ../doc/PTS.qch
  
  else
    echo "Error: could not find the Qt help file generator. If Qt is installed in a custom location, change this file accordingly."	
  fi

else
  echo "Platforms other than Mac OS X or Ubuntu are not supported!"
fi
