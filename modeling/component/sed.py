#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.component.sed Contains the SEDModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from .component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.data.sed import ObservedSED

# -----------------------------------------------------------------

class SEDModelingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDModelingComponent, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDModelingComponent, self).setup()

        # Load the environment
        # NO WE NEED TO DO IT IN THE BASE CLASS BECAUSE E.G. THE FITTINGCOMPONENT DIRECTLY INHERITS FROM THIS CLASS BUT ALSO NEEDS THE ENVIRONMENT
        #self.environment = SEDModelingEnvironment(self.config.path)

# -----------------------------------------------------------------

def get_observed_sed_file_path(modeling_path):

    """
    This function ...
    :return:
    """

    return fs.join(modeling_path, "sed.dat")

# -----------------------------------------------------------------

def get_observed_sed(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return ObservedSED.from_file(get_observed_sed_file_path(modeling_path))

# -----------------------------------------------------------------

def get_sed_plot_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "sed.pdf")

# -----------------------------------------------------------------

def get_ski_template_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "template.ski")

# -----------------------------------------------------------------

def get_ski_template(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """
    return SkiFile(get_ski_template_path(modeling_path))

# -----------------------------------------------------------------

def get_ski_input_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "input")

# -----------------------------------------------------------------

def get_ski_input_paths(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_ski_input_path(modeling_path)
    if not fs.is_directory(path): return None
    else: return fs.files_in_path(path)

# -----------------------------------------------------------------
