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

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import LabeledSkiFile
from ...core.data.sed import ObservedSED

# -----------------------------------------------------------------

class SEDModelingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SEDModelingComponent, self).__init__(config)

        # The SED file path
        self.sed_path = None

        # The observed SED plot file path
        self.sed_plot_path = None

        # The ski template path
        self.ski_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDModelingComponent, self).setup()

        # Set the SED path
        self.sed_path = get_observed_sed_file_path(self.config.path)

        # Set the SED plot path
        self.sed_plot_path = fs.join(self.config.path, "sed.pdf")

        # Set the ski template path
        self.ski_path = get_ski_template_path(self.config.path)

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

    return LabeledSkiFile(get_ski_template_path(modeling_path))

# -----------------------------------------------------------------
