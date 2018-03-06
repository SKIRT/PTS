#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.component.sed Contains the ImagesModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.io.fits import Header

# Import the relevant PTS classes and modules
from .component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...magic.basics.coordinatesystem import CoordinateSystem
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class ImagesModelingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ImagesModelingComponent, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ImagesModelingComponent, self).setup()

        # Load the environment
        # NO WE NEED TO DO IT IN THE BASE CLASS BECAUSE E.G. THE FITTINGCOMPONENT DIRECTLY INHERITS FROM THIS CLASS BUT ALSO NEEDS THE ENVIRONMENT
        #self.environment = ImagesModelingEnvironment(self.config.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_header(self):

        """
        This function ...
        :return:
        """

        return Header.fromtextfile(self.environment.images_header_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.environment.images_header_path)

# -----------------------------------------------------------------

def get_images_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "images")

# -----------------------------------------------------------------

def get_images_paths(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_images_path(modeling_path), extension="fits")

# -----------------------------------------------------------------

def get_images_header_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(get_images_path(modeling_path), "header.txt")

# -----------------------------------------------------------------

def get_images_header(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return Header.fromtextfile(modeling_path)

# -----------------------------------------------------------------

def get_images_wcs(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return CoordinateSystem.from_file(get_images_header_path(modeling_path))

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
