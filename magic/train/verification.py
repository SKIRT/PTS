#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.train.verification Contains the Verifier class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
from sklearn import svm
from sklearn.externals import joblib

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..core import Image, Source
from ...core.tools import filesystem, inspection

# -----------------------------------------------------------------

class Verifier(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(Verifier, self).__init__(config, "magic")

        # The classifier object
        self.classifier = None

        # Determine the path to the magic/classification user directory
        self.classification_user_path = os.path.join(inspection.pts_user_dir, "magic", "classification")

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Verifier instance
        verifier = cls()

        # -- Adjust the configuration settings according to the command-line arguments --

        verifier.config.mode = arguments.mode

        # Return the classifier
        return verifier
        
    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Verify
        self.verify()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup of the base class
        super(Classifier, self).setup()
        
    # -----------------------------------------------------------------

    def verify(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Starting the classification procedure ...")
        
    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Dump the classifier
        self.dump_classifier()
        
# -----------------------------------------------------------------
