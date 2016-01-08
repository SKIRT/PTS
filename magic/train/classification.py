#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.train.classification Contains the Classifier class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
from sklearn import svm
from sklearn import datasets
from sklearn.externals import joblib

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..core import Image, Source
from ...core.tools import filesystem, inspection

# -----------------------------------------------------------------

class Classifier(Configurable):

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
        super(Classifier, self).__init__(config, "magic")

        # Determine the path to the magic/classification user directory
        self.collection_user_path = os.path.join(inspection.pts_user_dir, "magic", "collection")
        self.classification_user_path = os.path.join(inspection.pts_user_dir, "magic", "classification")

        # Create the user classification directory
        filesystem.create_directory(self.classification_user_path)

        self.yes_path = None
        self.no_path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Classifier instance
        classifier = cls()

        # -- Adjust the configuration settings according to the command-line arguments --

        classifier.mode = arguments.mode

        # Return the classifier
        return classifier

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Classify
        self.classify()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup of the base class
        super(Classifier, self).setup()

        # Determine the paths to the star collection and saturation collection directories
        mode_path = os.path.join(self.collection_user_path, self.config.mode)

        # Determine the paths to the 'yes' and 'no' saturation collection directories
        self.yes_path = os.path.join(mode_path, "yes")
        self.no_path = os.path.join(mode_path, "no")

    # -----------------------------------------------------------------

    def classify(self):

        """
        This function ...
        :return:
        """

        # Loop over all FITS files found in the current directory
        for file_path in filesystem.files_in_path(os.getcwd(), extension="fits"):

            # Open the image
            image = Image(file_path, always_call_first_primary=False)
            image.frames.select_all()

            # Create a source
            source = Source.from_image(image)

            name = os.path.basename(file_path).split(".fits")[0]
            object_type, level, counter = name.split("_")

            # Skip galaxies for now
            if object_type == "galaxy": continue

            # Get description
            description = object_type + " (" + level + ")"

            # Get target class
            target = classes[description]

            plot_title = description + " (class " + str(target) + ")"

            # Plot the source
            source.plot(title=plot_title)

        # Load test data
        #digits = datasets.load_digits()

        classifier = svm.SVC(gamma=0.001, C=100.) # support vector classification

        data = digits.data[:-10]
        targets = digits.target[:-10]

        test_data = digits.data[-10:]
        test_targets = digits.target[-10:]



        # The classifier is fit to the model, or learns from the model: by passing the training set
        classifier.fit(data, targets)


        # Serialize and dump the classifier


        # Determine the path to the pickle file
        classifier_path = os.path.join(self.classification_user_path, "classifier.pkl")

        # Dump the classifier
        joblib.dump(classifier, classifier_path)

        # Load the classifier
        classifier2 = joblib.load(classifier_path)

        # Predict the targets of the test data
        result_targets = classifier2.predict(test_data)

        for test_target, result_target in zip(test_targets, result_targets):

            print(test_target, result_target)

    # -----------------------------------------------------------------

    def is_star(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        return False

    # -----------------------------------------------------------------

    def is_galaxy(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        return False

# -----------------------------------------------------------------
