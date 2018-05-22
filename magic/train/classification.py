#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.train.classification Contains the Classifier class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
from sklearn import svm
from sklearn.externals import joblib

# Import the relevant PTS classes and modules
from ..core.image import Image
from ..core.source import Source
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class Classifier(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(Classifier, self).__init__(*args, **kwargs)

        # The classifier object
        self.vector_classifier = None

        # Determine the path to the magic/classification user directory
        self.collection_user_path = os.path.join(introspection.pts_user_dir, "magic", "collection")
        self.classification_user_path = os.path.join(introspection.pts_user_dir, "magic", "classification")

        # Create the user classification directory
        fs.create_directory(self.classification_user_path)

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

        classifier.config.mode = arguments.mode

        # Return the classifier
        return classifier

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new classifier instance
        classifier = cls()

        # Load the classifier
        classifier.vector_classifier = joblib.load(path)

        # Return the classifier
        return classifier

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the training data
        self.load_data()

        # 3. Classify
        self.classify()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup of the base class
        super(Classifier, self).setup(**kwargs)

        # Create the vector classifier
        self.vector_classifier = svm.SVC(gamma=0.001, C=100.) # support vector classification

        # Determine the path to the collection directory for the current mode
        collection_mode_path = os.path.join(self.collection_user_path, self.config.mode)

        # Determine the paths to the 'yes' and 'no' saturation collection directories
        self.yes_path = os.path.join(collection_mode_path, "yes")
        self.no_path = os.path.join(collection_mode_path, "no")

        # Determine the path to the classification directory for the current mode
        self.classification_mode_path = os.path.join(self.classification_user_path, self.config.mode)

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Starting the classification procedure ...")

        # Get a list of the filepaths for every FITS file in the 'yes' directory
        yes_paths = fs.files_in_path(self.yes_path, extension="fits")

        # Get a list of the filepaths for every FITS file in the 'no' directory
        no_paths = fs.files_in_path(self.no_path, extension="fits")

        # Create a dictionary to contain all the file paths
        paths = {"yes": yes_paths, "no": no_paths}

        # Create the data structure
        self.rows = []
        self.targets = []

        # Inform the user
        log.info("Gathering files ...")

        # Loop over all FITS files
        for label in paths:

            # Target should be 1 for 'yes', 0 for 'no'
            if label == "yes": target = 1
            else: target = 0

            # Loop over the paths classified under the current label
            for path in paths[label]:

                # Open the image, select all frames
                image = Image.from_file(path, always_call_first_primary=False)
                image.frames.select_all()

                # Create a source
                source = Source.from_image(image)

                if source.has_background:

                    array = source.subtracted
                    array[source.background_mask] = 0.0

                else:

                    array = source.cutout
                    array[source.background_mask] = 0.0

                # Add the flattened data and target
                self.rows.append(array.flatten())
                self.targets.append(target)

        # Make the input data array
        self.data = np.vstack(self.rows)

        # Set nan values to zero
        self.data[np.isnan(self.data)] = 0.0

    # -----------------------------------------------------------------

    def classify(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting training data ...")

        # The classifier is fit to the model, or learns from the model: by passing the training set
        self.vector_classifier.fit(self.data, self.targets)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Dump the classifier
        self.dump_classifier()

    # -----------------------------------------------------------------

    def dump_classifier(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the pickle file
        classifier_path = os.path.join(self.classification_mode_path, "classifier.pkl")

        # Inform the user
        log.info("Writing the classifier to " + classifier_path)

        # Serialize and dump the classifier
        joblib.dump(self.vector_classifier, classifier_path)

    # -----------------------------------------------------------------

    def predict(self, data, single=False):

        """
        This function ...
        :return:
        """

        if single: data = data.reshape(1, -1)

        # Predict the targets of the test data
        result = self.vector_classifier.predict(data)

        # Return the result
        if single: return result[0]
        else: return result

    # -----------------------------------------------------------------

    def is_star(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        if source.has_background:

            data = source.subtracted
            data[source.background_mask] = 0.0
            data[np.isnan(data)] = 0.0

        else:

            data = source.cutout
            data[source.background_mask] = 0.0
            data[np.isnan(data)] = 0.0

        label = self.predict(data, single=True)
        return label == 1

    # -----------------------------------------------------------------

    def is_galaxy(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        return False

# -----------------------------------------------------------------
