#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.train Train ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse
from sklearn import svm
from sklearn import datasets
from sklearn.externals import joblib

# Import the relevant PTS classes and modules
from pts.magic.core import Image, Source
from pts.core.tools import filesystem, inspection

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("--debug", type=str, help="debug mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Loop over all FITS files found in the current directory
for file_path in filesystem.files_in_path(os.getcwd(), extension="fits"):

    # Open the image
    image = Image(file_path, always_call_first_primary=False)

    image.frames.select_all()

    # Create a source
    source = Source.from_image(image)

    name = os.path.basename(file_path).split(".fits")[0]
    object_type, level, counter = name.split("_")

    plot_title = object_type + " with " + level

    # Plot the source
    source.plot(title=plot_title)

exit()

# Load test data
digits = datasets.load_digits()

classifier = svm.SVC(gamma=0.001, C=100.) # support vector classification

data = digits.data[:-10]
targets = digits.target[:-10]

test_data = digits.data[-10:]
test_targets = digits.target[-10:]

# The classifier is fit to the model, or learns from the model: by passing the training set
classifier.fit(data, targets)


# Serialize and dump the classifier
magic_user_path = os.path.join(inspection.pts_user_dir, "magic")
filesystem.create_directory(magic_user_path)

classifier_path = os.path.join(magic_user_path, "classifier.pkl")

# Dump the classifier
joblib.dump(classifier, classifier_path)

# Load the classifier
classifier2 = joblib.load(classifier_path)

# Predict the targets of the test data
result_targets = classifier2.predict(test_data)

for test_target, result_target in zip(test_targets, result_targets):
    
    print(test_target, result_target)

# -----------------------------------------------------------------
