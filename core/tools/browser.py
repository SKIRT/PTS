#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.browser Provides functions for interacting with a web browser.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import webbrowser
webbrowser._tryorder = ["safari"]

# Import the relevant PTS classes and modules
from . import filesystem as fs
from . import introspection

# -----------------------------------------------------------------

def open_url(url):

    """
    This function ...
    :param url:
    :return:
    """

    return webbrowser.open(url, new=2)

# -----------------------------------------------------------------

def open_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    return webbrowser.open(path, new=2)

# -----------------------------------------------------------------

def open_html(html):

    """
    This function ...
    :param html:
    :return:
    """

    temp_path = fs.join(introspection.pts_temp_dir, "page.html")
    fs.write_text(temp_path, html)

# -----------------------------------------------------------------
