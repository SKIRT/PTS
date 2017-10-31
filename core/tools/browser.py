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
import threading
from SimpleHTTPServer import SimpleHTTPRequestHandler, BaseHTTPServer

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

    return webbrowser.open(path, new=1)

# -----------------------------------------------------------------

def open_html(html):

    """
    This function ...
    :param html:
    :return:
    """

    temp_path = fs.join(introspection.pts_temp_dir, "page.html")
    fs.write_text(temp_path, html)
    return open_path(temp_path)

# -----------------------------------------------------------------

def open_page(page):

    """
    This function ...
    :param page:
    :return:
    """

    html = str(page)
    return open_html(html)

# -----------------------------------------------------------------

def start_localhost(port=8000, protocol="HTTP/1.0"):

    """
    This function ...
    :param port:
    :param protocol:
    :return:
    """

    server_address = ('', port)
    SimpleHTTPRequestHandler.protocol_version = protocol
    httpd = BaseHTTPServer.HTTPServer(server_address, SimpleHTTPRequestHandler)

    # Serve from other thread
    # httpd.serve_forever()
    thread = threading.Thread(target=httpd.serve_forever, args=[])
    thread.start()

    #httpd.serve_forever()

    # Return the thread?
    return thread

# -----------------------------------------------------------------

class serve_local_host():

    """
    This class ...
    """

    def __init__(self, port=8000, protocol="HTTP/1.0"):

        """
        The constructor ...
        :param port:
        """

        # Set properties
        self.port = port
        self.protocol = protocol

        # The server
        self.httpd = None

        # Thread
        self.thread = None

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        Thisfunction ...
        :return:
        """

        server_address = ('', self.port)
        SimpleHTTPRequestHandler.protocol_version = self.protocol
        self.httpd = BaseHTTPServer.HTTPServer(server_address, SimpleHTTPRequestHandler)

        # Serve from other thread
        #self.httpd.serve_forever()
        self.thread = threading.Thread(target=self.httpd.serve_forever, args=[])
        self.thread.start()

    # -----------------------------------------------------------------

    @property
    def sa(self):

        """
        This function ...
        :return:
        """

        return self.httpd.socket.getsockname()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print("Serving HTTP on", self.sa[0], "port", self.sa[1], "...")

    # -----------------------------------------------------------------

    def __exit__(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        #print("Stopping ...")

        # Stop serving
        # SHOULD WORK ...?
        #self.httpd.shutdown()

        finish = raw_input("Press ENTER to finish ...")
        exit()

# -----------------------------------------------------------------
