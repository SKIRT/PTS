#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.calibration Contains the CalibrationErrors class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import urllib
import requests
import subprocess
from lxml import html

# Import the relevant PTS classes and modules
from ...core.tools import inspection
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.tools import tables
from ...core.tools import serialization

# -----------------------------------------------------------------

# The path to the PTS user/papers directory
papers_path = fs.join(inspection.pts_user_dir, "papers")

# -----------------------------------------------------------------

class Papers(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Create the session
        self.session = requests.session()

        # Create the papers directory if necessary
        if not fs.is_directory(papers_path): fs.create_directory(papers_path)

        # Determine the path to the table
        self.table_path = fs.join(papers_path, "papers.dat")

        #if not fs.is_file(self.table_path):
            ## Create a new table
            #names = ["Label", "Bibcode", "Title", "Authors", "Journal", "Date", "URL"]
            #data = [[] for _ in names]


        #self.table = dict()

        if not fs.is_file(self.table_path):

            # Create a new dictionary
            self.table = dict()

            serialization.dump(self.table, self.table_path)

        else:

            self.table = serialization.load(self.table_path)

    # -----------------------------------------------------------------

    def add_entry(self, label, url):

        """
        This function ...
        :param label:
        :param url:
        :return:
        """

        r = self.session.get(url)

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        refereed_link = None
        arxiv_link = None

        for e in tree.iter():

            if not e.tag == 'a': continue

            if e.text_content() == "Full Refereed Journal Article (PDF/Postscript)":

                refereed_link = e.get("href")

            elif e.text_content() == "arXiv e-print":

                arxiv_link = e.get("href")

        print("ref", refereed_link)
        print("arxiv", arxiv_link)

        table_list = [e for e in tree.iter() if e.tag == 'table']
        table = table_list[1]

        table_rows = [e for e in table.iter() if e.tag == 'tr']

        title = None
        authors = []
        journal = None
        date = None
        bibcode = None

        for row in table_rows:

            property = None

            for e in row.iter():
                if e.tag == "b":
                    property = e.text_content().split(":")[0]
                    break

            if property == "Title": title = row.getchildren()[2].text_content()
            elif property == "Authors":

                for ee in row.getchildren()[2].iter():

                    if ee.tag != "a": continue

                    authors.append(ee.text.encode('ascii','ignore'))

                    #authors.append(ee.values())

                    #print(ee.body)

            elif property == "Publication": journal = row.getchildren()[2].text_content()
            elif property == "Publication Date": date = row.getchildren()[2].text_content().split(",")[0]
            elif property == "Bibliographic Code": bibcode = row.getchildren()[2].text_content()

        local_path = fs.join(papers_path, bibcode + ".pdf")

        # Download the PDF
        if refereed_link is not None:
            urllib.urlretrieve(refereed_link, local_path)
        else: urllib.urlretrieve(arxiv_link, local_path)

        entry = dict()
        entry["title"] = title
        entry["Authors"] = authors
        entry["Journal"] = journal
        entry["Date"] = date
        entry["Bibcode"] = bibcode

        self.table[label] = entry

    # -----------------------------------------------------------------

    def path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        bibcode = self.table[label]["Bibcode"]

        path = fs.join(papers_path, bibcode + ".pdf")

        return path

    # -----------------------------------------------------------------

    def open(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        path = self.path(label)

        subprocess.call(["open", path])

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        serialization.dump(self.table, self.table_path)

    # -----------------------------------------------------------------

    def to_bibliography(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        pass

# -----------------------------------------------------------------

