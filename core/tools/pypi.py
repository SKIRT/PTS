#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
pypi_cli
~~~~~~~~

A command line interface to the Python Package Index.

:copyright: (c) 2014 by Steven Loria.
:license: MIT, see LICENSE for more details.
"""

#Copyright 2014 Steven Loria

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import division, print_function

# Import standard modules
import re
from collections import OrderedDict

# Import other modules
#PY2 = int(sys.version[0]) == 2
#if PY2:
from xmlrpclib import ServerProxy
#from urllib import quote as urlquote # eliminated
#else:
#    from xmlrpc.client import ServerProxy
#    from urllib.parse import quote as urlquote

import requests
#from dateutil.parser import parse as dateparse # eliminated
#import click
#from click import echo, style, echo_via_pager
#from click.termui import get_terminal_size

__version__ = '0.4.1'
__author__ = 'Steven Loria'
__license__ = 'MIT'

DATE_FORMAT = "%y/%m/%d"
MARGIN = 3
DEFAULT_SEARCH_RESULTS = 100

TICK = '*'
DEFAULT_PYPI = 'https://pypi.python.org/pypi'
PYPI_RE = re.compile('''^(?:(?P<pypi>https?://[^/]+/pypi)/)?
                        (?P<name>[-A-Za-z0-9_.]+)
                        (?:/(?P<version>[-A-Za-z0-9.]+))?$''', re.X)
SEARCH_URL = 'https://pypi.python.org/pypi?%3Aaction=search&term={query}'

# Number of characters added by bold formatting
_BOLD_LEN = 8
# Number of characters added by color formatting
_COLOR_LEN = 9

# -----------------------------------------------------------------

def get_package(name_or_url, client=None):

    m = PYPI_RE.match(name_or_url)
    if not m:
        return None
    pypi_url = m.group('pypi') or DEFAULT_PYPI
    name = m.group('name')
    return Package(name, pypi_url=pypi_url, client=client)

# -----------------------------------------------------------------

#def echo_download_summary(package):
#    echo('Last day:    {daily:12,}'.format(daily=package.downloads_last_day))
#    echo('Last week:   {weekly:12,}'.format(weekly=package.downloads_last_week))
#    echo('Last month:  {monthly:12,}'.format(monthly=package.downloads_last_month))

# -----------------------------------------------------------------

def browse(package, homepage):

    """
    Browse to a package's PyPI or project homepage.
    #@cli.command()
    #@click.option('--homepage', is_flag=True, default=False)
    #@click.argument('package', required=True)
    """

    p = Package(package)
    try:
        if homepage:
            #echof(u'Opening homepage for "{0}"...'.format(package), bold=True)
            print('Opening homepage for "{0}"...'.format(package))
            url = p.home_page
        else:
            #echof(u'Opening PyPI page for "{0}"...'.format(package), bold=True)
            print('Opening PyPI page for "{0}"...'.format(package))
            url = p.package_url
    except NotFoundError:
        #abort_not_found(package)
        raise Exception("Not found")
    #click.launch(url)

    print(url)

# -----------------------------------------------------------------

def search(query, n_results=100):

    """Search for a pypi package.

    #@cli.command()
    #@click.option('--web', '-w', is_flag=True, default=False,
    #    help='Open search results in your web browser.')
    #@click.option('--n-results', '-n', default=DEFAULT_SEARCH_RESULTS,
    #    help='Max number of results to show.')
    #@click.argument('query', required=True, type=str)

    \b
    Examples:
        \b
        pypi search requests
        pypi search 'requests oauth'
        pypi search requests -n 20
        pypi search 'requests toolbelt' --web

    """

    searcher = Searcher()
    results = searcher.search(query, n=n_results)
    return results

# -----------------------------------------------------------------

def info(package, long_description, classifiers, license):

    """
    Get info about a package or packages.

    #@cli.command()
    #@click.option('--license/--no-license',
    #    is_flag=True, default=True, help='Show license.')
    #@click.option('--classifiers', '-c',
    #    is_flag=True, default=False, help='Show classifiers.')
    #@click.option('--long-description', '-L',
    #    is_flag=True, default=False, help='Show long description.')
    #@click.argument('package', nargs=-1, required=True)
    """

    client = requests.Session()
    for name_or_url in package:
        package = get_package(name_or_url, client)
        if not package:
            #echo(style(
            #    u'Invalid name or URL: "{name}"'.format(name=name_or_url), fg='red'),
            #    file=sys.stderr)
            print('Invalid name or URL: "{name}"'.format(name=name_or_url))
            continue

        # Name and summary
        try:
            info = package.data['info']
        except NotFoundError:
            #echo(style(u'No versions found for "{0}". Skipping. . .'.format(package.name),
            #    fg='red'), file=sys.stderr)
            print('No versions found for "{0}". Skipping. . .'.format(package.name))
            continue
        #echo_header(name_or_url)
        if package.summary:
            #echo(package.summary)
            print(package.summary)

        # Version info
        #echo()
        print()
        #echo('Latest release:   {version:12}'.format(version=info['version']))
        print('Latest release:   {version:12}'.format(version=info['version']))

        # Long description
        if long_description:
            #echo()
            print()
            #echo(package.description)
            print(package.description)

        # Download info
        #echo()
        print()
        echo_download_summary(package)

        # Author info
        #echo()
        print()
        author, author_email = package.author, package.author_email
        if author:
            #echo(u'Author:   {author:12}'.format(**locals()))
            print('Author:   {author:12}'.format(**locals()))
        if author_email:
            #echo(u'Author email: {author_email:12}'.format(**locals()))
            print('Author email: {author_email:12}'.format(**locals()))

        # Maintainer info
        maintainer, maintainer_email = package.maintainer, package.maintainer_email
        if maintainer or maintainer_email:
            #echo()
            print()
        if maintainer:
            #echo(u'Maintainer:   {maintainer:12}'.format(**locals()))
            print('Maintainer:   {maintainer:12}'.format(**locals()))
        if maintainer_email:
            #echo(u'Maintainer email: {maintainer_email:12}'.format(**locals()))
            print('Maintainer email: {maintainer_email:12}'.format(**locals()))

        # URLS
        #echo()
        print()
        #echo(u'PyPI URL:  {pypi_url:12}'.format(pypi_url=package.package_url))
        print('PyPI URL:  {pypi_url:12}'.format(pypi_url=package.package_url))
        if package.home_page:
            #echo(u'Home Page: {home_page:12}'.format(home_page=package.home_page))
            print('Home Page: {home_page:12}'.format(home_page=package.home_page))
        if package.docs_url:
            #echo(u'Documentation: {docs_url:12}'.format(docs_url=package.docs_url))
            print('Documentation: {docs_url:12}'.format(docs_url=package.docs_url))

        # Classifiers
        if classifiers:
            #echo()
            print()
            #echo(u'Classifiers: ')
            print("Classifiers: ")
            for each in info.get('classifiers', []):
                #echo('\t' + each)
                print('\t' + each)

        if license and package.license:
            #echo()
            print()
            #echo(u'License: ', nl=False)
            print('License: ')
            # license may be just a name, e.g. 'BSD' or the full license text
            # If a new line is found in the text, print a new line
            if package.license.find('\n') >= 0 or len(package.license) > 80:
                #echo()
                print()
            #echo(package.license)
            print(package.license)
        #echo()
        print()


# Utilities
# #########

# -----------------------------------------------------------------

def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated."""
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazy_property

# -----------------------------------------------------------------

class PackageError(Exception):
    pass

# -----------------------------------------------------------------

class NotFoundError(PackageError):
    pass

# -----------------------------------------------------------------

# API Wrapper
# ###########

class Package(object):

    """
    This function ...
    """

    def __init__(self, name, client=None, pypi_url=DEFAULT_PYPI):

        self.client = client or requests.Session()
        self.name = name
        self.url = '{pypi_url}/{name}/json'.format(pypi_url=pypi_url, name=name)

    # -----------------------------------------------------------------

    @lazy_property
    def data(self):

        resp = self.client.get(self.url)
        if resp.status_code == 404:
            raise NotFoundError('Package not found')
        return resp.json()

    # -----------------------------------------------------------------

    @lazy_property
    def versions(self):

        """Return a list of versions, sorted by release datae."""

        return [k for k, v in self.release_info]

    # -----------------------------------------------------------------

    @lazy_property
    def version_downloads(self):

        """Return a dictionary of version:download_count pairs."""

        ret = OrderedDict()
        for release, info in self.release_info:
            download_count = sum(file_['downloads'] for file_ in info)
            ret[release] = download_count
        return ret

    # -----------------------------------------------------------------

    @property
    def release_info(self):

        release_info = self.data['releases']
        # filter out any versions that have no releases
        filtered = [(ver, releases) for ver, releases in release_info.items()
                    if len(releases) > 0]
        # sort by first upload date of each release
        return sorted(filtered, key=lambda x: x[1][0]['upload_time'])

    # -----------------------------------------------------------------

    @lazy_property
    def downloads(self):

        """Total download count.
        :return: A tuple of the form (version, n_downloads)
        """

        return sum(self.version_downloads.values())

    # -----------------------------------------------------------------

    @lazy_property
    def max_version(self):

        """Version with the most downloads.
        :return: A tuple of the form (version, n_downloads)
        """

        data = self.version_downloads
        if not data:
            return None, 0
        return max(data.items(), key=lambda item: item[1])

    # -----------------------------------------------------------------

    @lazy_property
    def min_version(self):

        """Version with the fewest downloads."""

        data = self.version_downloads
        if not data:
            return (None, 0)
        return min(data.items(), key=lambda item: item[1])

    # -----------------------------------------------------------------

    @lazy_property
    def average_downloads(self):

        """Average number of downloads."""

        return int(self.downloads / len(self.versions))

    # -----------------------------------------------------------------

    @property
    def author(self):
        return self.data['info'].get('author')

    # -----------------------------------------------------------------

    @property
    def description(self):
        return self.data['info'].get('description')

    # -----------------------------------------------------------------

    @property
    def summary(self):
        return self.data['info'].get('summary')

    # -----------------------------------------------------------------

    @property
    def author_email(self):
        return self.data['info'].get('author_email')

    # -----------------------------------------------------------------

    @property
    def maintainer(self):
        return self.data['info'].get('maintainer')

    # -----------------------------------------------------------------

    @property
    def maintainer_email(self):
        return self.data['info'].get('maintainer_email')

    # -----------------------------------------------------------------

    @property
    def license(self):
        return self.data['info'].get('license')

    # -----------------------------------------------------------------

    @property
    def downloads_last_day(self):
        return self.data['info']['downloads']['last_day']

    # -----------------------------------------------------------------

    @property
    def downloads_last_week(self):
        return self.data['info']['downloads']['last_week']

    # -----------------------------------------------------------------

    @property
    def downloads_last_month(self):
        return self.data['info']['downloads']['last_month']

    # -----------------------------------------------------------------

    @property
    def package_url(self):
        return self.data['info']['package_url']

    # -----------------------------------------------------------------

    @property
    def home_page(self):
        return self.data['info'].get('home_page')

    # -----------------------------------------------------------------

    @property
    def docs_url(self):
        return self.data['info'].get('docs_url')

    # -----------------------------------------------------------------

    def __repr__(self):
        return '<Package(name={0!r})>'.format(self.name)

# -----------------------------------------------------------------

class Searcher(object):

    """PyPI package search wrapper that uses the PyPI's XMLRPC API.
    Search algorithm adapted from Supreet Sethi's implementation (MIT Licensed).
    https://github.com/djinn/pypi-json/blob/master/LICENSE.md
    """

    STOP_WORDS = {
        "a", "and", "are", "as", "at", "be", "but", "by",
        "for", "if", "in", "into", "is", "it",
        "no", "not", "of", "on", "or", "such",
        "that", "the", "their", "then", "there", "these",
        "they", "this", "to", "was", "will",}

    NAME_MATCH_WEIGHT = 16
    CONTAINS_NAME_MULT = 4
    NAME_IN_SUMMARY_MULT = 2

    # -----------------------------------------------------------------

    def __init__(self, pypi_url=DEFAULT_PYPI, client=None):

        self.pypi_url = pypi_url
        self.client = client or ServerProxy(pypi_url)

    # -----------------------------------------------------------------

    def score(self, tokens, record):

        score = 0
        name, summary = record['name'].lower(), record['summary']
        for token in tokens:
            qtf = 0
            if token == name:
                qtf += self.NAME_MATCH_WEIGHT
            else:
                n_name_matches = len(re.compile(token).findall(name))
                qtf += self.CONTAINS_NAME_MULT * n_name_matches
            if record['summary'] is not None:
                summary_matches = re.compile(token).findall(summary.lower())
                qtf += 2 * len(summary_matches)
            score += qtf
        return score

    # -----------------------------------------------------------------

    def search(self, query, n=None):

        tokens = [each.strip() for each in query.strip().lower().split()
                  if each not in self.STOP_WORDS]
        results = self.client.search({'name': tokens}, 'and')
        visited = []
        nd = []
        for result in results:
            name = result['name']
            try:
                visited.index(name)
            except ValueError:
                nd.append(result)
                visited.append(name)
        ranked = [(self.score(tokens, result), result) for result in nd]
        sorted_results = sorted(ranked, reverse=True, key=lambda t: t[0])
        limited_results = sorted_results[:n] if n else sorted_results
        return (result for score, result in limited_results)

# -----------------------------------------------------------------
