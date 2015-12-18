#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

## \package pts.magic.tools.sesame Contains the Sesame class.

# -----------------------------------------------------------------

"""
Created on Mar 13, 2011
Sesame class to access Sesame name resolver service
Based on 2005-06-11 by Shui Hung Kwok 
See http://cdsweb.u-strasbg.fr/doc/sesame.htx for description of Sesame
@author: shkwok
"""

from urllib2 import urlopen
#from xparser.XParser import XParser
#from .. import XParser

# -----------------------------------------------------------------

class Sesame (object):

    """
    This class ...
    """

    CatalogOpt = "SNV" # S simbad, N ned, V vizier, A All
    OutputOpt = "oxp" # xp for xml as text/plain rather then text/xml (-ox)
    SesameURL = "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame"
    
    def __init__(self, urn=SesameURL, opt=CatalogOpt, opt1=OutputOpt):
        
        """
        Initializes Sesame URL and options
        Default options are SNV for CatalogOpt
        and -oxp for OutputOpt.
        SNV = Simbad + NED + Vizier and A for All
        The order indicates the order to search.
        "All" means search all services, otherwise stops when 
        first entry found.
        Output options start with -o 
        followed by 
        x : xml output
        p : plain text
        I : include all identifiers
        """
        
        self.catOpt = opt
        self.outOpt = opt1
        self.urn = urn
        # Sesame

    def getCoord(self, node):
        
        """
        Helper method to extract ra and dec from node 
        """
        
        res = node.getResource("/Sesame/Target");
        resolvers = res.getChildren ("Resolver")
        for r in resolvers:
            try:
                ra = float (r.getResourceContent("/Resolver/jradeg").strip())
                dec = float (r.getResourceContent("/Resolver/jdedeg").strip())
                return ra, dec
            except Exception:
                raise Exception, "invalid coordinates"
        else:
            raise Exception, "no ra/dec values found"
        # getCoord

    def getAliases(self):
        
        """
        Extracts aliases for the given target.
        Returns a list of names.
        """
        
        res = []
        for resolver in self.xml.root.Sesame.Resolver:
            try:
                for a in resolver.alias:
                    res.append (a.content)
            except:
                pass
        return res
  
    def buildQuery(self, name, all=True):
        
        """
        Builds query URL for use with HTTP GET
        If all is true, then all known identifiers shall be returned.
        """
        
        opt = self.catOpt
        opt1 = '-' + self.outOpt
        if all:
            opt += 'A'
            opt1 += 'I' # all identifiers
        queryURL = "%s/%s/%s?%s" % (self.urn, opt1, opt, name)
        return queryURL

    def resolveRaw(self, name, all=True):
        
        """
        Performs a raw query.
        Returns what the server returns.
        """
        
        query = self.buildQuery (name, all)
        print "query=", query
        hcon = urlopen (query)
        res = hcon.read ()
        hcon.close ()
        return res
    
    def resolve(self, name, all=True):
        
        """
        Performs a query.
        Returns ra and dec 
        """
        
        query = self.buildQuery(name, all)
        xp = XParser()
        xn = xp.parseFromFile(query)
        return self.getCoord(xn)