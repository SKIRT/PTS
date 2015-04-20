#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.starcatalog This module is used to fetch data from the USNOB/NOMAD star catalogues. It is adapted
#  from the 'getusnobn.py' script (version 1.1 2004/06/11) found at
#  https://github.com/cenko/python/blob/master/catalogs/getusnobn.py, written by D. Fox and A.J. Pickles and in
#  turn adapted from 'cda.py' (A. Ptak) and based on code developed for XAssist, see http://xassist.pha.jhu.edu

# -----------------------------------------------------------------

# Import standard modules
import sys, string, os, traceback
import getopt, re, time, glob
from types import *
from httplib import HTTP
import numpy as np

# Import astronomical modules
import pyregion

# -----------------------------------------------------------------

# Defaults / Globals
domain = "www.nofs.navy.mil" # USNO web site
#ajp: USNO updates the script occasionally. 
#     Check by running search on web-page, then view source info.
script = "/cgi-bin/tfch3tI.cgi"

# -----------------------------------------------------------------

class dummy_ui:
    def error_mesg(self, txt):
        print txt

    def writeln(self, txt):
        print txt

# -----------------------------------------------------------------
    
def download_web_page(domain, url):
    try:
        h = HTTP(domain)
        h.putrequest("GET", url)
        h.putheader('Accept', 'text/html')
        h.putheader('Accept', 'text/plain')
        h.endheaders()
    except:
        return(None)

    try:
        errcode, errmsg, headers = h.getreply()
    except:
        sys.stderr.write("Error in receiving response from " + domain + \
                         '\n')
        return None
    if errcode != 200:
        sys.stderr.write("Error in receiving response from " + domain + \
                         '\n')
        return None
    results = h.getfile().read()
    return(results)

# -----------------------------------------------------------------

def dec2sxg(coord,ra=0):

    coox=float(coord)

    if ra:
        rah=int(coox/15.0)
        ram=int(60*(coox/15.0-rah))
        ras=60*(60*(coox/15.0-rah)-ram)
        return ("%02d" % rah,"%02d" % ram,"%07.4f" % ras)

    dcsgn= '+' if coord > 0 else '-'

    absdec=abs(coox)
    dcd=int(absdec)
    dcm=int(60*(absdec-abs(dcd)))
    dcs=60*(60*(absdec-abs(dcd))-dcm)

    return (dcsgn+"%02d" % dcd,"%02d" % dcm,"%06.3f" % dcs)

# -----------------------------------------------------------------

## This function ...
# "Usage: %s [-c catusno] [-o orefcat] [-m magclr] [-b brtmag] [-f fntmag] [-e epoch] [-s flag] [-r retry] [-g regfile] <ra> <dec> [size]" % xname
#    print "    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg"
#    print "    <size> is length of a side in arcmin (%s)" % def_radius
#    print "    -c <usnob/nomad> gives alternate input catalog (%s)" % def_catusno
#    print "    -o <orefcat> gives alternate name for ouput reference catalog (%s)" % def_orefcat
#    print "    -m <R2/K> gives alternate USNO or NOMAD sort magnitude (%s)" % def_magclr
#    print "    -b <val> sets bright search magnitude (%s)" % def_bright
#    print "    -f <val> sets  faint search magnitude (%s)" % def_faint
#    print "    -e <epoch> gives the desired input coordinate epoch (%s)" % def_epoch
#    print "    -s <off/on> sets the silent output flag (%s)" % def_silent
#    print "    -r <val> sets the number of retries (%s)" % MAXTRY
#    print "    -g <regfile> gives alternate name for ouput region file (%s)" % def_orefcat.replace('.cat','.reg')
#
def fetch(center, radius=20, catalog="nomad", filename="asc_ref.txt", magnitude="R2", bright=5.0, faint=19.0, epoch="2000.0", retry=15, debug=False, qonly=False, text=False):

    # User-controllable values
    catusno = catalog
    magclr = magnitude
    bright = bright
    faint = faint
    epoch = epoch
    radius = radius
    retry = retry

    # Process details
    (xdir,xname)=os.path.split(sys.argv[0])
    pid=os.getpid()

    # Get the RA and DEC
    rain = center[0]
    dcin = center[1]

    # Parse coordinates (RA)
    (rah,ram,ras)=dec2sxg(rain, ra=1)

    raout='+'.join([rah,ram,ras])
    ranice=':'.join([rah,ram,ras])

    # Parse coordinates (Dec)
    (dcd,dcm,dcs)=dec2sxg(dcin, ra=0)

    dcout='+'.join([dcd,dcm,dcs])
    dcnice=':'.join([dcd,dcm,dcs])
        
    # First request: (colbits=cb_flg&) removed below, as hex output harder to read back in
    url = ("%s?ra=%s&dec=%s&equinox=J2000&epoch=%s&cextract=rect&" + \
           "rawid=%s&decwid=%s&wunits=Minutes&cat=%s&surims=None&" + \
           "getcat=yes&colbits=cb_id&colbits=cb_ra&slf=ddd.ddd/dd.ddd&" + \
           "colbits=cb_sigra&colbits=cb_mura&colbits=cb_smural&colbits=cb_fitpts&" + \
           "colbits=cb_mag&clr=%s&skey=mag&bri=%s&fai=%s&" + \
           "gzf=No&cftype=ASCII") % \
          (script,raout,dcout,epoch,radius,radius,catusno,magclr,bright,faint)

    print "Sending %s request for %s, %s (%s')" % (catusno,ranice,dcnice,radius)

    if qonly:
        print "Full URL of initial query:"
        print url
        return
    
    results = download_web_page(domain, url)

    if debug:
        test1="%s-%d-1.html" % (xname,pid)
        fil=open(test1, 'w')
        fil.write(results)

    if type(results)==NoneType or len(results)<100:
        sys.exit("Something went wrong, results too short\n")

    # Successful first request:  Get link for "progress" web page
    urltrk=""
    inrange=0
    lines=string.split(results,'\n')
    for line in lines:
        re1=re.search("refresh.+(http:.+\.html)",line)
        if re1:
            urltrk=re1.group(1)
            break

    if len(urltrk)>1:
        urlcat=urltrk.replace("fch.html",catusno)
        print "Tracking results through %s" % urltrk
        print "Catalog should appear as %s" % urlcat
    else:
        sys.exit("Trouble identifying tracking URL\n")

    re2=re.search("http://([^/]+)(/.+)$",urlcat)
    if re2:
        dom2=re2.group(1)
        url2=re2.group(2)
    else:
        sys.exit("Problem parsing target URL\n")

    # Wait a bit
    time.sleep(10)

    # Start trying to retrieve the catalog
    itry=0
    while 1:
        results2 = download_web_page(dom2, url2)

        if results2:
            lines=string.split(results2,'\n')
            if len(lines)>=26:
                break

        itry +=1
        if itry > retry:
            print "Giving up"
            return

        time.sleep(10)

    # Add comment flags to header
    for i in range(len(lines)):
        if re.search('^\s*\d\d\d\d',lines[i]):
            break
        lines[i]='#'+lines[i]
    nhead=i

    # Write the catalog to a text file, if requested
    if filename:
        fcat=open(filename, 'w')
        for line in lines:
            # Remove leading spaces
            while line.startswith(' '):
                line=line[1:]
            fcat.write(line+'\n')
        fcat.close()

    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"

    # For each line that represents a star
    for line in lines[nhead:]:

        # Split the information in this line
        els = line.split()

        # If the line is invalid, go to the next one
        if len(els) < 16: continue

        rstar=5.0*(20.0/float(els[16]))**2
        if text:
            regline=("fk5;circle(%s,%s,%.2f\") # text={%s} " +
                     "B1MAG={%s} R1MAG={%s} B2MAG={%s} R2MAG={%s} " +
                     "I2MAG={%s}\n") % (els[1],els[2],rstar,els[0],
                                        els[12],els[13],els[14],els[15],els[16])
        else:
            regline=("fk5;circle(%s,%s,%.2f\") # " +
                     "B1MAG={%s} R1MAG={%s} B2MAG={%s} R2MAG={%s} " +
                     "I2MAG={%s}\n") % (els[1],els[2],rstar,
                                        els[12],els[13],els[14],els[15],els[16])
        region_string += regline

    # Create a region
    region = pyregion.parse(region_string)

    # Return the region
    return region
