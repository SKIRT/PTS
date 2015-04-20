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
import sys, string, os
import re, time
from types import *
from httplib import HTTP

# Import astronomical modules
import pyregion

# Import relevant PTS modules
from pts.log import Log

# -----------------------------------------------------------------

# Useful adresses
domain = "www.nofs.navy.mil"      # USNO web site
script = "/cgi-bin/tfch3tI.cgi"   #ajp: USNO updates the script occasionally. Check by running search on web-page, then view source info.

# Create a logger
log = Log()

# -----------------------------------------------------------------

## This function downloads the source code of a web page and returns its contents as a list of lines
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

## This function converts right ascensions or declinations from degrees to hh::mm::ss or dd::mm::ss, respectively
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

## This function fetches data from the USNOB or NOMAD catalog in a box centered around a specified location (right
#  ascension, declination) and with a certain region. It takes the following parameters:
#
#  - center: the position (right ascension, declination) of the center of the box (in decimal degrees)
#  - radius: the radius for the (square) box
#  - catalog: defines which catalog to use (default = nomad)
#  - filename: gives the name of the ouput reference catalog (default = asc_ref.txt)
#  - magnitude: gives alternate USNO or NOMAD sort magnitude (default = R2)
#  - bright: sets bright search magnitude (default = 5.0)
#  - faint: sets faint search magnitude (default = 19.0)
#  - epoch: gives the desired input coordinate epoch (default = 2000.0)
#  - retry: sets the number of retries (default = 15)
#  - debug: set to True for debug mode (default = False)
#  - qonly: set to True if you only want to get the full URL for the query (default = False)
#  - text: set to True if the ID's of the stars have to be included as text in the region (default = False)
#
def fetch(center, radius=20, catalog="nomad", filename="asc_ref.txt", magnitude="R2", bright=5.0, faint=19.0, epoch="2000.0", retry=15, debug=False, qonly=False, text=False):

    # Get the RA and DEC
    rain = center[0]
    dcin = center[1]

    # Convert the right ascension to hours, minutes, seconds
    (rah,ram,ras)=dec2sxg(rain, ra=1)

    raout='+'.join([rah,ram,ras])
    ranice=':'.join([rah,ram,ras])

    # Convert the declination to degrees, minutes, seconds
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
          (script,raout,dcout,epoch,radius,radius,catalog,magnitude,bright,faint)

    log.info("Sending %s request for %s, %s (%s')" % (catalog,ranice,dcnice,radius))

    if qonly:

        log.info("Full URL of initial query: " + url)
        return

    # Download the results from the web
    results = download_web_page(domain, url)

    if debug:

        fil=open("debug.html", 'w')
        fil.write(results)

    # Check if something went wrong when getting the results
    if results is None or len(results) < 100:

        log.error("Something went wrong, results too short")
        exit()

    # Successful first request:  Get link for "progress" web page
    urltrk=""
    inrange=0
    lines=string.split(results,'\n')
    for line in lines:
        re1=re.search("refresh.+(http:.+\.html)",line)
        if re1:
            urltrk=re1.group(1)
            break

    if len(urltrk) > 1:

        urlcat=urltrk.replace("fch.html", catalog)

        log.info("Tracking results through %s" % urltrk)
        log.info("Catalog should appear as %s" % urlcat)

    else:

        log.error("Trouble identifying tracking URL")
        exit()

    re2=re.search("http://([^/]+)(/.+)$",urlcat)
    if re2:
        dom2=re2.group(1)
        url2=re2.group(2)
    else:

        log.error("Problem parsing target URL")
        exit()

    # Wait for successful processing of query
    time.sleep(10)

    # Start trying to retrieve the catalog
    itry=0
    while True:

        # Download the catalog
        results2 = download_web_page(dom2, url2)

        if results2:

            lines=string.split(results2,'\n')
            if len(lines)>=26: break

        itry +=1
        if itry > retry:

            log.error("Giving up")
            return

        # Wait before next try
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

        # Calculate a radius for this star
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

        # Add the parameters of this star to the region string
        region_string += regline

    # Create a region
    region = pyregion.parse(region_string)

    # Return the region
    return region
