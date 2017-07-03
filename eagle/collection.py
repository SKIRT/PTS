#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.collection Utilities to load and use data from EAGLE SKIRT-run result collections.
#
# The facilities in this module allow loading one or more EAGLE SKIRT-run result collections,
# extracting the desired galaxy properties, and calculating observables from these properties.

# ----------------------------------------------------------------------

import pickle
import types
import os.path
import numpy as np

from . import config
from ..core.filter.broad import BroadBandFilter
from ..core.basics.greybody import Bnu, GreyBody, kappa350_Cortese, kappa350_Zubko

# ----------------------------------------------------------------------

## An instance of the Collection class represents the contents of a particular collection of EAGLE SKIRT-run results.
# It offers various facilities for querying galaxy properties.
class Collection:

    ## The constructor loads the contents of the specified collection and optionally reads extra data fields
    # from column text files.
    # The collection name should \em not include the directory (which is taken from eagle.config) nor the
    # postfix "_info_collection.dat". The optional collection label can be provided as a short identifier,
    # for example for use in the legend of a plot. If a list of file names is provided, the constructor adds
    # information from these text data files to the collection. The first column in each file contains a GalaxyID,
    # and subsequent columns contain the extra fields for that galaxy in the order specified by the field names.
    # The specified files are simply concatenated.
    def __init__(self, collection_name, collection_label=None, file_names=None, field_names=None):
        self._name = collection_name
        self._label = collection_label if collection_label is not None else collection_name

        # load the collection
        infilepath = os.path.join(config.collections_path, collection_name+"_info_collection.dat")
        infile = open(infilepath, "r")
        self._info = pickle.load(infile)
        infile.close()

        # construct an index on galayid
        self._ids = { }       # key = galaxy id, value = index in collection
        index = 0
        for galaxyid in self._info["galaxy_id"]:
            self._ids[galaxyid] = index
            index+=1

        if (file_names is not None) and (field_names is not None):
            # read the extra fields from the input files, and store them keyed on galaxy id
            galaxies = {}
            for filename in file_names:
                for row in np.loadtxt(filename):
                    galaxies[ int(row[0]) ] = row[1:]

            # create appropriate collection entries for the new fields
            for field in field_names:
                self._info[field] = np.zeros_like(self._info["skirt_run_id"])

            # copy the field values for each galaxy in the collection
            index = 0
            for galaxyid in self._info["galaxy_id"]:
                for field, value in zip(field_names, galaxies[galaxyid]):
                    self._info[field][index] = value
                index+=1

        # replace infinities (signifying a non-detection) by NaNs
        for value in self._info.values():
            value[np.isinf(value)] = np.nan

    ## This function returns the collection name.
    def name(self):
        return self._name

    ## This function returns the collection label.
    def label(self):
        return self._label

    ## This function returns a set containing all galaxy ids in the collection.
    def galaxy_ids(self):
        return set(self._ids.keys())

    ## This function returns a set containing the names of the properties provided in the collection.
    def property_names(self):
        return set(self._info.keys())

    ## This function returns an array with the values of the specified property for each galaxy in the specified
    # list of galaxy ids, in the same order. An error is raised if the galaxy id and/or the property are not available.
    def property_values(self, property_name, galaxy_ids):
        result = np.zeros(len(galaxy_ids))
        index = 0
        for galaxy_id in galaxy_ids:
            result[index] = self._info[property_name][self._ids[galaxy_id]]
            index+=1
        return result

    ## This function returns a dictionary with key-value pairs for all properties in the collection.
    # The key is the property name, and the value is a single-dimensional array with
    # the values of that property for all galaxies in the collection, in arbitrary order.
    def all_property_values(self):
        return self._info

# ----------------------------------------------------------------------

## An instance of the CollectionSet class manages the contents of one or more EAGLE SKIRT-run result collections.
# It allows querying the properties in all collections for the set of common galaxies.
class CollectionSet:

    ## The constructor loads the specified collections, optionally reading extra data fields from text files,
    # and prints some statistics. The arguments are similar to those for the Collection constructor, except
    # that the first two should be lists of equal length (if provided).
    def __init__(self, collection_names, collection_labels=None, file_names=None, field_names=None):
        # load the collections
        if isinstance(collection_names,types.StringTypes):
            collection_names = [ collection_names ]
        if collection_labels is None:
            collection_labels = [ None ] * len(collection_names)
        self._collections = [ Collection(name,label,file_names,field_names) \
                                for name,label in zip(collection_names,collection_labels) ]

        # find the set of common galaxies and the set of common properties
        if len(self._collections) > 1:
            self._ids = sorted(reduce(lambda x,y: x&y, [ c.galaxy_ids() for c in self._collections ]))
            self._props = reduce(lambda x,y: x&y, [ c.property_names() for c in self._collections ])
        else:
            self._ids = sorted(self._collections[0].galaxy_ids())
            self._props = self._collections[0].property_names()

        # print the number of common galaxies
        print "Loaded a set of {} collections with {} common galaxies and {} common properties" \
                    .format(len(self._collections), len(self._ids), len(self._props))

    ## This function returns a two-dimensional array with the values of the specified property for all common galaxies
    # in all collections of the set. The index on the first axis iterates over the collections, the index on the last
    # axis iterates over the galaxies, in order of increasing galaxy id.
    def property_values(self, property_name):
        return np.vstack([ c.property_values(property_name,self._ids) for c in self._collections ])

    ## This function returns a dictionary with key-value pairs for all properties that are common to the
    # collections in the set. The key is the property name, and the value is a two-dimensional array with
    # the values of that property for all common galaxies in all collections of the set. The index on the
    # first axis iterates over the collections, the index on the last axis iterates over the galaxies,
    # in order of increasing galaxy id.
    def all_property_values(self):
        return { p:self.property_values(p) for p in self._props }

# ----------------------------------------------------------------------

## An instance of the CollectionData class encapsulates the information contained in a Collection or in a
# CollectionSet instance, offering the following benefits:
#  - optionally include 'stripped' property names for a specific instrument (i.e. omitting the instrument name);
#  - access (common) properties for all (common) galaxies using Python property syntax;
#  - calculate a number of extra observable galaxy properties from the basic properties.
class CollectionData:

    ## The constructor retrieves the relevant information from the specified Collection or CollectionSet instance,
    # and adds stripped property names for the specified instrument, if any.
    def __init__(self, collection, instrument=None):
        # get the info dictionary
        self._info = collection.all_property_values()

        # add stripped instrument keys if requested
        if instrument is not None:
            for key in self._info.keys():
                if key.startswith("instr_"):
                    segments = key.split("_")
                    if segments[1]==instrument:
                        segments.pop(1)
                        cleankey = "_".join(segments)
                        self._info[cleankey] = self._info[key]

    ## This function ensures that the data array for property "some_property" in CollectionData instance cd
    # can be accessed through the regular Python cd.some_property syntax. For a single collection, the function
    # returns a single-dimensional array with the values of that property for all galaxies in the collection.
    # For multiple collections, the function returns a two-dimensional array with the values of that property
    # for all common galaxies in all collections of the set. The index on the first axis iterates over the
    # collections, the index on the last axis iterates over the galaxies.
    def __getattr__(self, name):
        return self._info[name]

    ## This function returns log10 of stellar mass (in solar units) according to Zibetti et al 2009, table B1,
    # using color g-i and i-band luminosity
    def log_stellar_mass_as_zibetti(self):
        color = self.instr_magnitude_sdss_g - self.instr_magnitude_sdss_i   # AB color g - i
        logUpsi = -0.963 + 1.032*color    # Upsilon in solar units (coefficients a_i and b_i for color g-i in table B1)
        logLi = (4.54 - self.instr_magnitude_sdss_i) / 2.5   # solar AB magnitude in i band is 4.54
        return logUpsi + logLi

    ## This function returns dust temperature (in K) and mass (in Msun) for best fit with Herschel 160, 250, 350, 500 data points
    # of the specified flux type (default is 'limited'), using beta=2 and kappa=kappa350_Cortese
    def dust_temperature_and_mass_from_grey_body_fit(self, fluxtype='limited'):
        # get the Herschel 160, 250, 350, 500 wavelengths
        waves = np.array( [ BroadBandFilter(fs).pivotwavelength() for fs in ("Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW")] )
        sigmas = np.array(( 3,1,1,3 ))      # pacs is less sensitive; longer wavelength fluxes are harder to measure

        # get the Herschel 160, 250, 350, 500 datapoints
        fluxstring = '''[ self.instr_fluxdensity_pacs_red_{0}, self.instr_fluxdensity_spire_psw_{0},
                          self.instr_fluxdensity_spire_pmw_{0}, self.instr_fluxdensity_spire_plw_{0} ]'''.format(fluxtype)
        fluxes = eval(fluxstring)

        # setup an iterator over the galaxies, specifying two to-be-allocated output arrays for T and M
        it = np.nditer([None, None, self.setup_distance_instrument] + fluxes,
                       op_flags = [['writeonly','allocate'],['writeonly','allocate'],['readonly'],
                                   ['readonly'], ['readonly'], ['readonly'], ['readonly']])

        # do the fit, iterating over the galaxies
        for Ti,Mi,di,f160i,f250i,f350i,f500i in it:
            greybody = GreyBody(di, 2, kappa350_Cortese)
            #greybody = GreyBody(di, 2, kappa350_Zubko)
            it[0],it[1] = greybody.fit(waves, (f160i,f250i,f350i,f500i), sigmas)

        # return the two result arrays T and M allocated by the iterator
        return it.operands[0:2]

    ## This function returns dust temperature (in K) for best fit with Herschel 160, 250, 350, 500 data points
    def dust_temperature_from_grey_body_fit(self, fluxtype='limited'):
        return self.dust_temperature_and_mass_from_grey_body_fit(fluxtype)[0]

    ## This function returns log10 of dust mass (in Msun) for best fit with Herschel 160, 250, 350, 500 data points
    def log_dust_mass_from_grey_body_fit(self, fluxtype='limited'):
        return log_if_positive(self.dust_temperature_and_mass_from_grey_body_fit(fluxtype)[1])

    ## This function returns fraction of total observed dust mass contributed by HII regions (in range 0..1),
    # calculated from the continuum fluxes through a best fit with Herschel 160, 250, 350, 500 data points
    def dust_fraction_in_hii_regions(self):
        T,Mhii = self.dust_temperature_and_mass_from_grey_body_fit("hii_continuum")
        T,Mother = self.dust_temperature_and_mass_from_grey_body_fit("other_continuum")
        return divide_if_positive(Mhii, Mhii+Mother)

    ## This function returns log10 of dust mass (in Msun) according to Cortese et al 2012, appendix B, using beta=2 for extended sources
    def log_dust_mass_as_cortese(self):
        x = log_divide_if_positive(self.instr_fluxdensity_spire_psw_limited,self.instr_fluxdensity_spire_plw_limited)
        logMFD = 16.880 - 1.559*x + 0.160*x**2 - 0.079*x**3 - 0.363*x**4
        logD = np.log10(self.setup_distance_instrument/1e6)
        logF = log_if_positive(self.instr_fluxdensity_spire_pmw_limited)
        logDust = logMFD + 2*logD + logF - 11.32
        #logDust += np.log10( kappa350_Cortese / kappa350_Zubko )    # compensate for kappa assumed in Cortese vs Zubko
        return logDust

    ## This function returns log10 of dust mass (in Msun) based on the dust temperature probed in the dust grid and the 350 micron flux
    def log_dust_mass_from_grid_temperature(self):
        Msun = 1.9891e30     # solar mass in kg
        pc = 3.08567758e16   # parsec in m
        f350 = self.instr_fluxdensity_spire_pmw_limited * 1e-26      # W/m2
        D = self.setup_distance_instrument * pc                      # m
        T = self.probe_average_temperature_dust                      # K
        T[T<1] = 1
        return log_divide_if_positive(f350*D*D, kappa350_Cortese * Bnu(350,T) * Msun)

# ----------------------------------------------------------------------

## This function returns a CollectionData instance encapsulating the common properties of the common galaxies
# for the collections with the specified names, optionally including extra data fields from column text files,
# and providing stripped property names for a specific instrument.
# It is simply a convenience function combing the CollectionSet and CollectionData constructors.
def load_collections(collection_names, file_names=None, field_names=None, instrument=None):
    return CollectionData(CollectionSet(collection_names, file_names=file_names, field_names=field_names), instrument)

# ----------------------------------------------------------------------

# Some generic functions used for calculating observables

## This function returns log10(x), or NaN for x<=0
def log_if_positive(x):
    positive = x>0
    result = np.empty_like(x)
    result[positive] = np.log10(x[positive])
    result[~positive] = np.nan
    return result

## This function returns x/y, or NaN for y<=0
def divide_if_positive(x,y):
    positive = y>0
    result = np.empty_like(x)
    result[positive] = x[positive] / y[positive]
    result[~positive] = np.nan
    return result

## This function returns log10(x/y), or NaN for x<=0 or y<=0
def log_divide_if_positive(x,y):
    result = np.zeros_like(x)
    positive = y>0
    result[positive] = x[positive] / y[positive]
    positive = result>0
    result[positive] = np.log10(result[positive])
    result[~positive] = np.nan
    return result

# ----------------------------------------------------------------------
