#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.wcsfinder Contains the WCSFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class WCSFinder(Configurable):

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
        super(WCSFinder, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

# From: https://github.com/fred3m/astropyp/blob/master/astropyp/astrometry.py

import numpy as np
from collections import OrderedDict
import logging

from astropy.io import fits
import astropy.wcs
from astropy import coordinates
import astropy.units as apu
from astropy import table

logger = logging.getLogger('astropyp.astrometry')


class Astrometry:
    """
    ImageSolution can be used to calculate multivariate polynomial
    transformations from one coordinate to another
    """

    def __init__(self, transforms=None,
                 order=None, init_mean=[None, None], init_stddev=[None, None],
                 init_rms=[None, None], **kwargs):

        if transforms is None:
            transforms = OrderedDict()
        self.transforms = transforms
        self.order = order
        self.stats = OrderedDict()

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_solution(self, catalog1, catalog2,
                     src1, src2, dest1, dest2, order=None, match=False,
                     weights=None, method='statsmodels', separation=1 * apu.arcsec):
        """
        Get transformation from one set of pixel coordinates to another
        """
        from astropyp.catalog import Catalog

        if order is None:
            if self.order is None:
                raise Exception("You must give an order for the polynomial")
            order = self.order
        if match:
            if not isinstance(catalog1, Catalog):
                cat1 = Catalog(catalog1)
            else:
                cat1 = catalog1
            if not isinstance(catalog2, Catalog):
                cat2 = Catalog(catalog2)
            else:
                cat2 = Catalog(catalog2)
            # Match the catalog with the reference catalog
            cat_coords = coordinates.SkyCoord(cat1.ra, cat1.dec, unit='deg')
            refcat_coords = coordinates.SkyCoord(
                cat2.ra, cat2.dec, unit='deg')
            idx, d2, d3 = cat_coords.match_to_catalog_sky(refcat_coords)
            matched = d2 < separation
            result = OrderedDict()
            cat1 = catalog1[matched]
            cat2 = catalog2[idx][matched]
        else:
            if len(catalog1) != len(catalog2):
                raise Exception("Catalogs either must be the same length "
                                "or you must use matching")
            cat1 = catalog1
            cat2 = catalog2

        result1 = _map_coordinates(cat1[src1], cat1[src2], cat2[dest1],
                                   order, weights, method)
        result2 = _map_coordinates(cat1[src1], cat1[src2], cat2[dest2],
                                   order, weights, method)
        self.transforms[(src1, dest1)] = result1
        self.transforms[(src2, dest2)] = result2
        self.create_fit_stats(cat1, cat2, src1, src2, dest1, dest2)

        if match:
            result = (cat1, cat2, idx, matched)
        else:
            result = None
        return result

    def transform_coords(self, src1, src2, dest1, dest2, **kwargs):
        """
        Use the solution to transform the pixel coordinates

        Parameters
        ----------
        src1, src2: string
            Names of the source columns in the transformation
        dest1, dest2: string
            Names of the destination columns in the transformation
        kwargs: Dict of arrays
            Keyword Arguments for the transformation. The keys should
            always be the two source column names specified in ``src``
            and the values should be arrays for each column

        Result
        ------
        result: tuple of arrays
            Result of the astrometric solution transformation
        """
        result = _transform_coords(
            kwargs[src1], kwargs[src2], src1, src2,
            self.transforms[(src1, dest1)],
            self.transforms[(src2, dest2)])
        return result

    def create_fit_stats(self, catalog1, catalog2, src1, src2,
                         dest1, dest2):
        kwargs = {
            'src1': src1,
            'src2': src2,
            'dest1': dest1,
            'dest2': dest2,
            src1: catalog1[src1],
            src2: catalog1[src2]
        }
        coord1, coord2 = self.transform_coords(**kwargs)
        coord1_diff = coord1 - catalog2[dest1]
        coord2_diff = coord2 - catalog2[dest2]
        mean_diff = np.mean(coord1_diff)
        std_diff = np.std(coord1_diff)
        rms_diff = np.sqrt(mean_diff ** 2 + std_diff ** 2)
        self.stats[(src1, dest1)] = OrderedDict([
            ('mean', mean_diff),
            ('std', std_diff),
            ('rms', rms_diff)])
        mean_diff = np.mean(coord2_diff)
        std_diff = np.std(coord2_diff)
        rms_diff = np.sqrt(mean_diff ** 2 + std_diff ** 2)
        self.stats[(src2, dest2)] = OrderedDict([
            ('mean', mean_diff),
            ('std', std_diff),
            ('rms', rms_diff)])
        return self.stats


def _get_power(x_str, power):
    """
    Get the string for a power of a variable in `_build_formula`
    """
    if power > 1:
        return '{0}**{1}'.format(x_str, power)
    elif power == 1:
        return x_str
    return None


def _build_formula(x, y, order):
    """
    Build a polynomial formula of order ``order`` to
    use in statsmodels.OLS or WLS

    Parameters
    ----------
    x: array-like
        x Positions of sources
    y: array-like
        y Positions of sources
    order: int
        Order of the polynomial

    Result
    ------
    formula: str
        Formula to use in statsmodels OLS or WLS
    """
    formula = 'result ~ '
    first_term = True
    for l in range(1, order + 1):
        for m in range(l + 1):
            x_power = _get_power('x', l - m)
            y_power = _get_power('y', m)
            if not first_term:
                formula += '+'
            else:
                first_term = False
            if x_power is None:
                if m > 1:
                    formula += 'I({0})'.format(y_power)
                else:
                    formula += y_power
            elif y_power is None:
                if l - m > 1:
                    formula += 'I({0})'.format(x_power)
                else:
                    formula += x_power
            else:
                formula += 'I({0}*{1})'.format(x_power, y_power)
    logger.debug('formula: {0}'.format(formula))
    return formula


def _map_coordinates(coord1, coord2, ref_coord, order=3,
                     weights=None, method='statsmodels'):
    """
    Given a set of source coordinates and the reference coordinate,
    get the coefficients for the astrometric (or inverse) solution

    Parameters
    ----------
    coord1, coord2: array-like
        Arrays of source coordinates
    ref_coord: array-like
        Array of reference coordinates
    order: int
        Order of the solution polynomial
    param_name: string
        The name of the parameter we are solving for (for example:
        'x','y','ra','dec')
    prefix: string
        The prefix of the coefficient name
    fitpackage: string, optional
        Name of the fit package to use. Currently only the default
        'statsmodels' is available.
    """
    if method == 'statsmodels':
        try:
            import statsmodels.formula.api as smf
        except:
            raise Exception(
                "'statsmodels' package required for this type of "
                "astrometric solution")
        formula = _build_formula(coord1, coord2, order)
        tbl = table.Table([coord1, coord2, ref_coord], names=('x', 'y', 'result'))
        if weights is None:
            result = smf.OLS.from_formula(formula=formula,
                                          data=tbl).fit()
        else:
            result = smf.WLS.from_formula(formula=formula,
                                          data=tbl, weights=weights).fit()
        result = {'statsmodels': result}
    elif method == 'lst_sq':
        result = _lst_sq_fit(coord1, coord2, refcoord, order, weights)
        result = {'lst_sq': result}
    return result


def _transform_coords(coord1, coord2, coord1_str, coord2_str,
                      coord1_solution, coord2_solution):
    """
    Transform a set of coordinates to another frame

    Parameters
    ----------
    coord1: float or array of floats
        1st coordinate of each point
    coord2: float or array of floats
        2nd coordinate of each point
    coord1_str: string
        Name of the first coordinate (usually either 'x' or 'ra')
    coord2_str: string
        Name of the second coordinate (usually either 'y' or 'dec')
    coord1_solution: dict
        Result from _map_coordinates
    coord2_solution: dict
        Result from _map_coordinates

    Results
    -------
    coord1_result, coord2_result: float or array of floats
        Transformed coordinates
    """
    # from astropyp.catalog import Catalog
    ungroup = False
    try:
        len(coord1)
    except TypeError:
        coord1 = [coord1]
        coord2 = [coord2]
        ungroup = True
    catalog = table.Table([coord1, coord2],
                          names=('x', 'y'))
    # names=(coord1_str,coord2_str))
    # catalog = Catalog(catalog)

    if 'statsmodels' in coord1_solution:
        solution1 = coord1_solution['statsmodels'].predict
        solution2 = coord2_solution['statsmodels'].predict
    else:
        raise ValueError("Unrecognized astrometric solution")

    coord1_result = solution1(catalog)
    coord2_result = solution2(catalog)
    if ungroup:
        coord1_result = coord1_result[0]
        coord2_result = coord2_result[0]
    return coord1_result, coord2_result


#########################################################
# All Code below will be depreciated in futue releases
#########################################################

class AstrometricSolution:
    def __init__(self, x2ra=None, y2dec=None, ra2x=None, dec2y=None,
                 crpix=[0, 0], crval=[0, 0], order=3, **kwargs):
        self.x2ra = x2ra
        self.y2dec = y2dec
        self.ra2x = ra2x
        self.dec2y = dec2y
        self.crpix = crpix
        self.crval = crval
        self.order = order

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_solution(self, cat, refcat, order=None,
                     match=True, separation=1 * apu.arcsec, weights=None,
                     ref_weights=None):

        if order is None:
            if self.order is None:
                raise Exception("You must give an order for the polynomial")
            order = self.order

        if match:
            # Match the catalog with the reference catalog
            cat_coords = coordinates.SkyCoord(cat.ra, cat.dec, unit='deg')
            refcat_coords = coordinates.SkyCoord(
                refcat.ra, refcat.dec, unit='deg')
            idx, d2, d3 = cat_coords.match_to_catalog_sky(refcat_coords)
            matched = d2 < separation
            result = OrderedDict()

        x = cat.x[matched] - self.crpix[0]
        y = cat.y[matched] - self.crpix[1]
        ra = np.array(refcat.ra[idx][matched]) - self.crval[0]
        dec = np.array(refcat.dec[idx][matched]) - self.crval[1]
        if weights is not None:
            weights = weights[matched]
        if ref_weights is not None:
            ref_weights = ref_weights[idx][matched]
            if weights is None:
                weights = ref_weights
            else:
                weights = 1. / np.sqrt(1. / weights ** 2 + 1. / ref_weights ** 2)

        # Get RA fit
        result = _get_solution(x, y, ra, order, 'ra', 'A', weights=weights)
        self.x2ra = result['params']
        del result['params']
        self.x2ra_info = result
        # Get DEC fit
        result = _get_solution(y, x, dec, order, 'dec', 'B', weights=weights)
        self.y2dec = result['params']
        del result['params']
        self.y2dec_info = result

        # Get x fit
        ra, dec = self.pix2world(x + self.crpix[0], y + self.crpix[1])
        ra = ra - self.crval[0]
        dec = dec - self.crval[1]
        result = _get_solution(ra, dec, x, order, 'x', 'Ap', weights=weights)
        self.ra2x = result['params']
        del result['params']
        self.ra2x_info = result
        # Get y fit
        result = _get_solution(dec, ra, y, order, 'y', 'Bp', weights=weights)
        self.dec2y = result['params']
        del result['params']
        self.dec2y_info = result

        return idx, matched

    def pix2world(self, x, y):
        """
        Convert pixel coordinates to wcs using
        an astropyp astrometric solution.

        Parameters
        ----------
        x,y: array-like
            Pixel coordinates of the sources

        Result
        ------
        ra, dec
        """
        ra_tbl = get_transform_tbl(x - self.crpix[0], y - self.crpix[1],
                                   self.order, 'A')
        dec_tbl = get_transform_tbl(y - self.crpix[1], x - self.crpix[0],
                                    self.order, 'B')
        ra_tbl['Intercept'] = 1
        dec_tbl['Intercept'] = 1
        ra = 0
        dec = 0
        for column, value in self.x2ra.items():
            ra += ra_tbl[column] * value
        for column, value in self.y2dec.items():
            dec += dec_tbl[column] * value
        return ra + self.crval[0], dec + self.crval[1]

    def world2pix(self, ra, dec):
        """
        Convert pixel coordinates to wcs using
        an astropyp astrometric solution.

        Parameters
        ----------
        x,y: array-like
            Pixel coordinates of the sources

        Result
        ------
        ra, dec
        """
        x_tbl = get_transform_tbl(ra - self.crval[0], dec - self.crval[1],
                                  self.order, 'Ap')
        y_tbl = get_transform_tbl(dec - self.crval[1], ra - self.crval[0],
                                  self.order, 'Bp')
        x_tbl['Intercept'] = 1
        y_tbl['Intercept'] = 1
        x = 0
        y = 0
        for column, value in self.ra2x.items():
            x += x_tbl[column] * value
        for column, value in self.dec2y.items():
            y += y_tbl[column] * value
        return x + self.crpix[0], y + self.crpix[1]


class ImageSolution:
    """
    ImageSolution can be used to calculate transformations
    from one image to another (ie x coordinates to y coordinates)
    """

    def __init__(self, x_tx=None, y_tx=None, order=None,
                 init_mean=[None, None], init_stddev=[None, None],
                 init_rms=[None, None], **kwargs):
        self.x_tx = x_tx
        self.y_tx = y_tx
        self.order = order
        self.mean = init_mean
        self.stddev = init_stddev
        self.rms = init_rms

        for k, v in kwargs.items():
            setattr(self, k, v)

    @property
    def stats(self):
        stats = self.get_fit_stats()
        result = {'mean': stats[0], 'stddev': stats[1], 'rms': stats[2]}
        return result

    def get_solution(self, catalog1, catalog2, order=None, match=False,
                     weights=None, fit_package='statsmodels'):
        """
        Get transformation from one set of pixel coordinates to another
        """
        from astropyp.catalog import Catalog

        if order is None:
            if self.order is None:
                raise Exception("You must give an order for the polynomial")
            order = self.order
        if match:
            # Match the catalog with the reference catalog
            cat_coords = coordinates.SkyCoord(cat.ra, cat.dec, unit='deg')
            refcat_coords = coordinates.SkyCoord(
                refcat.ra, refcat.dec, unit='deg')
            idx, d2, d3 = cat_coords.match_to_catalog_sky(refcat_coords)
            matched = d2 < separation
            result = OrderedDict()
            cat1 = catalog1[matched]
            cat2 = catalog2[idx][matched]
        else:
            if len(catalog1) != len(catalog2):
                raise Exception("Catalogs either must be the same length "
                                "or you must use matching")
            cat1 = catalog1
            cat2 = catalog2
        # Make sure that you are using Catalogs
        if not isinstance(cat1, Catalog):
            cat1 = Catalog(cat1)
        if not isinstance(cat2, Catalog):
            cat2 = Catalog(cat2)

        x_result = _get_solution(cat1.x, cat1.y, cat2.x,
                                 order, 'ref_x', prefix='A', weights=weights,
                                 fit_package=fit_package)
        y_result = _get_solution(cat1.y, cat1.x, cat2.y,
                                 order, 'ref_y', prefix='B', weights=weights,
                                 fit_package=fit_package)
        self.x_tx = x_result['params']
        self.y_tx = y_result['params']
        del x_result['params']
        del y_result['params']
        self.x_info = x_result
        self.y_info = y_result
        if match:
            result = (cat1, cat2, idx, matched)
        else:
            result = None
        self.create_fit_stats(cat1, cat2)
        return result

    def transform_coords(self, catalog=None, x=None, y=None):
        """
        Use the solution to transform the pixel coordinates
        """
        from astropyp.catalog import Catalog
        ungroup = False
        if catalog is None:
            try:
                len(x)
            except TypeError:
                x = [x]
                y = [y]
                ungroup = True
            catalog = table.Table([x, y], names=('x', 'y'))
        if not isinstance(catalog, Catalog):
            cat = Catalog(catalog)
        else:
            cat = catalog
        x_tbl = get_transform_tbl(cat.x, cat.y, self.order, 'A')
        x_tbl['Intercept'] = 1
        x = 0
        for column, value in self.x_tx.items():
            x += x_tbl[column] * value
        y_tbl = get_transform_tbl(cat.y, cat.x, self.order, 'B')
        y_tbl['Intercept'] = 1
        y = 0
        for column, value in self.y_tx.items():
            y += y_tbl[column] * value
        if ungroup:
            x = x[0]
            y = y[0]
        return x, y

    def create_fit_stats(self, catalog1, catalog2):
        from astropyp.catalog import Catalog
        if not isinstance(catalog2, Catalog):
            cat2 = Catalog(catalog2)
        else:
            cat2 = catalog2
        x, y = self.transform_coords(catalog1)
        x_diff = x - cat2.x
        y_diff = y - cat2.y
        self.mean = [np.mean(x_diff), np.mean(y_diff)]
        self.stddev = [np.std(x_diff), np.std(y_diff)]
        self.rms = [np.sqrt(np.sum(x_diff ** 2) / len(x_diff)),
                    np.sqrt(np.sum(y_diff ** 2) / len(y_diff))]
        return self.get_fit_stats()

    def get_fit_stats(self):
        return self.mean, self.stddev, self.rms


def get_transform_tbl(x, y, order, prefix='A'):
    """
    Given a set of x and y positions, create a table with the different
    polynomial terms [x,y,x**2,x*y,y**2,...] as columns

    Parameters
    ----------
    x: array-like
        x Positions of sources
    y: array-like
        y Positions of sources
    order: int
        Order of the polynomial
    prefix: string, optional
        Name of the coefficient. *Default='A'*

    Result
    ------
    tbl: `~astropy.table.Table`
        Table with coefficient names as column names and
        polynomial terms as columns
    """
    vec = []
    columns = []
    for l in range(1, order + 1):
        for m in range(l + 1):
            vec.append(x ** (l - m) * y ** m)
            columns.append('{0}_{1}_{2}'.format(prefix, l - m, m))
    tbl = table.Table(vec, names=columns)
    tbl.meta['order'] = order
    return tbl


def build_formula(param_name, transform_tbl):
    """
    Build a formula for statsmodels to use to calculate the
    coefficients.

    Parameters
    ----------
    param_name: string
        Name of the reference column in ``transform_tbl``
        to solve for.
    transform_tbl: `~astropy.table.Table`
        Table of polynomial terms with coefficients as
        column names created by
        `~astropyp.astrometry.get_transform_tbl`.

    Result
    ------
    formula: string
        formula to use in statsmodels
    """
    formula = "{0} ~ ".format(param_name)
    formula += '+'.join(transform_tbl.columns.keys())
    logger.debug('formula: "{0}"'.format(formula))
    return formula


def _statsmodels_to_result(params):
    result = OrderedDict(
        [(column, value)
         for column, value in params.iteritems()])
    return result


def _lst_sq_fit(coord1, coord2, ref_coord, order, weights=None, prefix='A'):
    from astropyp.astrometry import get_transform_tbl
    A = get_transform_tbl(coord1, coord2, order, prefix)
    columns = ['Intercept'] + A.columns.keys()
    tbl_dtype = A[A.columns.keys()[0]].dtype
    A['Intercept'] = np.ones((len(A)), dtype=tbl_dtype)
    A = A[columns]
    A = A.as_array().view(tbl_dtype).reshape((len(A), len(A.columns)))
    y = ref_coord
    if weights is not None:
        w = np.vstack([weights for n in range(A.shape[1])]).T
        A = A * w
        y = y * weights
    result = np.linalg.lstsq(A, y)[0]
    solution = OrderedDict([
                               (columns[n], result[n]) for n in range(len(columns))])
    return solution


def _get_solution(coord1, coord2, ref_coord, order, param_name,
                  prefix='A', fit_package='statsmodels', weights=None):
    """
    Given a set of source coordinates and the reference coordinate,
    get the coefficients for the astrometric (or inverse) solution

    Parameters
    ----------
    coord1, coord2: array-like
        Arrays of source coordinates
    ref_coord: array-like
        Array of reference coordinates
    order: int
        Order of the solution polynomial
    param_name: string
        The name of the parameter we are solving for (for example:
        'x','y','ra','dec')
    prefix: string
        The prefix of the coefficient name
    fitpackage: string, optional
        Name of the fit package to use. Currently only the default
        'statsmodels' is available.
    """
    if fit_package == 'statsmodels':
        try:
            import statsmodels.formula.api as smf
        except:
            raise Exception(
                "'statsmodels' package required for this type of "
                "astrometric solution")
        tbl = get_transform_tbl(coord1, coord2, order, prefix)
        formula = build_formula(param_name, tbl)
        tbl[param_name] = ref_coord
        if weights is None:
            result = smf.OLS.from_formula(formula=formula,
                                          data=tbl).fit()
        else:
            result = smf.WLS.from_formula(formula=formula,
                                          data=tbl, weights=weights).fit()
        result = {
            'params': _statsmodels_to_result(result.params),
            'statsmodels': result
        }
    elif fit_package == 'lst_sq':
        result = _lst_sq_fit(coord1, coord2, ref_coord, order, weights, prefix)
        result = {'params': result}
    else:
        raise ValueError("fit_package must be either statsmodels' or"
                         " 'lst_sq'")
    return result


def clear_header_wcs(header):
    wcs_keys = ['CTYPE', 'CRVAL', 'CRPIX', 'CD1_', 'CD2_',
                'PV', 'A_', 'AP_', 'B_', 'BP_']
    for k in header.keys():
        if any([k.startswith[key] for key in wcs_keys]):
            del header[k]
    return header


def update_sip_header(header, order, result):
    """
    Update a header with new WCS information
    """
    crval = [header['CRVAL1'], header['CRVAL2']]
    crpix = [header['CRPIX1'], header['CRPIX2']]
    header = clear_header_wcs(header)
    header['CTYPE1'] = 'RA---TAN-SIP'
    header['CTYPE2'] = 'DEC--TAN-SIP'
    header['CRPIX1'], header['CRPIX2'] = crpix
    header['CRVAL1'], header['CRVAL2'] = crval
    if 'ra' in result:
        header['A_ORDER'] = order
        header['B_ORDER'] = order
        for k, v in result['ra'].params.iteritems():
            if k not in ['Intercept']:
                header[k] = v
        for k, v in result['dec'].params.iteritems():
            if k not in ['Intercept']:
                header[k] = v
        header['CD1_1'] = result['ra'].params.A_1_0
        header['CD1_2'] = result['ra'].params.A_0_1

        # Temporary, delete this later
        header['CRPIX1'] = 0
        header['CRPIX2'] = 0
        header['CRVAL1'] = 0  # result['ra'].Intercept
        header['CRVAL2'] = 0  # result['dec'].Intercept
    return header


def get_astrometric_solution(cat, refcat, order, direction='both'):
    """
    Get the astrometric solution from an image catalog to a reference
    catalog.

    Parameters
    ----------
    cat: `~astropyp.catalog.Catalog`
        Point source catalog with x,y,ra,dec defined.
    refcat: `~astropyp.catalog.Catalog`
        Reference catalog with ra,dec defined for 'pix2wcs'
        and ra,dec,x,y defined to 'world2pix'
    order: int
        Order of the polynomial solution
    direction: string, optional
        Direction to calculate solution. This can be
        'pix2wcs', 'wcs2pix', 'both'. *Default='both'*
    header: astropy FITS Header, optional
        If a header is given, WCS header parameters are calculated and
        the header is updated. If no header is given then this
        step is skipped. *Default=None*
    """
    calc_world = (direction == 'pix2wcs' or direction == 'both')
    calc_pix = (direction == 'wcs2pix' or direction == 'both')

    cat_coords = coordinates.SkyCoord(cat.ra, cat.dec, unit='deg')
    refcat_coords = coordinates.SkyCoord(refcat.ra, refcat.dec, unit='deg')
    idx, d2, d3 = cat_coords.match_to_catalog_sky(refcat_coords)
    matched = d2 < 1 * apu.arcsec
    result = OrderedDict()
    result['idx'] = idx
    result['matched'] = matched
    if calc_world:
        x = cat.x[matched]
        y = cat.y[matched]
        ra = refcat.ra[idx][matched]
        dec = refcat.dec[idx][matched]
        result['ra'] = _get_solution(x, y, ra, order, 'ra', 'A')
        result['dec'] = _get_solution(y, x, dec, order, 'dec', 'B')
    if calc_pix:
        ra = cat.ra[matched]
        dec = cat.dec[matched]
        x = refcat.x[idx][matched]
        y = refcat.y[idx][matched]
        result['x'] = _get_solution(ra, dec, x, order, 'x', 'Ap')
        result['y'] = _get_solution(dec, ra, y, order, 'y', 'Bp')
    return result

# -----------------------------------------------------------------
