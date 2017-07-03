##############################
# Exact deprojections of Sersic surface brightness profiles following 
# Baes and Gentile, Astronomy and Astrophysics, Volume 525:A136 (2011)
# http://arxiv.org/abs/1009.4713
# 
# Very likely the only function any user will care about is:
#
# def luminosity(pp, qq, reff=1.0, lum=1.0) where pp and qq are the
# numerator and denominator of the Sersic index (both integers) so
# that nn=pp/qq, reff is the projected half-light radius, and
# luminosity is the total luminosity.  This returns a function that
# takes a radius and returns a luminosity density.
# 
# >>> lum = luminosity(5,3)
# >>> lum(1.1)
# >>> lum([1.1, 2.2, 3.3])
# 
# All of the other functions in the module are just to help make sure
# that the Baes + Gentile expressions are implemented correctly.
# 
# If you find that the calculations are taking a long time, you can
# try reducing the precision of the calculation of the Meijer G
# function.  The mpmath library is designed to be able to work to
# arbitrary precision, and is by default set to ~double precision.
# For information on changing the precision, see:
# http://mpmath.googlecode.com/svn/trunk/doc/build/basics.html#setting-the-precision
# Keep in mind that small values of the numerator and denominator are
# faster than large ones (e.g. luminosity(2,1) is fast,
# luminosity(200,99) is slow.
#
# The Baes + Gentile expressions use the Meijer G function to compute
# exact deprojectsions for rational Sersic indicies.  They provide
# expresssions for irrational Sersic indicies using the Fox H
# function, but I couldn't find an implementation of the Fox H
# function so the present code is restricted to rational Sersic
# indicies.  I haven't found this to be a problem since I can find a
# rational number arbitrarly close to any desired value (although
# computation of the Meijer G function starts to take a long time as
# the numerator and denominator become large.  It would be a problem
# if one were, for instance, trying to do a fit to some 3D luminosity
# density using the expressions provided.  However, generally one fits
# to the surface brightness distribution (where the restriction to
# integers doens't apply) and then evaluates deprojected luminosity
# for the specific desired Sersic index.
#
# Naming conventions
#
# I wanted to be rather pedantic about making sure that I had coded
# this up correctly, and that the results had the expected properties
# (that the projection integral of the 3d luminosity density is indeed
# the 2d surface brightness that it's supposed to be).  Therefore I
# end up using three different parameterizations of the Sersic profile
# and cross checking them.  I also wanted the notation of the code to
# directly reflect the notation used in the papers.  This led to some
# abbreviations referring to symbols in defined in the papers.
#
# The bg, bm, and ml prefixes refer to the three parameterizations:
# bg: Baes and Gentile (http://arxiv.org/abs/1009.4713)
# bm: Binney + Merrifield (Galactic Astronomy, Princeton Univ. Press)
# ml: Mamon + Lokas (http://arxiv.org/abs/astro-ph/0405466)
#
# Short cryptic letter combinations that refer to the papers are:
# ie: I_e, surface brightness at the effective radius, defined by B+M
# bn: b_n, appears in the exponent of the Sersic profile in the B+M
#     parameterization and must be adusted as a function of the Sersic
#     index n to ensure that the half-light radius is what you want it
#     to be
# bb: b appearing in the exponent of the Sersic profile in the B+G
#     parameterization.  It's doubled just to avoid a single-character
#     variable name.

import operator

import scipy.integrate, scipy.optimize
#import mpmath ## temporarily disable this import statement for 'pts depends'
import numpy as np

# Maximum multiple of effective radius to which to carry integrals.
# In principle should be able to set this to inf, but that causes the
# integrals to take a long time.
#reff_max = np.inf
reff_max = 1000.0

##############################
## Baes + Gentile formulas work for rational numbers, so it's
## convenient to have a few functions to reduce fractions, etc.
##############################

def euclid_gcf(a,b):
    """Find greatest common factor by Euclid's algorithm"""
    def gcf_internal(a,b):
        if a < b: return gcf_internal(b,a)
        if b==0: return a
        return gcf_internal(b, a%b)

    if a<=0 or b<=0: raise RuntimeError
    if int(a)!=a or int(b)!=b: raise RuntimeError
    return gcf_internal(a,b)

def all_rationals(qmax, rmin=None, rmax=None, sort=False):
    """Return tuples representing all rational numbers with denom < qmax"""
    # Seems to scale as qmax^3, takes ~40 sec for qmax=2000, generates
    # about a million numbers.
    if rmin is None and rmax is None: rmin, rmax = 0, 1
    elif rmax is None: rmin, rmax = 0, rmin
    elif rmin is None: raise SyntaxError
    
    result = [ [(pp,qq) for pp in range(1,int(rmax*qq+1))
                if euclid_gcf(pp,qq) == 1 and pp/(1.0*qq) >= rmin and pp/(1.0*qq) <= rmax]
               for qq in range(1,qmax+1)]
    # flatten result into a single list
    result = reduce(operator.add, result)
    # Sort if requested
    if sort: result = sorted(result, key=lambda (pp,qq): pp/(1.0*qq))
    return result

##############################
## General properties: projection integrals, etc
##############################

def lum_int_2d(ir, rpi, rpo):
    # Tested 2011-08-30
    "Total projected lum between rpi and rpo"    
    return scipy.integrate.quad(lambda x: 2*np.pi*x*ir(x), rpi, rpo)[0]

def lum_int_3d(lr, ri, ro):
    # Tested 2011-08-30
    "Total lum between ri and ro"
    return scipy.integrate.quad(lambda x: 4*np.pi*x**2*lr(x), ri, ro)[0]

def lum_proj(lr, rp, ro=reff_max):
    # Tested 2011-08-30
    "Project a 3d lum dens to a 2d surface brightness"
    return 2*scipy.integrate.quad(lambda xx: lr(np.sqrt(rp**2 + xx**2)), 0, ro)[0]

def half_light_2d(ir,ro=reff_max):
    # Tested 2011-08-30
    "Find half light radius of surf brightness profile"
    ltot = lum_int_2d(ir, 0.0, ro)
    return scipy.optimize.bisect(lambda xx: lum_int_2d(ir, 0.0, xx)-0.5*ltot,
                                 0.01, 100.0)

def half_light_3d(lr, ro=reff_max):
    # Tested 2011-08-30
    "Find half light radius of 3d luminosity density profile"
    ltot = lum_int_3d(lr, 0.0, ro)
    return scipy.optimize.bisect(lambda xx: lum_int_3d(lr, 0.0, xx)-0.5*ltot,
                                 0.01, 100.0)
    
##############################
## Binney + Merrifield parameterization.
##############################

def bm(rp, nn, ie, bn):
    # Tested 2011-08-30 through bm_func
    "Projected sersic profile using parameterization from B+M"
    return ie*10**(-bn*(rp**(1.0/nn)-1))

def bm_proj_lum_int(rp, nn, ie, bn):
    # Tested 2011-08-30 through bm_func
    "Total lum of projected sersic profile interior to rp"    
    return lum_int_2d(lambda rr: bm(rr, nn, ie, bn), 0.0, rp)

def bm_bn_estimate(nn):
    "Guess for bn constant defined by B+M"
    if np.iterable(nn): return np.array([bm_bn_estimate(n) for n in nn])
    est = 0.87*nn-0.15
    if nn < 0.2: est = 0.01*(nn/0.1)**4
    return est

def bm_bn(nn, ro=reff_max):    
    # Tested 2011-08-30 through bm_func
    "Find value of bn to make reff actually reff"
    assert nn >= 0.014 and nn <= 29.99  # Looks like numerical noise outside this range

    def ff(xx):
        return bm_proj_lum_int(1.0, nn, 1.0, xx) - 0.5*bm_proj_lum_int(ro, nn, 1.0, xx)
    # Rather convoluted estimates of b and boundaries of search.
    bbi = bm_bn_estimate(nn)
    if nn < 0.1: bl, bh = 1e-30, 10*bbi
    elif nn < 0.45: bl, bh = 1e-10, 10*bbi
    elif nn < 5.0: bl, bh = bbi/10, 10*bbi
    elif nn < 30.0: bl, bh = bbi/2.0, 2*bbi
    else: bl, bh = bbi/1.2, 1.2*bbi

    return scipy.optimize.bisect(ff, bl, bbi*10.0)

def bm_ie(nn, bn, ro=reff_max):
    # Tested 2011-08-30 through bm_func
    "Find value of ie to make luminosity unity"
    return 1.0/bm_proj_lum_int(ro, nn, 1.0, bn)

def bm_func(nn, ro=reff_max):
    # Tested 2011-08-30
    "Projected sersic profile"
    bn = bm_bn(nn, ro=ro)
    ie = bm_ie(nn, bn, ro=ro)    
    return lambda xx: bm(xx, nn, ie, bn)

##############################
## Mamon + Lokas (astro-ph/0405466) parameterization.
##############################

def ml(rp, nn, i0, rs):
    # Tested 2011-08-30
    "Mamon + Lokas parameterization of projected sersic profile."
    return i0*np.exp(-(rp/rs)**(1.0/nn))

def ml_constants(nn, ro=reff_max):
    "Compute constants for M+L parameterization of Sersic profile"
    # Tested through ml
    bn = bm_bn(nn, ro=ro)
    ie = bm_ie(nn, bn, ro=ro)    
    reff = 1.0
    i0 = ie*10**bn
    rs = reff/(bn*np.log(10))**nn
    return i0, rs

def ml_lum_func(nn, ro=reff_max):
    # Tested 2011-08-30
    """Approx to lum density from progniel + Simien 97 as modified by
    Lima Neto et al 99, nicely written by Mamon + Lokas 05"""
    i0, rs = ml_constants(nn, ro=ro)
    gamma = scipy.special.gamma    
    pp = 1.0-0.6097/nn + 0.05463/nn**2
    l1 = (gamma(2.0*nn)/gamma((3.0-pp)*nn)) * 0.5*i0/rs

    def lum(rr):
        xx = rr/rs
        return l1*xx**(-pp)*np.exp(-xx**(1.0/nn))
        
    return lum

##############################
## Baes + Gentile (arxiv:1009.4713) parameterization.
##############################

def bg(rp, nn, i0, bb):
    # Tested 2011-08-30
    """Sersic profile as parameterized by Baes + Gentile"""
    return i0*np.exp(-bb*rp**(1.0/nn))

def bg_constants_from_bm(nn, ro=reff_max):
    # Tested 2011-08-30 through bg
    """Find Baes + Gentile constants by first finding B+M constants numerically"""
    bn = bm_bn(nn, ro=ro)
    ie = bm_ie(nn, bn, ro=ro)    
    reff = 1.0
    i0 = ie*10**bn
    bb = np.log(10)*bn
    return i0,bb

def bg_bb_estimate(pp,qq=None):
    """Estimate b constant defined by B+G"""
    # Allow input as both float or as num/denom of rational number.
    if qq is None: nn = pp
    else: nn = pp/(1.0*qq)

    # Seems to capture asymptotic behavior very well with only a little wiggle around n~0.8
    est = np.exp(-1/(2.7*nn))
    est = np.where(nn < 0.7, est, est + 2*nn-0.33-1)
    return est
    
def bg_constants(pp, qq):
    # Tested 2011-08-31
    """Find Baes + Gentile constants directly from their expressions"""
    nn = (1.0*pp)/qq

    assert 1.0/30 <= nn and nn <= 71.0
    # Initial guesses
    bbi = bg_bb_estimate(pp,qq)
    i0i = 1.0

    # find bb
    reff=1.0
    def ff(xx): return (lum_int_2d(lambda rp: bg(rp, nn, i0i, xx), 0.0, reff) - 
                        0.5*bg_lum_tot(pp,qq,i0i,xx))
    bb = scipy.optimize.newton(ff, bbi)

    # find i0
    i0 = i0i/bg_lum_tot(pp,qq,i0i,bb)
    return i0,bb
                                                  
def bg_3d_lum_int_func(pp,qq):
    """Return a function that gives deprojected sersic profile"""
    i0,bb = bg_constants(pp, qq)
    def ff(xx):
        if np.iterable(xx):
            return [ff(x) for x in xx]
        return float(bg_3d_lum_int(xx,pp,qq,i0,bb))
    return ff

def bg_3d_lum_int(rr,pp,qq,i0,bb):
    # Tested 2011-08-31
    """Fully analytic calculation of lum internal to rr"""
    if not (pp==int(pp) and qq==int(qq)): raise RuntimeError
    pp, qq = int(pp), int(qq)

    reff = 1.0
    ss = rr/reff
    avect = [[1-1.0/qq], [xx/(1.0*qq) for xx in range(1,qq)]]
    bvect = [[xx/(2.0*pp) for xx in range(1,2*pp)] +
             [xx/(2.0*qq) for xx in range(1,2*qq,2)], [-1.0/qq]]    
    factor = 2*i0*reff**2*np.sqrt(pp)/((2*np.pi)**(pp-1)*np.sqrt(qq))    
    zz = (bb/(2*pp))**(2*pp) * ss**(2*qq)
    return factor*ss**2*mpmath.meijerg(avect, bvect, zz)

def bg_lum_tot(pp,qq,i0,bb):
    # Tested 2011-08-30
    """Fully analytic calculation of total lum of sersic profile"""
    reff = 1.0
    mm = (1.0*pp)/qq
    return np.pi*i0*reff**2*scipy.special.gamma(2*mm+1) / (bb**(2*mm))

def bg_lum_hi(pp, qq):
    # Tested 2011-08-31
    """Exact deprojection of Sersic profile using Meijer G function as
    described by Baes + Gentile 1009.4713.  Use formula valid for half
    integers"""

    if not (pp==int(pp) and qq==int(qq)): raise RuntimeError
    if not (qq == 1 or qq == 2): raise RuntimeError

    pp, qq = int(pp), int(qq)    
    mm = (1.0*pp)/qq
    i0, bb = bg_constants(pp, qq)

    # a and b vectors are specified: [[num1, num2, num3],[denom1,denom2,denom3]]    
    avect = [[], []]
    nums = range(1,2*pp/qq)
    bvect = [[xx/(2.0*mm) for xx in nums] + [0.5], []]
    reff = 1.0

    factor = 2*i0*np.sqrt(mm)/(reff*(2*np.pi)**mm)

    def lum(rr):
        if np.iterable(rr): return np.array([lum(r) for r in rr])
        ss = rr/reff
        zz = (bb/(2*mm))**(2*mm) * ss**2
        return (factor/ss)*mpmath.meijerg(avect, bvect, zz)

    return lum

##################################################
## The big money function is defined right here!  The name conforms to
## naming conventions used throughout this file, but is nearly useless
## to users.  So make it available as luminosity() to users via a line
## in __init__.py since it's likely the only function any user will
## care about.
def bg_lum(pp, qq, reff=1.0, lum=1.0):
    # Tested 2011-08-31
    """Exact deprojection of Sersic profile using Meijer G function as
    described by Baes + Gentile arxiv:1009.4713.

    pp and qq are the numerator and denominator of the Sersic index
      (both integers) so that n=pp/qq, 
    reff is the projected half-light radius
    lum is the total luminosity.  

    This returns a function that takes a radius and returns a
    luminosity density.

    >>> lum = luminosity(5,3)
    >>> lum(1.1)
    >>> lum([1.1, 2.2, 3.3])

    """
    if not (pp==int(pp) and qq==int(qq)): raise RuntimeError
    pp, qq = int(pp), int(qq)
    i0, bb = bg_constants(pp, qq)

    # Solution gets slow for larger p,q, so make sure that fraction is reduced
    the_gcf = euclid_gcf(pp,qq)    
    if the_gcf != 1: return bg_lum(pp/the_gcf, qq/the_gcf)        

    # a and b vectors are specified: [[num1, num2, num3],[denom1,denom2,denom3]]    
    avect = [[], [xx/(1.0*qq) for xx in range(1,qq)]]
    bvect = [[xx/(2.0*pp) for xx in range(1,2*pp)] +
             [xx/(2.0*qq) for xx in range(1,2*qq,2)], []]
    factor = 2*i0*np.sqrt(pp*qq)/(reff*(2*np.pi)**pp)

    def luminosity(rr):
        if np.iterable(rr): return np.array([luminosity(r) for r in rr])
        ss = rr/reff
        zz = (bb/(2*pp))**(2*pp) * ss**(2*qq)
        return lum*((factor/ss)*mpmath.meijerg(avect, bvect, zz))

    return luminosity
