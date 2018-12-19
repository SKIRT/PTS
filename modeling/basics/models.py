#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.models Contains the SersicModel, ExponentialDiskModel and DeprojectionModel classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from abc import ABCMeta, abstractmethod, abstractproperty
from scipy.special import gammaincinv
from scipy import ndimage
import itertools
import warnings

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite
from ...core.units.parsing import parse_unit as u
from ...magic.core.frame import Frame
from ...core.tools import filesystem as fs
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from ..basics.instruments import FrameInstrument, SEDInstrument, SimpleInstrument, FullInstrument
from ...magic.basics.pixelscale import Pixelscale
from ...core.tools import numbers
from ...magic.basics.vector import PixelShape
from ...magic.basics.coordinate import PixelCoordinate
from ...core.basics.range import QuantityRange
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

default_truncation = 4.

# -----------------------------------------------------------------

sersic = "sersic"
exponential = "exponential"
deprojection = "deprojection"

# -----------------------------------------------------------------

models_2D = [sersic, exponential]
models_3D = [sersic, exponential, deprojection]

# -----------------------------------------------------------------

class Model(SimplePropertyComposite):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @abstractproperty
    def xmin(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def xmax(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def xrange(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def x_range(self):

        """
        This function ...
        :return:
        """

        return self.xrange

    # -----------------------------------------------------------------

    @abstractproperty
    def ymin(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def ymax(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def yrange(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def y_range(self):

        """
        This function ...
        :return:
        """

        return self.yrange

# -----------------------------------------------------------------

class Model2D(Model):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

# -----------------------------------------------------------------

class Model3D(Model):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @abstractproperty
    def zmin(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def zmax(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def zrange(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def z_range(self):

        """
        This function ...
        :return:
        """

        return self.zrange

# -----------------------------------------------------------------
#
# 3D MODELS
#
# -----------------------------------------------------------------

def load_3d_model(path):

    """
    This function ...
    :param path:
    :return:
    """

    first_line = fs.get_first_line(path)

    # Create the appropriate model
    if "SersicModel3D" in first_line: return SersicModel3D.from_file(path)
    elif "ExponentialDiskModel3D" in first_line: return ExponentialDiskModel3D.from_file(path)
    elif "DeprojectionModel3D" in first_line: return DeprojectionModel3D.from_file(path)
    else: raise ValueError("Unrecognized model file")

# -----------------------------------------------------------------

def create_3d_model(model_type, **properties):

    """
    This function ...
    :param model_type:
    :param properties:
    :return:
    """

    # Create the appropriate model
    if model_type == sersic: return SersicModel3D(**properties)
    elif model_type == exponential: return ExponentialDiskModel3D(**properties)
    elif model_type == deprojection: return DeprojectionModel3D(**properties)
    else: raise ValueError("Unrecognized model type: " + model_type)

# -----------------------------------------------------------------

# FROM SKIRT:
#SersicFunction::SersicFunction(double n)

    # if (n<0.5 || n>10.0)
    #     throw FATALERROR("The Sersic parameter should be between 0.5 and 10 (n = " + QString::number(n) + ")");
    # double b = 2.0*n - 1.0/3.0 + 4.0/405.0/n + 46.0/25515.0/(n*n) + 131.0/1148175.0/(n*n*n);
    # double I0 = pow(b,2.0*n) / (M_PI*SpecialFunctions::gamma(2.0*n+1));
    # int Ns = 101;
    # _sv.resize(Ns);
    # _Sv.resize(Ns);
    # _Mv.resize(Ns);
    # double logsmin = -6.0;
    # double logsmax = 4.0;
    # double dlogs = (logsmax-logsmin)/(Ns-1.0);
    # for (int i=0; i<Ns; i++)
    # {
    #     double logs = logsmin + i*dlogs;
    #     double s = pow(10.0,logs);
    #     _sv[i] = s;
    #     double alpha = b*pow(s,1.0/n);
    #     double sum = 0.0;
    #     int Nu = 10000;
    #     double tmax = 100.0;
    #     double umax = sqrt((tmax+1.0)*(tmax-1.0));
    #     double du = umax/Nu;
    #     for (int j=0; j<=Nu; j++)
    #     {
    #         double weight = 1.0;
    #         if (j==0 || j==Nu) weight=0.5;
    #         double u = j*du;
    #         double u2 = u*u;
    #         double w;
    #         if (u>1e-3)
    #             w = (pow(1.0+u2,2.0*n)-1.0)/u2;
    #         else
    #             w = 2.0*n + n*(2.0*n-1.0)*u2 + 2.0/3.0*n*(2.0*n-1.0)*(n-1.0)*u2*u2;
    #         double integrandum = 2.0*exp(-alpha*(1.0+u2))/sqrt(w);
    #         sum += weight*integrandum;
    #     }
    #     _Sv[i] = I0 * pow(b,n) * pow(alpha,1.0-n) / M_PI * du * sum;
    # }
    #
    # // calculate the cumulative mass
    #
    # for (int i=1; i<Ns; i++)
    # {
    #     double sum = 0.0;
    #     for (int j=0; j<=32; j++)
    #     {
    #         double weight=1.0;
    #         if (j==0 || j==32) weight=0.5;
    #         double ds = (_sv[i]-_sv[i-1])/32.0;
    #         double s = _sv[i-1] + j*ds;
    #         double S = operator()(s);
    #         sum += weight*S*s*s*ds;
    #     }
    #     double dM = 4.0 * M_PI * sum;
    #     _Mv[i] = _Mv[i-1] + dM;
    # }
    # for (int i=0; i<Ns; i++) _Mv[i] /= _Mv[Ns-1];

# -----------------------------------------------------------------

# FROM GalSIM:
# https://github.com/GalSim-developers/GalSim/blob/4a2a281f96d11eab63c7584b34e2598a978e5b86/src/SBSersic.cpp

# // Find what radius encloses (1-missing_flux_frac) of the total flux in a Sersic profile,
# // in units of half-light radius re.
# double SBSersic::SersicInfo::findMaxR(double missing_flux_frac, double gamma2n)
# {
#     // int(exp(-b r^1/n) r, r=R..inf) = x * int(exp(-b r^1/n) r, r=0..inf)
#     //                                = x n b^-2n Gamma(2n)    [x == missing_flux_frac]
#     // Change variables: u = b r^1/n,
#     // du = b/n r^(1-n)/n dr
#     //    = b/n r^1/n dr/r
#     //    = u/n dr/r
#     // r dr = n du r^2 / u
#     //      = n du (u/b)^2n / u
#     // n b^-2n int(u^(2n-1) exp(-u), u=bR^1/n..inf) = x n b^-2n Gamma(2n)
#     // Let z = b R^1/n
#     //
#     // int(u^(2n-1) exp(-u), u=z..inf) = x Gamma(2n)
#     //
#     // The lhs is an incomplete gamma function: Gamma(2n,z), which according to
#     // Abramowitz & Stegun (6.5.32) has a high-z asymptotic form of:
#     // Gamma(2n,z) ~= z^(2n-1) exp(-z) (1 + (2n-2)/z + (2n-2)(2n-3)/z^2 + ... )
#     // ln(x Gamma(2n)) = (2n-1) ln(z) - z + 2(n-1)/z + 2(n-1)(n-2)/z^2
#     // z = -ln(x Gamma(2n) + (2n-1) ln(z) + 2(n-1)/z + 2(n-1)(n-2)/z^2
#     // Iterate this until it converges.  Should be quick.
#     dbg<<"Find maxR for missing_flux_frac = "<<missing_flux_frac<<std::endl;
#     double z0 = -std::log(missing_flux_frac * gamma2n);
#     // Successive approximation method:
#     double z = 4.*(_n+1.);  // A decent starting guess for a range of n.
#     double oldz = 0.;
#     const int MAXIT = 15;
#     dbg<<"Start with z = "<<z<<std::endl;
#     for(int niter=0; niter < MAXIT; ++niter) {
#         oldz = z;
#         z = z0 + (2.*_n-1.) * std::log(z) + 2.*(_n-1.)/z + 2.*(_n-1.)*(_n-2.)/(z*z);
#         dbg<<"z = "<<z<<", dz = "<<z-oldz<<std::endl;
#         if (std::abs(z-oldz) < 0.01) break;
#     }
#     dbg<<"Converged at z = "<<z<<std::endl;
#     double R=std::pow(z/_b, _n);
#     dbg<<"R = (z/b)^n = "<<R<<std::endl;
#     return R;
# }

# -----------------------------------------------------------------

class SersicModel3D(Model3D):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SersicModel3D, self).__init__()

        # Define properties
        self.add_property("effective_radius", "quantity", "effective radius")
        self.add_property("index", "real", "sersic index")
        self.add_property("y_flattening", "real", "flattening along y direction", 1.)
        self.add_property("z_flattening", "real", "flattening along z direction", 1.)
        self.add_property("azimuth", "angle", "azimuth angle", Angle(0.0, "deg"))
        self.add_property("tilt", "angle", "tilt angle", Angle(0.0, "deg"))

        # Truncation
        self.add_property("truncation", "real", "truncation limit (in units of effective_radius)", default_truncation)

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_2d(cls, sersic2d, inclination, position_angle, azimuth_or_tilt="azimuth"):

        """
        :param sersic2d:
        :param inclination:
        :param position_angle:
        :param azimuth_or_tilt:
        :return:
        """

        # Get effective radius and Sersic index
        effective_radius = sersic2d.effective_radius
        index = sersic2d.index

        # Tilt angle and z flattening
        if azimuth_or_tilt == "tilt":

            # Calculate the intrinsic flattening
            y_flattening = 1.
            z_flattening = intrinsic_z_flattening(sersic2d.axial_ratio, inclination)

            # Set azimuth
            azimuth = 0. * u("deg")

            # Calculate the tilt angle of the bulge (tilt w.r.t. the x-axis)
            tilt = deproject_pa_to_tilt(sersic2d.position_angle - position_angle, inclination)
            tilt = Angle(90., "deg") - tilt

        # Azimuth angle and y flattening
        elif azimuth_or_tilt == "azimuth":

            # TODO: this is not working right and therefore unusable for the moment

            # Calculate the intrinsic flattening
            y_flattening = intrinsic_z_flattening(sersic2d.axial_ratio, inclination)
            #z_flattening = intrinsic_z_flattening(sersic2d.axial_ratio, inclination)
            z_flattening = 1.

            # Calculate the azimuth angle of the bulge
            azimuth = deproject_pa_to_azimuth(sersic2d.position_angle - position_angle, inclination)

            # Set tilt
            #tilt = Angle(90., "deg")
            tilt = Angle(0., "deg")

        # Other input
        else: raise ValueError("Incorrect value for 'azimuth_or_tilt'")

        # Create a new Sersic model and return it
        return cls(effective_radius=effective_radius, index=index, y_flattening=y_flattening, z_flattening=z_flattening,
                   azimuth=azimuth, tilt=tilt)

    # -----------------------------------------------------------------

    @property
    def xmin(self):

        """
        This function ...
        :return:
        """

        return - self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def xmax(self):

        """
        This function ...
        :return:
        """

        return self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def xrange(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.xmin, self.xmax)

    # -----------------------------------------------------------------

    @property
    def ymin(self):

        """
        This function ...
        :return:
        """

        return - self.truncation * self.effective_radius * self.y_flattening

    # -----------------------------------------------------------------

    @property
    def ymax(self):

        """
        This function ...
        :return:
        """

        return self.truncation * self.effective_radius * self.y_flattening

    # -----------------------------------------------------------------

    @property
    def yrange(self):

        """
        THis function ...
        :return:
        """

        return QuantityRange(self.ymin, self.ymax)

    # -----------------------------------------------------------------

    @property
    def zmin(self):

        """
        This function ...
        :return:
        """

        return - self.truncation * self.effective_radius * self.z_flattening

    # -----------------------------------------------------------------

    @property
    def zmax(self):

        """
        This function ...
        :return:
        """

        return self.truncation * self.effective_radius * self.z_flattening

    # -----------------------------------------------------------------

    @property
    def zrange(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.zmin, self.zmax)

    # -----------------------------------------------------------------

    @property
    def rho0(self):

        """
        This function ...
        :return:
        """

        # FROM SKIRT

        return 1.0 / self.effective_radius**3

    # -----------------------------------------------------------------

    @property
    def b(self):

        """
        This function ...
        :return:
        """

        # FROM SKIRT

        return 2.0 * self.index - 1.0 / 3.0 + 4.0 / 405.0 / self.index + 46.0 / 25515.0 / self.index**2 + 131.0 / 1148175.0 / (self.index**3)

    # -----------------------------------------------------------------

    def _sersic_function(self, s):

        """
        This function ...
        :param s:
        :return:
        """

        # FROM SKIRT

        # Ns = _sv.size()
        # if s <= _sv[0]: return _Sv[0]
        # elif (s >= _sv[Ns - 1])
        #     return _Sv[Ns - 1]
        # else:
        #     i = NR::locate_clip(_sv, s)
        #     return NR::interpolate_loglog(s, _sv[i], _sv[i + 1], _Sv[i], _Sv[i + 1]);

        pass

    # -----------------------------------------------------------------

    def _spherical_density(self, radius):

        """
        This function ...
        :param radius:
        :return:
        """

        # FROM SKIRT
        #amplitude = 1.
        #s = radius / self.effective_radius
        #return self.rho0 * sersic_function(s) #(*_sersicfunction)(s)

        # ADAPTED FROM 2D SERSIC (ASTROPY):

        #amplitude = self.rho0
        #return amplitude * np.exp(-bn * (z ** (1 / n) - 1))

    # -----------------------------------------------------------------

    def density(self, x, y, height):

        """
        This function ...
        :param x:
        :param y:
        :param height:
        :return:
        """

        # FROM SKIRT

        #m = np.sqrt(radius**2 + height**2 / self.z_flattening**2)
        #return 1.0 / self.z_flattening * self._spherical_density(m)

        # ADAPTED FROM 2D SERSIC (ASTROPY)

        bn = gammaincinv(2. * self.index, 0.5)

        x0 = 0.0
        y0 = 0.0
        z0 = 0.0

        #a, b = self.effective_radius, (1 - ellip) * self.effective_radius
        a, b = self.effective_radius, self.y_flattening * self.effective_radius
        c = self.z_flattening * self.effective_radius

        #print("a", a)
        #print("b", b)
        #print("c", c)

        # TODO: for now, we ignored the tilt angle

        #cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        #x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
        #x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta

        x_maj = x - x0
        x_min = y - y0
        zz = height - z0

        amplitude = self.rho0

        z = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2 + (zz / c) **2)
        #print("z", z)
        return amplitude * np.exp(-bn * (z ** (1 / self.index) - 1))

    # -----------------------------------------------------------------

    def density_function(self, unit="pc", normalize=False):

        """
        This function ...
        :param unit:
        :param normalize:
        :return:
        """

        # Parse the unit
        unit = u(unit)
        density_unit = unit**(-3)

        amplitude = self.rho0.to(density_unit).value

        bn = gammaincinv(2. * self.index, 0.5)

        # Calcualte unitless
        a = self.effective_radius.to(unit).value
        b = (self.y_flattening * self.effective_radius).to(unit).value
        c = (self.z_flattening * self.effective_radius).to(unit).value

        # Define the function
        def sersic(x, y, z):

            """
            This function ...
            :param x:
            :param y:
            :param z:
            :return:
            """

            t = np.sqrt((x / a) ** 2 + (y / b) ** 2 + (z / c) ** 2)
            result = amplitude * np.exp(-bn * (t ** (1 / self.index) - 1))
            if normalize: result /= np.sum(result)
            return result

        # Return the function
        return sersic

# -----------------------------------------------------------------

# IN SKIFILE:
# # Create and add the new exponential disk geometry
#attrs = {"radialScale": str(radial_scale), "axialScale": str(axial_scale), "radialTrunc": str(radial_truncation), "axialTrunc": str(axial_truncation), "innerRadius": str(inner_radius)}

# _Rmax = radialTrunc
# _zmax = axialTrunc
# _Rmin = innerRadius

# -----------------------------------------------------------------

class ExponentialDiskModel3D(Model3D):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ExponentialDiskModel3D, self).__init__()

        # Define properties
        self.add_property("radial_scale", "quantity", "radial scale")
        self.add_property("axial_scale", "quantity", "axial scale")
        self.add_property("radial_truncation", "real", "radial truncation", 0)
        self.add_property("axial_truncation", "real", "axial truncation", 0)
        self.add_property("inner_radius", "real", "inner radius", 0)
        self.add_property("tilt", "angle", "tilt", Angle(0.0, "deg"))

        # Truncation
        self.add_property("truncation", "real", "truncation limit (in units of the radial and axial scales", default_truncation)

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_2d(cls, exponentialdiskmodel2d, inclination, position_angle):

        """
        This function ...
        :param exponentialdiskmodel2d:
        :param inclination:
        :param position_angle:
        :return:
        """

        # Get the radial scale
        radial_scale = exponentialdiskmodel2d.scalelength

        # Calculate the intrinsic flattening
        flattening = intrinsic_z_flattening(exponentialdiskmodel2d.axial_ratio, inclination)

        # Calculate the axial scale
        axial_scale = flattening * radial_scale

        # Calculate the tilt angle of the disk (tilt w.r.t. the x-axis)
        #tilt = deproject_pa_to_tilt(parameters.PA - position_angle, inclination)
        tilt = deproject_pa_to_tilt(position_angle - exponentialdiskmodel2d.position_angle, inclination)
        tilt = Angle(90., "deg") - tilt

        # Create a new exponential disk model and return it
        return cls(radial_scale=radial_scale, axial_scale=axial_scale, tilt=tilt)

    # -----------------------------------------------------------------

    @property
    def xmin(self):

        """
        This function ...
        :return:
        """

        if self.radial_truncation > 0: return - self.radial_truncation
        else: return - self.truncation * self.radial_scale

    # -----------------------------------------------------------------

    @property
    def xmax(self):

        """
        This function ...
        :return:
        """

        if self.radial_truncation > 0: return self.radial_truncation
        else: return self.truncation * self.radial_scale

    # -----------------------------------------------------------------

    @property
    def xrange(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.xmin, self.xmax)

    # -----------------------------------------------------------------

    @property
    def ymin(self):

        """
        This function ...
        :return:
        """

        if self.radial_truncation > 0: return - self.radial_truncation
        else: return - self.truncation * self.radial_scale

    # -----------------------------------------------------------------

    @property
    def ymax(self):

        """
        This function ...
        :return:
        """

        if self.radial_truncation > 0: return self.radial_truncation
        else: return self.truncation * self.radial_scale

    # -----------------------------------------------------------------

    @property
    def yrange(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.ymin, self.ymax)

    # -----------------------------------------------------------------

    @property
    def zmin(self):

        """
        This function ...
        :return:
        """

        if self.axial_truncation > 0: return - self.axial_truncation
        else: return - self.truncation * self.axial_scale

    # -----------------------------------------------------------------

    @property
    def zmax(self):

        """
        This function ...
        :return:
        """

        if self.axial_truncation > 0: return self.axial_truncation
        else: return self.truncation * self.axial_scale

    # -----------------------------------------------------------------

    @property
    def zrange(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.zmin, self.zmax)

    # -----------------------------------------------------------------

    @property
    def intz(self):

        """
        This function ...
        :return: 
        """

        if self.axial_truncation > 0: return -2.0 * self.axial_scale * (np.exp(- self.axial_truncation / self.axial_scale) - 1.)
        else: return 2.0 * self.axial_scale

    # -----------------------------------------------------------------

    @property
    def tmin(self):

        """
        This function ...
        :return:
        """

        if self.inner_radius > 0: return np.exp(- self.inner_radius / self.radial_scale) * (1.0 + self.inner_radius / self.radial_scale)
        else: return 1.0

    # -----------------------------------------------------------------

    @property
    def tmax(self):

        """
        This function ...
        :return:
        """

        if self.radial_truncation > 0: return np.exp(- self.radial_truncation / self.radial_scale) * (1.0 + self.radial_truncation / self.radial_scale)
        else: return 0.0

    # -----------------------------------------------------------------

    @property
    def intR(self):

        """
        Thisf unction ...
        :return:
        """

        return self.radial_scale**2 * (self.tmin - self.tmax)

    # -----------------------------------------------------------------

    @property
    def intphi(self):

        """
        This function ...
        :return:
        """

        return 2.0 * np.pi

    # -----------------------------------------------------------------

    @property
    def rho0(self):

        """
        This function ...
        :return:
        """

        return 1.0 / (self.intR * self.intphi * self.intz)

    # -----------------------------------------------------------------

    def surface_density(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        radius = np.sqrt(x ** 2 + y ** 2)

        # Checks
        if self.radial_truncation > 0.0 and radius > self.radial_truncation: return 0.0
        if radius < self.inner_radius: return 0.0

        # Return the surface density
        return self.rho0 * np.exp(- radius / self.radial_scale)

    # -----------------------------------------------------------------

    def density(self, x, y, z):

        """
        This function ...
        :param x:
        :param y:
        :param z:
        :return:
        """

        absz = abs(z)

        # Check
        if self.axial_truncation > 0.0 and absz > self.axial_truncation: return 0.0

        # Return the density
        return self.surface_density(x, y) * np.exp(- absz / self.axial_scale)

    # -----------------------------------------------------------------

    def surface_density_function(self, unit="pc", normalize=False):

        """
        This function ...
        :param unit:
        :param normalize:
        :return:
        """

        # Parse the unit
        unit = u(unit)
        density_unit = unit ** -3

        amplitude = self.rho0.to(density_unit).value
        radialtrunc = self.radial_truncation.to(unit).value if self.radial_truncation > 0.0 else None
        innerradius = self.inner_radius.to(unit).value if self.inner_radius > 0.0 else None

        radialscale = self.radial_scale.to(unit).value

        # Define the function
        def exponentialdisk(x, y):

            radius = np.sqrt(x ** 2 + y ** 2)

            density = amplitude * np.exp(- radius / radialscale)

            if radialtrunc is not None: density[radius > radialtrunc] = 0.0
            if innerradius is not None: density[radius < innerradius] = 0.0

            if normalize: density /= np.sum(density)
            return density

        # Return the function
        return exponentialdisk

    # -----------------------------------------------------------------

    def density_function(self, unit="pc", normalize=False):

        """
        This function ...
        :param unit:
        :param normalize:
        :return:
        """

        # Parse the unit
        unit = u(unit)
        density_unit = unit**-3

        amplitude = self.rho0.to(density_unit).value
        radialtrunc = self.radial_truncation.to(unit).value if self.radial_truncation > 0.0 else None
        axialtrunc = self.axial_truncation.to(unit).value if self.axial_truncation > 0.0 else None
        innerradius = self.inner_radius.to(unit).value if self.inner_radius > 0.0 else None

        radialscale = self.radial_scale.to(unit).value
        axialscale = self.axial_scale.to(unit).value

        # Define the function
        def exponentialdisk(x, y, z):

            """
            This function ...
            :param x:
            :param y:
            :param z:
            :return:
            """

            radius = np.sqrt(x ** 2 + y ** 2)

            height = abs(z)
            density = amplitude * np.exp(- radius / radialscale) * np.exp(- height / axialscale)

            if radialtrunc is not None: density[radius > radialtrunc] = 0.0
            if axialtrunc is not None: density[height > axialtrunc] = 0.0
            if innerradius is not None: density[radius < innerradius] = 0.0

            if normalize: density /= np.sum(density)
            return density

        # Return the function
        return exponentialdisk

# -----------------------------------------------------------------

class RingModel3D(Model):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RingModel3D, self).__init__()

        # Define the properties
        self.add_property("radius", "quantity", "radius of the ring")
        self.add_property("width", "quantity", "width of the ring")
        self.add_property("height", "quantity", "height of the ring")

        # Set the properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class DeprojectionModel3D(Model3D):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(DeprojectionModel3D, self).__init__()

        # position angle: -360 deg to 360 deg
        # inclination: 0 deg to 90 deg
        # center in pixel coordinates!

        # Define the properties
        self.add_property("filename", "string", "name of the input FITS file")
        self.add_property("pixelscale", "quantity", "physical pixelscale of the FITS image")
        self.add_property("position_angle", "angle", "position angle")
        self.add_property("inclination", "angle", "inclination")
        self.add_property("x_size", "positive_integer", "number of x pixels in the map")
        self.add_property("y_size", "positive_integer", "number of y pixels in the map")
        self.add_property("x_center", "real", "x center in image coordinates")
        self.add_property("y_center", "real", "y center in image coordinates")
        self.add_property("scale_height", "quantity", "scale height")

        # Path of the directory containing the map
        self.add_property("dirpath", "string", "path to the directory where the map may be present", fs.cwd())

        # Distance (if not defined in the map)
        self.add_property("distance", "quantity", "distance to the galaxy")

        # Truncation
        self.add_property("axial_truncation", "real", "truncation in units of the scale height", default_truncation)

        # Truncate zero values around the luminosity of the map
        self.add_property("truncate_zeros", "boolean", "truncate in the xy plane around where the map is above zero", True)

        # Set the properties
        self.set_properties(kwargs)

        # The map that is loaded
        self._map = None

    # -----------------------------------------------------------------

    @classmethod
    def from_wcs(cls, wcs, galaxy_center, distance, pa, inclination, filepath, scale_height):

        """
        This function ...
        :param wcs:
        :param galaxy_center:
        :param distance:
        :param pa:
        :param inclination:
        :param filepath
        :param scale_height:
        :return:
        """

        # Get the center pixel
        pixel_center = galaxy_center.to_pixel(wcs)
        xc = pixel_center.x
        yc = pixel_center.y

        # Get the pixelscale in physical units
        pixelscale_angular = wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the number of x and y pixels
        x_size = wcs.xsize
        y_size = wcs.ysize

        # Create the deprojection model
        deprojection = cls(filename=filepath, pixelscale=pixelscale, position_angle=pa, inclination=inclination,
                           x_size=x_size, y_size=y_size, x_center=xc, y_center=yc, scale_height=scale_height, distance=distance)

        # Return the deprojection model
        return deprojection

    # -----------------------------------------------------------------

    @property
    def center(self):
        return PixelCoordinate(x=self.x_center, y=self.y_center)

    # -----------------------------------------------------------------

    @property
    def filepath(self):
        return fs.absolute_or_in(self.filename, self.dirpath)

    # -----------------------------------------------------------------

    @property
    def has_map(self):

        """
        has_map == True means that the map can be accessed via self.map
        :return:
        """

        return self._map is not None or fs.is_file(self.filepath)

    # -----------------------------------------------------------------

    @property
    def map_is_loaded(self):
        return self._map is not None

    # -----------------------------------------------------------------

    @property
    def map(self):

        """
        This function ...
        :return:
        """

        if self._map is not None: return self._map
        elif fs.is_file(self.filepath):

            # Load the frame
            frame = Frame.from_file(self.filepath)

            # Verify the image properties
            if self.x_size != frame.xsize: raise ValueError("Number of x pixels does not correspond with the number of x pixels of the image")
            if self.y_size != frame.ysize: raise ValueError("Number of y pixels does not correspond with the number of y pixels of the image")

            # Normalize the map
            frame.normalize()

            # Set the map
            self._map = frame

            # Return the map
            return self._map

        # Error
        else: raise ValueError("Map cannot be loaded from '" + self.filepath + "'")

    # -----------------------------------------------------------------

    @map.setter
    def map(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Make copy
        frame = value.copy()

        # Verify properties
        if self.x_size != frame.xsize: raise ValueError("Number of x pixels does not correspond with the number of x pixels of the image")
        if self.y_size != frame.ysize: raise ValueError("Number of y pixels does not correspond with the number of y pixels of the image")

        # Normalize
        frame.normalize()

        # Set the map
        self._map = frame

    # -----------------------------------------------------------------

    @lazyproperty
    def map_zeros(self):
        return self.map.zeroes

    # -----------------------------------------------------------------

    @lazyproperty
    def map_nonzeroes(self):
        return self.map.nonzeroes

    # -----------------------------------------------------------------

    @lazyproperty
    def map_zeros_pixels(self):
        return self.map.zeroes_pixels

    # -----------------------------------------------------------------

    @lazyproperty
    def map_nonzeroes_pixels(self):
        return self.map.nonzeroes_pixels

    # -----------------------------------------------------------------

    @lazyproperty
    def map_nonzeroes_x(self):
        return self.map.nonzeroes_x

    # -----------------------------------------------------------------

    @lazyproperty
    def map_nonzeroes_y(self):
        return self.map.nonzeroes_y

    # -----------------------------------------------------------------

    @lazyproperty
    def min_nonzero_x(self):
        return np.min(self.map_nonzeroes_x)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_nonzero_x(self):
        return np.max(self.map_nonzeroes_x)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_nonzero_y(self):
        return np.min(self.map_nonzeroes_y)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_nonzero_y(self):
        return np.max(self.map_nonzeroes_y)

    # -----------------------------------------------------------------

    @lazyproperty
    def ntruncated_min_x(self):
        if self.truncate_zeros: return self.min_nonzero_x
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def ntruncated_max_x(self):
        if self.truncate_zeros: return self.xsize - self.max_nonzero_x
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def ntruncated_min_y(self):
        if self.truncate_zeros: return self.min_nonzero_y
        else: return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def ntruncated_max_y(self):
        if self.truncate_zeros: return self.ysize - self.max_nonzero_y
        else: return 0.0

    # -----------------------------------------------------------------

    @property
    def xsize(self):
        ## NO: don't return the map xsize, return the deprojected xsize? -> NO: not needed so far
        assert self.map.xsize == self.x_size
        return self.map.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):
        ## NO: don't return the map ysize, return the deprojected ysize? -> NO: not needed so far
        assert self.map.ysize == self.y_size
        return self.map.ysize

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):
        if self.distance is not None: return self.distance
        else: return self.map.distance

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_galaxy_distance(self):
        return self.has_distance or (self.has_map and self.map.distance is not None)

    # -----------------------------------------------------------------

    @property
    def angular_pixelscale(self):

        """
        This function ....
        :return:
        """

        # Check whether the distance is defined
        if not self.has_galaxy_distance: raise ValueError("Distance of the map is not defined")

        # Calculate the pixelscale in degrees
        pixelscale = self.pixelscale / self.galaxy_distance * u("rad")
        return Pixelscale(pixelscale.to("arcsec"))

    # -----------------------------------------------------------------

    @lazyproperty
    def projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the distance is defined
        if not self.has_galaxy_distance: raise ValueError("Distance of the map is not defined")

        azimuth = 0.0

        # Create the 'earth' projection system
        return GalaxyProjection.from_deprojection(self, self.galaxy_distance, azimuth)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the distance is defined
        if self.galaxy_distance is None: raise ValueError("Distance of the map is not defined")

        # Create the face-on projection system
        return FaceOnProjection.from_deprojection(self, self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the distance is defined
        if self.galaxy_distance is None: raise ValueError("Distance of the map is not defined")

        # Create the edge-on projection system
        return EdgeOnProjection.from_deprojection(self, self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument(self):
        return FrameInstrument.from_projection(self.projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_instrument(self):
        return SimpleInstrument.from_projection(self.projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):
        return SEDInstrument.from_projection(self.projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument(self):
        return FullInstrument.from_projection(self.projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_frame_instrument(self):
        return FrameInstrument.from_projection(self.faceon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_simple_instrument(self):
        return SimpleInstrument.from_projection(self.faceon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_sed_instrument(self):
        return SEDInstrument.from_projection(self.faceon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_full_instrument(self):
        return FullInstrument.from_projection(self.faceon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_frame_instrument(self):
        return FrameInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_simple_instrument(self):
        return SimpleInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_sed_instrument(self):
        return SEDInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_full_instrument(self):
        return FullInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @property
    def xmax(self):
        #print("ntruncated max x: " + str(self.ntruncated_max_x) + " of " + str(self.xsize))
        return (self.x_size - self.x_center - self.ntruncated_max_x) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def map_xmax(self):
        return (self.x_size - self.x_center) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def xmin(self):
        #print("ntruncated min x: " + str(self.ntruncated_min_x) + " of " + str(self.xsize))
        return - (self.x_center - self.ntruncated_min_x) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def map_xmin(self):
        return - self.x_center * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def xrange(self):
        return QuantityRange(self.xmin, self.xmax)

    # -----------------------------------------------------------------

    @property
    def ymax(self):
        #print("ntruncated max y: " + str(self.ntruncated_max_y) + " of " + str(self.ysize))
        return (self.y_size - self.y_center - self.ntruncated_max_y) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def map_ymax(self):
        return (self.y_size - self.y_center) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def ymin(self):
        #print("ntruncated min y: " + str(self.ntruncated_min_y) + " of " + str(self.ysize))
        return - (self.y_center - self.ntruncated_min_y) * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def map_ymin(self):
        return - self.y_center * self.pixelscale

    # -----------------------------------------------------------------

    @property
    def yrange(self):
        return QuantityRange(self.ymin, self.ymax)

    # -----------------------------------------------------------------

    @property
    def zmin(self):
        return - self.axial_truncation * self.scale_height

    # -----------------------------------------------------------------

    @property
    def zmax(self):
        return self.axial_truncation * self.scale_height

    # -----------------------------------------------------------------

    @property
    def zrange(self):
        return QuantityRange(self.zmin, self.zmax)

    # -----------------------------------------------------------------

    @property
    def position_angle_radians(self):
        return self.position_angle.to("rad").value

    # -----------------------------------------------------------------

    @property
    def inclination_radians(self):
        return self.inclination.to("rad").value

    # -----------------------------------------------------------------

    @property
    def cospa(self):
        # Calculate the sines and cosines of the position angle and inclination
        return np.cos(self.position_angle_radians)

    # -----------------------------------------------------------------

    @property
    def sinpa(self):
        return np.sin(self.position_angle_radians)

    # -----------------------------------------------------------------

    @property
    def cosi(self):
        return np.cos(self.inclination_radians)

    # -----------------------------------------------------------------

    @property
    def sini(self):
        return np.sin(self.inclination_radians)

    # -----------------------------------------------------------------

    @property
    def deltay(self):
        return self.pixelscale

    # -----------------------------------------------------------------

    @property
    def deltax(self):
        # Calculate the physical pixel size in the x direction of the galactic plane
        return self.pixelscale / self.cosi

    # -----------------------------------------------------------------

    @lazyproperty
    def corner1(self):
        return self.derotate(self.xmax, self.ymax)

    # -----------------------------------------------------------------

    @lazyproperty
    def corner2(self):
        return self.derotate(self.xmin, self.ymax)

    # -----------------------------------------------------------------

    @lazyproperty
    def corner3(self):
        return self.derotate(self.xmin, self.ymin)

    # -----------------------------------------------------------------

    @lazyproperty
    def corner4(self):
        return self.derotate(self.xmax, self.ymin)

    # -----------------------------------------------------------------

    def rotate(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        """

        # Cache the original values of x and y
        xorig = x
        yorig = y

        # Calculate the coordinates in the plane of the image
        x = self.sinpa * xorig  + self.cospa * yorig
        y = -self.cospa * xorig + self.sinpa * yorig

        # Return
        return x, y

    # -----------------------------------------------------------------

    def rotate_arrays(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        xorig = np.copy(x)
        yorig = np.copy(y)

        # Calculate the coordinates in the plane of the image
        x = self.sinpa * xorig + self.cospa * yorig
        y = -self.cospa * xorig + self.sinpa * yorig

        # Return
        return x, y

    # -----------------------------------------------------------------

    def derotate(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Cache the original values of x and y
        xorig = x
        yorig = y

        #// Calculate the coordinates in the rotated plane
        x = (self.sinpa * xorig) - (self.cospa * yorig)
        y = (self.cospa * xorig) + (self.sinpa * yorig)

        # Return
        return x, y

    # -----------------------------------------------------------------

    def project(self, x):

        """
        This function ...
        :return:
        """

        return x * self.cosi

    # -----------------------------------------------------------------

    def project_array(self, x):

        """
        This function ...
        :param x:
        :return:
        """

        return x * self.cosi

    # -----------------------------------------------------------------

    def deproject(self, x):

        """
        This function ...
        :param x:
        :return:
        """

        return x / self.cosi

    # -----------------------------------------------------------------

    def surface_density(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Project and rotate the x and y coordinates
        x = self.project(x)
        x, y = self.rotate(x, y)

        # Find the corresponding pixel in the image
        i = numbers.round_down_to_int((x - self.map_xmin) / self.deltay)
        j = numbers.round_down_to_int((y - self.map_ymin) / self.deltay)

        # Not on the image
        if i < 0 or i >= self.x_size or j < 0 or j >= self.y_size: return 0.0

        # Return the density
        return self.map[j, i]

    # -----------------------------------------------------------------

    @property
    def density_normalization(self):
        return 1. / (2. * self.scale_height) / (self.deltax * self.deltay)

    # -----------------------------------------------------------------

    def density(self, x, y, z):

        """
        This function ...
        :param x:
        :param y:
        :param z:
        :return:
        """

        # Return the density
        z = abs(z)
        return self.surface_density(x, y) * np.exp(- z / self.scale_height) * self.density_normalization

    # -----------------------------------------------------------------

    @property
    def shape(self):
        return PixelShape(self.ysize, self.xsize)

    # -----------------------------------------------------------------

    @property
    def npixels(self):
        return self.shape.ntotalpixels

    # -----------------------------------------------------------------

    def surface_density_function(self, unit="pc", normalize=False):

        """
        This function ...
        :param unit:
        :param normalize:
        :return:
        """

        # Get properties as scalars (no units)
        xmin_scalar = self.map_xmin.to(unit).value
        ymin_scalar = self.map_ymin.to(unit).value
        deltay_scalar = self.deltay.to(unit).value
        deltax_scalar = self.deltax.to(unit).value

        # Define the density function
        def deprojection(x, y):

            # Project and rotate coordinates
            x = self.project_array(x)
            x, y = self.rotate_arrays(x, y)

            # Determine the coordinate mapping
            x_mapping = ((x - xmin_scalar) / deltay_scalar - 0.5).astype(int)
            y_mapping = ((y - ymin_scalar) / deltay_scalar - 0.5).astype(int)

            xy = x_mapping + self.xsize * y_mapping

            nx = xy.shape[0]
            ny = xy.shape[1]

            #xy = xy[:, :, 0]

            arrangement_x = range(self.xsize) * self.ysize
            lists = [[value] * self.xsize for value in range(self.ysize)]

            arrangement_y = list(itertools.chain.from_iterable(lists))

            xy_flattened = xy.flatten()

            x_mapping = [arrangement_x[k] if 0 <= k < self.npixels else -1 for k in xy_flattened]
            y_mapping = [arrangement_y[k] if 0 <= k < self.npixels else -1 for k in xy_flattened]
            mapping = [y_mapping, x_mapping]

            cval = 0.0
            data = ndimage.map_coordinates(self.map.data, mapping, mode='constant', cval=cval)

            # Reshape into the
            deprojected = data.reshape((ny, nx))

            # Normalize?
            if normalize: deprojected /= np.sum(deprojected)
            return deprojected

        # Return the function
        return deprojection

    # -----------------------------------------------------------------

    def density_function(self, unit="pc", normalize=False):

        """
        This function ...
        :param unit:
        :param normalize:
        :return:
        """

        # Get properties as scalars (no units)
        xmin_scalar = self.map_xmin.to(unit).value
        ymin_scalar = self.map_ymin.to(unit).value
        deltay_scalar = self.deltay.to(unit).value
        deltax_scalar = self.deltax.to(unit).value
        scale_height_scalar = self.scale_height.to(unit).value

        # Define the density function
        def deprojection(x, y, z):

            """
            Thisf unction ...
            :param x:
            :param y:
            :param z:
            :return:
            """

            # Project and rotate coordinates
            x = self.project_array(x)
            x, y = self.rotate_arrays(x, y)

            #i_array = numbers.round_down_to_int((x - xmin_scalar) / deltay_scalar)
            #j_array = numbers.round_down_to_int((y - ymin_scalar) / deltay_scalar)

            # grid = mapping
            #grid = np.array([yy1.reshape(outshape), xx1.reshape(outshape)])
            # Use Scipy to create the new image
            #data = scipy.ndimage.map_coordinates(frame, mapping, **kwargs)

            # Determine the coordinate mapping
            #x_mapping = (x - xmin_scalar) / deltay_scalar - 0.5
            #y_mapping = (y - ymin_scalar) / deltay_scalar - 0.5
            x_mapping = ((x - xmin_scalar) / deltay_scalar - 0.5).astype(int)
            y_mapping = ((y - ymin_scalar) / deltay_scalar - 0.5).astype(int)

            #x_mapping = x_mapping.flatten()
            #y_mapping = y_mapping.flatten()

            xy = x_mapping + self.xsize * y_mapping

            #print("xy", xy)

            #print("x_mapping", x_mapping)
            #print("y_mapping", y_mapping)

            #mapping = np.array([np.array([xi, yi]) for xi, yi in zip(x_mapping.flat, y_mapping.flat)])
            #mapping = np.array([x_mapping, y_mapping])
            #mapping = np.vstack((x_mapping,y_mapping))

            #print("mapping", mapping)
            #print(mapping.shape)

            #print("output shape", mapping.shape[1:])
            #print("input.ndim", self.map.data.ndim)

            # output_shape = coordinates.shape[1:]
            # input.ndim < 1 or len(output_shape) < 1:

            #print(xy)
            #print(xy.shape)

            nx = xy.shape[0]
            ny = xy.shape[1]
            #nz = xy.shape[2]

            #print(xy[:,:,0])
            xy = xy[:, :, 0]
            #print(xy)

            arrangement_x = range(self.xsize) * self.ysize
            lists = [[value] * self.xsize for value in range(self.ysize)]
            #print(lists)
            #arrangement_y = itertools.chain(*lists)
            arrangement_y = list(itertools.chain.from_iterable(lists))

            #print(xy.shape)

            #print("output shape", xy.shape[1:])
            #print("input.ndim", self.map.data.flatten().ndim)

            # Create mapping
            #mapping = [[arrangement_y[yi], arrangement_x[xi]] for xi, yi in xy]
            #print("mapping", mapping)

            #print([k for k in xy])

            xy_flattened = xy.flatten()
            #print(arrangement_x)
            #print(arrangement_y.shape)
            x_mapping = [arrangement_x[k] if 0<=k<self.npixels else -1 for k in xy_flattened]
            y_mapping = [arrangement_y[k] if 0<=k<self.npixels else -1 for k in xy_flattened]
            #mapping = [x_mapping, y_mapping]
            mapping = [y_mapping, x_mapping]

            #cval = float('nan')
            cval = 0.0
            data = ndimage.map_coordinates(self.map.data, mapping, mode='constant', cval=cval)
            #data = ndimage.map_coordinates(self.map.data.flatten(), xy, mode="constant", cval=cval)

            #nx = x.shape[2]
            #ny = y.shape[1]

            #print("nx", nx)
            #print("")

            # Reshape into the
            deprojected = data.reshape((ny, nx, 1))

            # Return
            z = abs(z)
            result = deprojected * np.exp(- z / scale_height_scalar) / (2. * scale_height_scalar) / (deltax_scalar * deltay_scalar)

            # Normalize?
            if normalize: result /= np.sum(result)
            return result

        # Return the function
        return deprojection

# -----------------------------------------------------------------

def intrinsic_z_flattening(qprime, inclination):

    """
    This function ...
    :param qprime:
    :param inclination:
    """

    # Get the inclination angle in radians
    i = inclination.to("radian").value

    # Calculate the intrinsic flattening
    difference = qprime**2 - math.cos(i)**2

    if difference < 0:
        # Give warning
        warnings.warn("Could not convert an apparent flattening of " + str(qprime) + " with an inclination of " + str(inclination) + " to a vertical flattening: using the apparent flattening")
        q = qprime
    else: q = math.sqrt(difference/math.sin(i)**2)

    # Return the intrinsic flattening
    return q

# -----------------------------------------------------------------

def intrinsic_y_flattening(qprime, inclination):

    """
    This function ...
    :param qprime:
    :param inclination:
    :return:
    """

    # Get the inclination angle in radians
    i = inclination.to("radian").value

    # Calculate the 'inclination' w.r.t. the y axis
    #i_wrt_y = 0.5 * math.pi - i

    print(inclination)
    print(qprime)
    print(i_wrt_y)
    print(math.cos(i_wrt_y))
    print(math.sin(i_wrt_y))

    #qprime =

    # Calculate the intrinsic flattening
    q = math.sqrt((qprime ** 2 - math.cos(i_wrt_y) ** 2) / math.sin(i_wrt_y) ** 2)

    # Return the intrinsic flattening
    return q

# -----------------------------------------------------------------

def project_azimuth_to_pa(azimuth, inclination):

    """
    This function ...
    :param azimuth:
    :param inclination:
    """

    # Get the azimuth angle and inclination in radians
    azimuth_radian = azimuth.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.cos(azimuth_radian)**2 * math.cos(i_radian)**2 + math.sin(azimuth_radian)**2)

    cos_pa = math.cos(azimuth_radian) * math.cos(i_radian) / denominator
    sin_pa = math.sin(azimuth_radian) / denominator

    pa_radian = math.atan2(sin_pa, cos_pa) * u("radian")

    return pa_radian.to("deg")

# -----------------------------------------------------------------

def deproject_pa_to_azimuth(pa, inclination):

    """
    This function ...
    :param pa:
    :param inclination:
    :return:
    """

    # Get the PA and inclination in radians
    pa_radian = pa.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.cos(pa_radian)**2 + math.sin(pa_radian)**2 * math.cos(i_radian)**2)

    cos_azimuth = math.cos(pa_radian) / denominator
    sin_azimuth = math.sin(pa_radian) * math.cos(i_radian) / denominator

    azimuth_radian = math.atan2(sin_azimuth, cos_azimuth) * u("radian")

    return azimuth_radian.to("deg")

# -----------------------------------------------------------------

def project_tilt_to_pa(tilt, inclination):

    """
    This function ...
    :param tilt:
    :param inclination:
    :return:
    """

    # Get the tilt angle and inclination in radians
    tilt_radian = tilt.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.sin(tilt_radian)**2 * math.sin(i_radian)**2 + math.cos(tilt_radian)**2)

    cos_pa = math.sin(tilt_radian) * math.sin(i_radian) / denominator
    sin_pa = math.cos(tilt_radian) / denominator

    pa_radian = math.atan2(sin_pa, cos_pa) * u("radian")

    return pa_radian.to("deg")

# -----------------------------------------------------------------

def deproject_pa_to_tilt(pa, inclination):

    """
    This function ...
    :param pa:
    :param inclination:
    :return:
    """

    # Get the PA and inclination in radians
    pa_radian = pa.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.sin(pa_radian)**2 * math.sin(i_radian)**2 + math.cos(pa_radian)**2)

    cos_tilt = math.sin(pa_radian) * math.sin(i_radian) / denominator
    sin_tilt = math.cos(pa_radian) / denominator

    tilt_radian = math.atan2(sin_tilt, cos_tilt) * u("radian")

    return tilt_radian.to("deg")

# -----------------------------------------------------------------

# Test the deprojection functions:

#test_angles = []
#test_angles.append(Angle(33, "deg"))
#test_angles.append(Angle(45, "deg"))
#test_angles.append(Angle(0, "deg"))
#test_angles.append(Angle(90, "deg"))
#test_angles.append(Angle(189, "deg"))

#for test_angle in test_angles:

#    result = deproject_pa_to_tilt(test_angle, Angle(0.0, "deg"))
#    print("Should fail/be zero:", Angle(90., "deg") - result, test_angle)

#    result = deproject_pa_to_tilt(test_angle, Angle(90., "deg"))
#    print("Should be the same:", Angle(90., "deg") - result, test_angle)

    #result = project_tilt_to_pa(test_angle, Angle(0.0, "deg"))
    #print("Should fail/be zero:", Angle(90., "deg") - result, test_angle)

    #result = project_tilt_to_pa(test_angle, Angle(90., "deg"))
    #print("Should be the same:", Angle(90., "deg") - result, test_angle)

# -----------------------------------------------------------------
#
# 2D MODELS
#
# -----------------------------------------------------------------

def load_2d_model(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create the appropriate model
    if "SersicModel2D" in first_line: return SersicModel2D.from_file(path)
    elif "ExponentialDiskModel2D" in first_line: return ExponentialDiskModel2D.from_file(path)
    else: raise ValueError("Unrecognized model file")

# -----------------------------------------------------------------

class SersicModel2D(Model2D):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SersicModel2D, self).__init__()

        # Define properties
        self.add_property("rel_contribution", "real", "relative contribution")
        self.add_property("fluxdensity", "quantity", "flux density")
        self.add_property("axial_ratio", "real", "axial ratio")
        self.add_property("position_angle", "angle", "position angle") # (degrees ccw from North)
        self.add_property("effective_radius", "quantity", "effective radius")
        self.add_property("index", "real", "sersic index")

        # Truncation
        self.add_property("truncation", "real", "truncation radius relative to the effective radius", default_truncation)

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @property
    def xmin(self):
        return - self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def xmax(self):
        return self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def xrange(self):
        return QuantityRange(self.xmin, self.xmax)

    # -----------------------------------------------------------------

    @property
    def ymin(self):
        return - self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def ymax(self):
        return self.truncation * self.effective_radius

    # -----------------------------------------------------------------

    @property
    def yrange(self):
        return QuantityRange(self.ymin, self.ymax)

    # -----------------------------------------------------------------

    @property
    def position_angle_radian(self):
        return self.position_angle.to("rad").value

    # -----------------------------------------------------------------

    def density(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # FROM ASTROPY, SERSIC2D FUNCTION

        bn = gammaincinv(2. * self.index, 0.5)

        # 1 - ELLIP = AXIAL RATIO?

        x_0 = y_0 = 0.0

        #a, b = self.effective_radius, (1. - ellip) * self.effective_radius
        a, b = self.effective_radius, self.axial_ratio * self.effective_radius
        cos_theta, sin_theta = np.cos(self.position_angle_radian), np.sin(self.position_angle_radian)
        x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
        x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
        z = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)

        amplitude = 1.

        return amplitude * np.exp( -bn * (z ** (1 / self.index) - 1))

# -----------------------------------------------------------------

class ExponentialDiskModel2D(Model2D):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ExponentialDiskModel2D, self).__init__()

        # Define properties
        self.add_property("rel_contribution", "real", "relative contribution")
        self.add_property("fluxdensity", "quantity", "flux density")
        self.add_property("axial_ratio", "real", "axial ratio")
        self.add_property("position_angle", "angle", "position_angle") # (degrees ccw from North)
        self.add_property("mu0", "quantity", "surface brightness at center")
        self.add_property("scalelength", "quantity", "scale length")

        # The truncation
        self.add_property("truncation", "real", "truncation length in units of the scalelength", default_truncation)

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @property
    def xmin(self):
        return - self.truncation * self.scalelength

    # -----------------------------------------------------------------

    @property
    def xmax(self):
        return self.truncation * self.scalelength

    # -----------------------------------------------------------------

    @property
    def xrange(self):
        return QuantityRange(self.xmin, self.xmax)

    # -----------------------------------------------------------------

    @property
    def ymin(self):
        return - self.truncation * self.scalelength

    # -----------------------------------------------------------------

    @property
    def ymax(self):
        return self.truncation * self.scalelength

    # -----------------------------------------------------------------

    @property
    def yrange(self):
        return QuantityRange(self.ymin, self.ymax)

# -----------------------------------------------------------------
