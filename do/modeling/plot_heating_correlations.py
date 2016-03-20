#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_heating_correlations

# -----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import matplotlib as mpl
from scipy import optimize
rc('text', usetex=True)
from scipy.stats import gaussian_kde
import math

def main():

    path  = "modelChecks/"
    Lsun = 3.846e26 # Watts
    
    input = np.loadtxt(path+"iteration4_J14/Ultimatrix")
    ID        = input[:,0]      # pixel number
    ra        = input[:,1]
    dec       = input[:,2]
    sSFR      = input[:,62]     # fundamental physical parameters
    sSFR_16   = input[:,63]
    sSFR_84   = input[:,64]
    Mstars    = input[:,65]
    Mstars_16 = input[:,66]
    Mstars_84 = input[:,67]
    Ldust     = input[:,68]
    Ldust_16  = input[:,69]
    Ldust_84  = input[:,70]
    Tw_BC     = input[:,71]
    Tw_BC_16  = input[:,72]
    Tw_BC_84  = input[:,73]
    Tc_ISM    = input[:,74]
    Tc_ISM_16 = input[:,75]
    Tc_ISM_84 = input[:,76]
    xi_PAH    = input[:,80]
    xi_PAH_16 = input[:,80]
    xi_PAH_84 = input[:,80]
    Mdust     = input[:,92]
    Mdust_16  = input[:,93]
    Mdust_84  = input[:,94]
    SFR       = input[:,95]
    SFR_16    = input[:,96]
    SFR_84    = input[:,97]
    F500      = input[:,47]     # observed fluxes
    eF500     = input[:,48]
    F350      = input[:,45]
    eF350     = input[:,46]
    F250      = input[:,43]
    eF250     = input[:,44]
    F160      = input[:,42]
    eF160     = input[:,41]
    F100      = input[:,39]
    eF100     = input[:,40]
    F70       = input[:,37]
    eF70      = input[:,38]
    Fr        = input[:,13]
    eFr       = input[:,14]
    FNUV      = input[:,7]
    eFNUV     = input[:,8]
    
    # Dust heating parameters
    input = np.loadtxt(path+"iteration5_J14/pixelHeating.dat")
    Fold    = input[:,1]
    Fyoung  = input[:,2]
    Lold    = input[:,3]
    Lyoung  = input[:,4]
    Ltot    = input[:,5]
    Lstar   = input[:,6]
    
    # compute derived quantities
    LPAH = Ldust*xi_PAH
    MdMs = Mdust/Mstars
    
    # Filter unphysical heating values
    radius = getRadius(ra,dec)
    radiusCut = (radius < 18)

    idx  = (Fyoung > 0)*(Fyoung < 1)
    
    
    
#eMdMs16 = np.absolute(all_MdMs - np.log10(MdMs16))
#eMdMs84 = np.absolute(np.log10(MdMs84) - all_MdMs)
#mean_MdMs16 = nanmean(eMdMs16)
#mean_MdMs84 = nanmean(eMdMs84)

    # radial profiles
    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_rad_Fyoung(fig_a, ra[idx],dec[idx], Fyoung[idx])
        fig.savefig(path+"plot_rad_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_rad_TcISM(fig_a, ra[idx],dec[idx], Tc_ISM[idx])
        fig.savefig(path+"plot_rad_TcISM.png",format='png')

    if 0:
        idx = (FWarmYoung > 0)*(FWarmYoung < 1)

        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_rad_FWarmYoung(fig_a, ra[idx],dec[idx], FWarmYoung[idx])
        fig.savefig(path+"plot_rad_FWarmYoung.png",format='png')

    if 0:
        idx = (FColdYoung > 0)*(FColdYoung < 1)
        
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_rad_FColdYoung(fig_a, ra[idx],dec[idx], FColdYoung[idx])
        fig.savefig(path+"plot_rad_FColdYoung.png",format='png')

    # Parameters
    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_Lstar_Fyoung(fig_a, Lstar[idx], Fyoung[idx])
        fig.savefig(path+"plot_Lstar_Fyoung.pdf",format='pdf')
    
    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_sSFR_Fyoung(fig_a, sSFR[idx], Fyoung[idx])
        fig.savefig(path+"plot_sSFR_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_Mstars_Fyoung(fig_a, Mstars[idx], Fyoung[idx])
        fig.savefig(path+"plot_Mstars_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_Ldust_Fyoung(fig_a, Ldust[idx], Fyoung[idx])
        fig.savefig(path+"plot_Ldust_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_TwBC_Fyoung(fig_a, Tw_BC[idx], Fyoung[idx])
        fig.savefig(path+"plot_TwBC_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_TcISM_Fyoung(fig_a, Tc_ISM[idx], Fyoung[idx])
        fig.savefig(path+"plot_TcISM_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_LPAH_Fyoung(fig_a, LPAH[idx], Fyoung[idx])
        fig.savefig(path+"plot_LPAH_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_Mdust_Fyoung(fig_a, Mdust[idx], Fyoung[idx])
        fig.savefig(path+"plot_Mdust_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_SFR_Fyoung(fig_a, SFR[idx], Fyoung[idx])
        fig.savefig(path+"plot_SFR_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_MdMs_Fyoung(fig_a, MdMs[idx], Fyoung[idx])
        fig.savefig(path+"plot_MdMs_Fyoung.pdf",format='pdf')

    # Colours

    # Absolute heating
    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_FNUV_r_Lyoung(fig_a, FNUV[idx], Fr[idx], Lyoung[idx])
        fig.savefig(path+"plot_FNUV_r_Lyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_F160_F250_Lyoung(fig_a, F160[idx], F250[idx], Lyoung[idx])
        fig.savefig(path+"plot_F160_F250_Lyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_F250_F350_Lyoung(fig_a, F250[idx], F350[idx], Lyoung[idx])
        fig.savefig(path+"plot_F250_F350_Lyoung.png",format='png')


    # Heating fractions
    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_FNUV_r_Fyoung(fig_a, FNUV[idx], Fr[idx], Fyoung[idx])
        fig.savefig(path+"plot_FNUV_r_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F70_F100_Fyoung(fig_a, F70[idx], F100[idx], Fyoung[idx])
        fig.savefig(path+"plot_F70_F100_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F70_F250_Fyoung(fig_a, F70[idx], F250[idx], Fyoung[idx])
        fig.savefig(path+"plot_F70_F250_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F100_F250_Fyoung(fig_a, F100[idx], F250[idx], Fyoung[idx])
        fig.savefig(path+"plot_F100_F250_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F100_F500_Fyoung(fig_a, F100[idx], F500[idx], Fyoung[idx])
        fig.savefig(path+"plot_F100_F500_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_F160_F250_Fyoung(fig_a, F160[idx], F250[idx], Fyoung[idx])
        fig.savefig(path+"plot_F160_F250_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F160_F500_Fyoung(fig_a, F160[idx], F500[idx], Fyoung[idx])
        fig.savefig(path+"plot_F160_F500_Fyoung.pdf",format='pdf')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.12,0.11,0.87,0.87])
        plot_F250_F350_Fyoung(fig_a, F250[idx], F350[idx], Fyoung[idx])
        fig.savefig(path+"plot_F250_F350_Fyoung.png",format='png')

    if 0:
        fig  = plt.figure(figsize=(5,5))
        fig_a = plt.axes([0.11,0.11,0.88,0.88])
        plot_F250_F500_Fyoung(fig_a, F250[idx], F500[idx], Fyoung[idx])
        fig.savefig(path+"plot_F250_F500_Fyoung.pdf",format='pdf')


    # Important parameters vs Fyoung plot
    if 0:
    
        locplot = [[0.08,0.59,0.30,0.40],[0.38,0.59,0.30,0.40],[0.68,0.59,0.30,0.40],
                   [0.08,0.10,0.30,0.40],[0.38,0.10,0.30,0.40],[0.68,0.10,0.30,0.40]]
        
        fig  = plt.figure(figsize=(11,8))


        fig_a = plt.axes(locplot[0])
        plot_Mdust_Fyoung(fig_a,Mdust[idx], Fyoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_Mstars_Fyoung(fig_b,Mstars[idx], Fyoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_MdMs_Fyoung(fig_c,MdMs[idx], Fyoung[idx])
        fig_c.get_yaxis().set_visible(False)


        fig_d = plt.axes(locplot[3])
        plot_SFR_Fyoung(fig_d,SFR[idx], Fyoung[idx])

        fig_e = plt.axes(locplot[4])
        plot_sSFR_Fyoung(fig_e,sSFR[idx], Fyoung[idx])
        fig_e.get_yaxis().set_visible(False)

        fig_f = plt.axes(locplot[5])
        plot_Ldust_Fyoung(fig_f,Ldust[idx], Fyoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        fig.savefig("paperFigures/heatingParameters.png",format='png')

    # combined parameters vs Fyoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        
        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_Mdust_Fyoung(fig_a,Mdust[idx], Fyoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_Mstars_Fyoung(fig_b,Mstars[idx], Fyoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_MdMs_Fyoung(fig_c,MdMs[idx], Fyoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_SFR_Fyoung(fig_d,SFR[idx], Fyoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_sSFR_Fyoung(fig_e,sSFR[idx], Fyoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_TcISM_Fyoung(fig_f,Tc_ISM[idx], Fyoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_Ldust_Fyoung(fig_g,Ldust[idx], Fyoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_LPAH_Fyoung(fig_h,LPAH[idx], Fyoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_TwBC_Fyoung(fig_i,Tw_BC[idx], Fyoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        
        fig.savefig(path+"plot_HeatingParams_combo.png",format='png')

    # combined parameters vs FWarmyoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        idx = (FWarmYoung > 0)*(FWarmYoung < 1)

        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_Mdust_FWarmYoung(fig_a,Mdust[idx], FWarmYoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_Mstars_FWarmYoung(fig_b,Mstars[idx], FWarmYoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_MdMs_FWarmYoung(fig_c,MdMs[idx], FWarmYoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_SFR_FWarmYoung(fig_d,SFR[idx], FWarmYoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_sSFR_FWarmYoung(fig_e,sSFR[idx], FWarmYoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_TcISM_FWarmYoung(fig_f,Tc_ISM[idx], FWarmYoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_Ldust_FWarmYoung(fig_g,Ldust[idx], FWarmYoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_LPAH_FWarmYoung(fig_h,LPAH[idx], FWarmYoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_TwBC_FWarmYoung(fig_i,Tw_BC[idx], FWarmYoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        
        fig.savefig(path+"plot_WarmHeatingParams_combo.png",format='png')


    # combined parameters vs FColdyoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        
        idx = (FColdYoung > 0)*(FColdYoung < 1)
        
        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_Mdust_FColdYoung(fig_a,Mdust[idx], FColdYoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_Mstars_FColdYoung(fig_b,Mstars[idx], FColdYoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_MdMs_FColdYoung(fig_c,MdMs[idx], FColdYoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_SFR_FColdYoung(fig_d,SFR[idx], FColdYoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_sSFR_FColdYoung(fig_e,sSFR[idx], FColdYoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_TcISM_FColdYoung(fig_f,Tc_ISM[idx], FColdYoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_Ldust_FColdYoung(fig_g,Ldust[idx], FColdYoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_LPAH_FColdYoung(fig_h,LPAH[idx], FColdYoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_TwBC_FColdYoung(fig_i,Tw_BC[idx], FColdYoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        
        fig.savefig(path+"plot_ColdHeatingParams_combo.png",format='png')

    # important colours vs Fyoung plot
    if 1:
    
        locplot = [[0.08,0.15,0.30,0.83],[0.38,0.15,0.30,0.83],[0.68,0.15,0.30,0.83]]
        
        fig  = plt.figure(figsize=(11,4))
        
        fig_a = plt.axes(locplot[0])
        plot_FNUV_r_Fyoung(fig_a, FNUV[idx], Fr[idx], Fyoung[idx],radiusCut[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_F160_F250_Fyoung(fig_b, F160[idx], F250[idx], Fyoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_F250_F350_Fyoung(fig_c, F250[idx], F350[idx], Fyoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        fig.savefig("paperFigures/heatingColours.png",format='png', dpi=150)

    # important colours vs Lyoung plot
    if 1:
    
        locplot = [[0.08,0.15,0.30,0.83],[0.38,0.15,0.30,0.83],[0.68,0.15,0.30,0.83]]
        
        fig  = plt.figure(figsize=(11,4))
        
        fig_a = plt.axes(locplot[0])
        plot_FNUV_r_Lyoung(fig_a, FNUV[idx], Fr[idx], Lyoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_F160_F250_Lyoung(fig_b, F160[idx], F250[idx], Lyoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_F250_F350_Lyoung(fig_c, F250[idx], F350[idx], Lyoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        fig.savefig("paperFigures/absoluteHeatingColours.png",format='png', dpi=150)


    # combined colours vs Fyoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        
        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_FNUV_r_Fyoung(fig_a, FNUV[idx], Fr[idx], Fyoung[idx])
    
        fig_b = plt.axes(locplot[1])
        plot_F70_F100_Fyoung(fig_b, F70[idx], F100[idx], Fyoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_F70_F250_Fyoung(fig_c, F70[idx], F250[idx], Fyoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_F100_F250_Fyoung(fig_d, F100[idx], F250[idx], Fyoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_F100_F500_Fyoung(fig_e, F100[idx], F500[idx], Fyoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_F160_F250_Fyoung(fig_f, F160[idx], F250[idx], Fyoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_F160_F500_Fyoung(fig_g, F160[idx], F500[idx], Fyoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_F250_F350_Fyoung(fig_h, F250[idx], F350[idx], Fyoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_F250_F500_Fyoung(fig_i, F250[idx], F500[idx], Fyoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        fig.savefig(path+"plot_HeatingColours_combo.png",format='png')

    # combined colours vs FWarmYoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        idx = (FWarmYoung > 0)*(FWarmYoung < 1)
        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_FNUV_r_FWarmYoung(fig_a, FNUV[idx], Fr[idx], FWarmYoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_F70_F100_FWarmYoung(fig_b, F70[idx], F100[idx], FWarmYoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_F70_F250_FWarmYoung(fig_c, F70[idx], F250[idx], FWarmYoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_F100_F250_FWarmYoung(fig_d, F100[idx], F250[idx], FWarmYoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_F100_F500_FWarmYoung(fig_e, F100[idx], F500[idx], FWarmYoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_F160_F250_FWarmYoung(fig_f, F160[idx], F250[idx], FWarmYoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_F160_F500_FWarmYoung(fig_g, F160[idx], F500[idx], FWarmYoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_F250_F350_FWarmYoung(fig_h, F250[idx], F350[idx], FWarmYoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_F250_F500_FWarmYoung(fig_i, F250[idx], F500[idx], FWarmYoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        fig.savefig(path+"plot_WarmHeatingColours_combo.png",format='png')

    # combined colours vs FColdYoung plot
    if 0:
    
        locplot = [[0.08,0.72,0.30,0.27],[0.38,0.72,0.30,0.27],[0.68,0.72,0.30,0.27],
                   [0.08,0.39,0.30,0.27],[0.38,0.39,0.30,0.27],[0.68,0.39,0.30,0.27],
                   [0.08,0.06,0.30,0.27],[0.38,0.06,0.30,0.27],[0.68,0.06,0.30,0.27]]
        
        idx = (FColdYoung > 0)*(FColdYoung < 1)
        
        fig  = plt.figure(figsize=(11,12))
        
        
        fig_a = plt.axes(locplot[0])
        plot_FNUV_r_FColdYoung(fig_a, FNUV[idx], Fr[idx], FColdYoung[idx])
        
        fig_b = plt.axes(locplot[1])
        plot_F70_F100_FColdYoung(fig_b, F70[idx], F100[idx], FColdYoung[idx])
        fig_b.get_yaxis().set_visible(False)
        
        fig_c = plt.axes(locplot[2])
        plot_F70_F250_FColdYoung(fig_c, F70[idx], F250[idx], FColdYoung[idx])
        fig_c.get_yaxis().set_visible(False)
        
        
        fig_d = plt.axes(locplot[3])
        plot_F100_F250_FColdYoung(fig_d, F100[idx], F250[idx], FColdYoung[idx])
        
        fig_e = plt.axes(locplot[4])
        plot_F100_F500_FColdYoung(fig_e, F100[idx], F500[idx], FColdYoung[idx])
        fig_e.get_yaxis().set_visible(False)
        
        fig_f = plt.axes(locplot[5])
        plot_F160_F250_FColdYoung(fig_f, F160[idx], F250[idx], FColdYoung[idx])
        fig_f.get_yaxis().set_visible(False)
        
        
        fig_g = plt.axes(locplot[6])
        plot_F160_F500_FColdYoung(fig_g, F160[idx], F500[idx], FColdYoung[idx])
        
        fig_h = plt.axes(locplot[7])
        plot_F250_F350_FColdYoung(fig_h, F250[idx], F350[idx], FColdYoung[idx])
        fig_h.get_yaxis().set_visible(False)
        
        fig_i = plt.axes(locplot[8])
        plot_F250_F500_FColdYoung(fig_i, F250[idx], F500[idx], FColdYoung[idx])
        fig_i.get_yaxis().set_visible(False)
        
        fig.savefig(path+"plot_ColdHeatingColours_combo.png",format='png')


def getDensityKernel(x,y):
    
    # filter NaNs and infinite numbers
    idx = (np.isfinite(x) * np.isfinite(y))
    x = x[idx]
    y = y[idx]
    
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    return x[idx], y[idx], z[idx]

# RADIAL PROFILES
def plot_rad_Fyoung(fig, ra, dec, Fyoung):

    radius = getRadius(ra,dec)
    x, y, z = getDensityKernel(radius,100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('radius (kpc)',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(0,25)
    fig.set_ylim(-10,90)

def plot_rad_TcISM(fig, ra, dec, TcISM):
    
    radius = getRadius(ra,dec)
    x, y, z = getDensityKernel(radius,TcISM)
    
    fig.set_ylabel('$T_\mathrm{C}^\mathrm{ISM} / K$',fontsize=18)
    fig.set_xlabel('radius (kpc)',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(0,25)
    fig.set_ylim(7,33)

def plot_rad_FWarmYoung(fig, ra, dec, FWarmYoung):
    
    radius = getRadius(ra,dec)
    x, y, z = getDensityKernel(radius,100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('radius (kpc)',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(0,25)
    fig.set_ylim(-10,90)

def plot_rad_FColdYoung(fig, ra, dec, FColdYoung):
    
    radius = getRadius(ra,dec)
    x, y, z = getDensityKernel(radius,100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('radius (kpc)',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(0,25)
    fig.set_ylim(-10,90)


# PARAMETERS vs Fyoung
def plot_Lstar_Fyoung(fig, Lstar, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(Lstar),100*Fyoung)
    
    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\star^\mathrm{bol}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(3,9)
    fig.set_ylim(-10,90)

def plot_sSFR_Fyoung(fig, sSFR, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(sSFR),100*Fyoung)
    
    DL14_sSFR, DL14_Fyoung = np.loadtxt("/Users/saviaene/Documents/Research/M31/SKIRT/paperFigures/DeLooze2014sSFRHeating.dat", usecols=(0,1), unpack=True)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{sSFR}/\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    # Plot De Looze 2014 relation between F'young and sSFR
    fig.plot(DL14_sSFR,100.*10**DL14_Fyoung,'g+')
    xrange = np.linspace(-14,-8,100)
    fig.plot(xrange,100.*10**(0.415*xrange+4.045), 'k-')
    fig.set_xlim(-13.5,-8.5)
    fig.set_ylim(-10,90)
    # Plot first order polynomial fit
    #solution = FitPolynomial(fig,x,y,1)
    #solution = FitPolynomialLog(fig,x,y,1)
    
    #print solution
    #fig3.errorbar(6.8,-5.5, [[mean_MdMs16],[mean_MdMs84]],[[mean_Mskpc216],[mean_Mskpc284]], 'k.')

def plot_Mstars_Fyoung(fig, Mstars, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(Mstars),100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\star/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.9,8.9)
    fig.set_ylim(-10,90)

def plot_Ldust_Fyoung(fig, Ldust, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(Ldust),100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{dust}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.4,6.9)
    fig.set_ylim(-10,90)

def plot_TwBC_Fyoung(fig, TwBC, Fyoung):

    x, y, z = getDensityKernel(TwBC,100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{W}^\mathrm{BC} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(29,71)
    fig.set_ylim(-10,90)

def plot_TcISM_Fyoung(fig, TcISM, Fyoung):
    
    x, y, z = getDensityKernel(TcISM,100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{C}^\mathrm{ISM} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(10,31)
    fig.set_ylim(-10,90)

def plot_LPAH_Fyoung(fig, LPAH, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(LPAH),100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{PAH}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(3.6,6.1)
    fig.set_ylim(-10,90)

def plot_Mdust_Fyoung(fig, Mdust, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(Mdust),100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(1.9,4.6)
    fig.set_ylim(-10,90)

def plot_SFR_Fyoung(fig, SFR, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(SFR),100*Fyoung)

    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{SFR}/M_\odot\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-7.0,-3.0)
    fig.set_ylim(-10,90)

def plot_MdMs_Fyoung(fig, MdMs, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(MdMs),100*Fyoung)
    
    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\star)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-6.0,-1.0)
    fig.set_ylim(-10,90)


# PARAMETERS vs FWarmYoung
def plot_sSFR_FWarmYoung(fig, sSFR, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(sSFR),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{sSFR}/\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    # Plot De Looze 2014 relation between F'young and sSFR
    xrange = np.linspace(-14,-8,100)
    fig.plot(xrange,100.*10**(0.42*xrange+4.14), 'g--')
    fig.set_xlim(-13.5,-8.5)
    fig.set_ylim(-10,90)
    # Plot first order polynomial fit
    solution = FitPolynomial(fig,x,y,1)
    #solution = FitPolynomialLog(fig,x,y,1)
    
    print solution
    #fig3.errorbar(6.8,-5.5, [[mean_MdMs16],[mean_MdMs84]],[[mean_Mskpc216],[mean_Mskpc284]], 'k.')

def plot_Mstars_FWarmYoung(fig, Mstars, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(Mstars),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\star/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.9,8.9)
    fig.set_ylim(-10,90)

def plot_Ldust_FWarmYoung(fig, Ldust, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(Ldust),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{dust}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.4,6.9)
    fig.set_ylim(-10,90)

def plot_TwBC_FWarmYoung(fig, TwBC, FWarmYoung):
    
    x, y, z = getDensityKernel(TwBC,100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{W}^\mathrm{BC} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(29,71)
    fig.set_ylim(-10,90)

def plot_TcISM_FWarmYoung(fig, TcISM, FWarmYoung):
    
    x, y, z = getDensityKernel(TcISM,100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{C}^\mathrm{ISM} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(10,31)
    fig.set_ylim(-10,90)

def plot_LPAH_FWarmYoung(fig, LPAH, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(LPAH),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{PAH}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(3.6,6.1)
    fig.set_ylim(-10,90)

def plot_Mdust_FWarmYoung(fig, Mdust, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(Mdust),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(1.9,4.6)
    fig.set_ylim(-10,90)

def plot_SFR_FWarmYoung(fig, SFR, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(SFR),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{SFR}/M_\odot\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-7.0,-3.0)
    fig.set_ylim(-10,90)

def plot_MdMs_FWarmYoung(fig, MdMs, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(MdMs),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\star)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-6.0,-1.0)
    fig.set_ylim(-10,90)


# PARAMETERS vs FColdYoung
def plot_sSFR_FColdYoung(fig, sSFR, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(sSFR),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{sSFR}/\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    # Plot De Looze 2014 relation between F'young and sSFR
    xrange = np.linspace(-14,-8,100)
    fig.plot(xrange,100.*10**(0.42*xrange+4.14), 'g--')
    fig.set_xlim(-13.5,-8.5)
    fig.set_ylim(-10,90)
    # Plot first order polynomial fit
    solution = FitPolynomial(fig,x,y,1)
    #solution = FitPolynomialLog(fig,x,y,1)
    
    print solution
#fig3.errorbar(6.8,-5.5, [[mean_MdMs16],[mean_MdMs84]],[[mean_Mskpc216],[mean_Mskpc284]], 'k.')

def plot_Mstars_FColdYoung(fig, Mstars, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(Mstars),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\star/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.9,8.9)
    fig.set_ylim(-10,90)

def plot_Ldust_FColdYoung(fig, Ldust, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(Ldust),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{dust}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(4.4,6.9)
    fig.set_ylim(-10,90)

def plot_TwBC_FColdYoung(fig, TwBC, FColdYoung):
    
    x, y, z = getDensityKernel(TwBC,100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{W}^\mathrm{BC} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(29,71)
    fig.set_ylim(-10,90)

def plot_TcISM_FColdYoung(fig, TcISM, FColdYoung):
    
    x, y, z = getDensityKernel(TcISM,100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$T_\mathrm{C}^\mathrm{ISM} / K$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(10,31)
    fig.set_ylim(-10,90)

def plot_LPAH_FColdYoung(fig, LPAH, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(LPAH),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(L_\mathrm{PAH}/L_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(3.6,6.1)
    fig.set_ylim(-10,90)

def plot_Mdust_FColdYoung(fig, Mdust, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(Mdust),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\odot)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(1.9,4.6)
    fig.set_ylim(-10,90)

def plot_SFR_FColdYoung(fig, SFR, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(SFR),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(\mathrm{SFR}/M_\odot\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-7.0,-3.0)
    fig.set_ylim(-10,90)

def plot_MdMs_FColdYoung(fig, MdMs, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(MdMs),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(M_\mathrm{dust}/M_\star)$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-6.0,-1.0)
    fig.set_ylim(-10,90)



# COLOURS vs Fyoung
def plot_FNUV_r_Fyoung(fig, nuv, r, Fyoung, radiusCut):
    
    x_r = np.log10(nuv[~radiusCut]/r[~radiusCut])
    y_r = 100*Fyoung[~radiusCut]

    x_r2 = np.log10(nuv[radiusCut]/r[radiusCut])
    y_r2 = 100*Fyoung[radiusCut]

    x, y, z = getDensityKernel(np.log10(nuv[radiusCut]/r[radiusCut]),100*Fyoung[radiusCut])
    
    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$NUV-r$',fontsize=18)
    fig.scatter(x_r,y_r, c='black',s=1 , alpha=0.1)
    #fig.scatter(x_r2,y_r2, c='red',s=1 , alpha=0.1)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-2.7,-0.51)
    fig.set_ylim(-10,110)
    #solution = FitPolynomial(fig,x,y,1)
    #print solution

def plot_F70_F100_Fyoung(fig, F70, F100, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F100),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{100})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,4)
    fig.set_ylim(-10,110)

def plot_F70_F250_Fyoung(fig, F70, F250, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F250),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,1)
    fig.set_ylim(-10,110)

def plot_F100_F250_Fyoung(fig, F100, F250, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F250),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,2)
    fig.set_ylim(-10,110)

def plot_F100_F500_Fyoung(fig, F100, F500, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F500),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-3,2)
    fig.set_ylim(-10,110)

def plot_F160_F250_Fyoung(fig, F160, F250, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F250),100*Fyoung)
    
    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.6,2.5)
    fig.set_ylim(-10,110)
    fig.set_xticks([-1.5,-1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5])

def plot_F160_F500_Fyoung(fig, F160, F500, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F500),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.5)
    fig.set_ylim(-10,110)

def plot_F250_F350_Fyoung(fig, F250, F350, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F350),100*Fyoung)
    
    fig.set_ylabel('$F_\mathrm{unev.} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{350})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.4,1.8)
    fig.set_ylim(-10,110)
    fig.set_xticks([-1.0, -0.5, 0., 0.5, 1.0, 1.5])


def plot_F250_F500_Fyoung(fig, F250, F500, Fyoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F500),100*Fyoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.3)
    fig.set_ylim(-10,110)



# COLOURS vs Lyoung
def plot_FNUV_r_Lyoung(fig, nuv, r, Lyoung):
    
    x, y, z = getDensityKernel(np.log10(nuv/r),np.log10(Lyoung))

    fig.set_ylabel('$\log(L_\mathrm{unev.}/L_\odot) $',fontsize=18)
    fig.set_xlabel('$NUV-r$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-2.7,-0.51)
    fig.set_ylim(0,7)
    #solution = FitPolynomial(fig,x,y,1)
    #print solution

def plot_F160_F250_Lyoung(fig, F160, F250, Lyoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F250),np.log10(Lyoung))
    
    fig.set_ylabel('$\log(L_\mathrm{unev.}/L_\odot) $',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.6,2.5)
    fig.set_ylim(0,7)
    fig.set_xticks([-1.5,-1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5])

def plot_F250_F350_Lyoung(fig, F250, F350, Lyoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F350),np.log10(Lyoung))
    
    fig.set_ylabel('$\log(L_\mathrm{unev.}/L_\odot) $',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{350})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.4,1.8)
    fig.set_ylim(0,7)
    fig.set_xticks([-1.0, -0.5, 0., 0.5, 1.0, 1.5])




# COLOURS vs FWarmYoung
def plot_FNUV_r_FWarmYoung(fig, nuv, r, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(nuv/r),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$NUV-r$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-2.7,-0.5)
    fig.set_ylim(-10,110)
    solution = FitPolynomial(fig,x,y,1)
    print solution

def plot_F70_F100_FWarmYoung(fig, F70, F100, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F100),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{100})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,4)
    fig.set_ylim(-10,110)

def plot_F70_F250_FWarmYoung(fig, F70, F250, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F250),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,1)
    fig.set_ylim(-10,110)

def plot_F100_F250_FWarmYoung(fig, F100, F250, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F250),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,2)
    fig.set_ylim(-10,110)

def plot_F100_F500_FWarmYoung(fig, F100, F500, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F500),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-3,2)
    fig.set_ylim(-10,110)

def plot_F160_F250_FWarmYoung(fig, F160, F250, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F250),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.6,1.4)
    fig.set_ylim(-10,110)

def plot_F160_F500_FWarmYoung(fig, F160, F500, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F500),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.5)
    fig.set_ylim(-10,110)

def plot_F250_F350_FWarmYoung(fig, F250, F350, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F350),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{350})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-0.5,0.7)
    fig.set_ylim(-10,110)

def plot_F250_F500_FWarmYoung(fig, F250, F500, FWarmYoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F500),100*FWarmYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{w} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.3)
    fig.set_ylim(-10,110)



# COLOURS vs FColdYoung
def plot_FNUV_r_FColdYoung(fig, nuv, r, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(nuv/r),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$NUV-r$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-2.7,-0.5)
    fig.set_ylim(-10,110)
    solution = FitPolynomial(fig,x,y,1)
    print solution

def plot_F70_F100_FColdYoung(fig, F70, F100, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F100),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{100})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,4)
    fig.set_ylim(-10,110)

def plot_F70_F250_FColdYoung(fig, F70, F250, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F70/F250),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{70}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,1)
    fig.set_ylim(-10,110)

def plot_F100_F250_FColdYoung(fig, F100, F250, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F250),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-4,2)
    fig.set_ylim(-10,110)

def plot_F100_F500_FColdYoung(fig, F100, F500, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F100/F500),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{100}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-3,2)
    fig.set_ylim(-10,110)

def plot_F160_F250_FColdYoung(fig, F160, F250, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F250),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{250})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1.6,1.4)
    fig.set_ylim(-10,110)

def plot_F160_F500_FColdYoung(fig, F160, F500, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F160/F500),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{160}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.5)
    fig.set_ylim(-10,110)

def plot_F250_F350_FColdYoung(fig, F250, F350, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F350),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{350})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-0.5,0.7)
    fig.set_ylim(-10,110)

def plot_F250_F500_FColdYoung(fig, F250, F500, FColdYoung):
    
    x, y, z = getDensityKernel(np.log10(F250/F500),100*FColdYoung)
    
    fig.set_ylabel('$F^\prime_\mathrm{young}^\mathrm{c} [\%]$',fontsize=18)
    fig.set_xlabel('$\log(F_{250}/F_{500})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-1,2.3)
    fig.set_ylim(-10,110)


#OTHER
def plot_fUV_sSFR(fig, fUV, sSFR, flag):

    fig.set_ylabel('UV heating [$\%$]',fontsize=18)
    fig.set_xlabel('$\log \mathrm{sSFR}/\mathrm{yr}^{-1})$',fontsize=18)
    fig.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    fig.set_xlim(-13.5,-8.5)
    fig.set_ylim(0,100)
    
    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)

    logx = np.log10(sSFR[flag==0])
    logy = np.log10(100.*fUV[flag==0])

    p1, succes, infodict,mesg,ier = optimize.leastsq(PowerLawErrorFunction, [1,1], args=(logx,logy), full_output=1)
    rms = np.sqrt(((infodict['fvec']**2).sum())/(len(logx)-1))
    print "power law best fit a*x^b, rms"
    print p1, rms
    
  
    xrange = fig.get_xlim()
    plotx = np.linspace(xrange[0],xrange[1],100)
    ploty = 10.**PowerLaw(p1,plotx)

    fig.plot(plotx,ploty, 'k-',linewidth=2)
    #fig.plot(plotx,72537.3507523*np.power(10**plotx,0.30174404), 'r-',linewidth=2)


def FitPolynomial(fig, x, y, n):

#solution, res = np.polynomial.polynomial.polyfit(x,y,order,full=True)
    
    solution, C_p = np.polyfit(x, y, n, cov=True)  # C_z is estimated covariance matrix

    # Do the interpolation for plotting:
    xrange = fig.get_xlim()
    t = np.linspace(xrange[0],xrange[1],100)
    # Matrix with rows 1, t, t**2, ...:
    TT = np.vstack([t**(n-i) for i in range(n+1)]).T
    yi = np.dot(TT, solution)  # matrix multiplication calculates the polynomial values
    C_yi = np.dot(TT, np.dot(C_p, TT.T)) # C_y = TT*C_z*TT.T
    sig_yi = np.sqrt(np.diag(C_yi))  # Standard deviations are sqrt of diagonal


    fig.plot(t,yi, 'k-')
    # fig.plot(t,yi+sig_yi, 'k--')
    # fig.plot(t,yi-sig_yi, 'k--')

    return solution

# Take the log of the y value
def FitPolynomialLog(fig, x, y, order):
    
    y = np.log10(y)

    solution, res, other, output, stuff = np.polynomial.polynomial.polyfit(x,y,order,full=True)
    print res
    xrange = fig.get_xlim()
    plotRangex = np.linspace(xrange[0],xrange[1],100)
    
    fig.plot(plotRangex,10**polynomial(plotRangex,solution), 'k-')
    
    return solution

def PowerLaw(params,logx):
    return params[0] + params[1] * logx

def PowerLawErrorFunction(params,logx,logy):
    return PowerLaw(params,logx) - logy


def polynomial(x, solution):
    
    x = np.array(x)
    solution = np.array(solution)
    
    y = np.zeros(len(x))
    for i in range(0,len(solution)):
        y += solution[i] * x**i
    
    return y

def getRadius(ra,dec):
    
    # M31 properties
    PA = 38.1
    inclin = 77.5
    centre = [10.77615,41.353394] # Center of our current FOV.
    #centre = [10.612508,41.208711] # Center from Viaene et al 2014
    dist = 0.785
    
    #Deproject the pixels to a physical radial distance.
    # convert angles to radians
    PA = ((90.-PA) / 180.0 * math.pi)
    inclin = inclin / 180.0 * math.pi
    
    radius = np.zeros(len(ra))
    for i in range(0,len(ra)):
        Xsquare = ((ra[i] - centre[0])*math.cos(dec[i] / 180.0 * math.pi)*math.cos(PA) + (dec[i] - centre[1])*math.sin(PA))**2
        Ysquare = (-(ra[i] - centre[0])*math.cos(dec[i] / 180.0 * math.pi)*math.sin(PA) + (dec[i] - centre[1])*math.cos(PA))**2
        radius[i] = math.sqrt(Xsquare + Ysquare / math.cos(inclin)**2.0)
        radius[i] = 2.0 * dist * 1000.0 * math.tan(radius[i]*math.pi/(180.0*2.0))
    
    return radius



if __name__ == '__main__':
    main()