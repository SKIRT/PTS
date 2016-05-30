#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.make_pdfs

# -----------------------------------------------------------------

pro getHeatedPixels

; M31Full.reg
ac = 10.664615
dc = 41.241692
;major axis in arcsec!
maxis = 5781.739
ratio = maxis / 1600.868
pangle = 38.
;pixel scale in arcsec!
pixelscale = 36.
outfile = "modelChecks/iteration5_J14/pixelHeating.dat"

; Heating maps
Fyoung = readfits("modelChecks/iteration5_J14/heatingTotYoung.fits",hdrYoung)
Fold   = readfits("modelChecks/iteration5_J14/heatingTotOld.fits",hdrOld)
Lold   = readfits("modelChecks/iteration5_J14/Ldust_old.fits",hdrLOld)
Lyoung = readfits("modelChecks/iteration5_J14/Ldust_young.fits",hdrLYoung)
Ltot   = readfits("modelChecks/iteration5_J14/Ldust_tot.fits",hdrLTot)

; Bolometric stellar luminosity map
Lstar   = readfits("modelChecks/iteration5_J14/nodust/M31_212full_Lstar.fits",hdrLstar)


;FWarmYoung = readfits("modelChecks/iteration4_J14/mod270/heatingTotWarmYoung.fits",hdrWarmYoung)
;FWarmOld   = readfits("modelChecks/iteration4_J14/mod270/heatingTotWarmOld.fits",hdrWarmOld)
;FColdYoung = readfits("modelChecks/iteration4_J14/mod270/heatingTotColdYoung.fits",hdrColdYoung)
;FColdOld   = readfits("modelChecks/iteration4_J14/mod270/heatingTotColdOld.fits",hdrColdOld)

; reference image
F36    = readfits("Images36/im3.6Jy.fits",hdr36)

xaxis36 = sxpar(hdr36, "NAXIS1")
yaxis36 = sxpar(hdr36, "NAXIS2")
mapFyoung = fltarr(xaxis36,yaxis36)*sqrt(-1)
mapFold   = fltarr(xaxis36,yaxis36)*sqrt(-1)
mapLold   = fltarr(xaxis36,yaxis36)*sqrt(-1)
mapLyoung = fltarr(xaxis36,yaxis36)*sqrt(-1)
mapLtot   = fltarr(xaxis36,yaxis36)*sqrt(-1)
mapLstar  = fltarr(xaxis36,yaxis36)*sqrt(-1)

;mapFWarmYoung = fltarr(xaxis36,yaxis36)*sqrt(-1)
;mapFWarmOld   = fltarr(xaxis36,yaxis36)*sqrt(-1)
;mapFColdYoung = fltarr(xaxis36,yaxis36)*sqrt(-1)
;mapFColdOld   = fltarr(xaxis36,yaxis36)*sqrt(-1)

xaxis = sxpar(hdrOld, "NAXIS1")
yaxis = sxpar(hdrOld, "NAXIS2")
dx = 46
dy = 31

for l=0L,xaxis -1 ,1 do begin
    for k=0L,yaxis -1 ,1 do begin

        mapFyoung[dx+l,dy+k] = Fyoung[l,k]
        mapFold[dx+l,dy+k]   = Fold[l,k]

        mapLold[dx+l,dy+k]   = Lold[l,k]
        mapLyoung[dx+l,dy+k] = Lyoung[l,k]
        mapLtot[dx+l,dy+k]   = Ltot[l,k]

        mapLstar[dx+l,dy+k]   = Lstar[l,k]

        ;mapFWarmYoung[dx+l,dy+k] = FWarmYoung[l,k]
        ;mapFWarmOld[dx+l,dy+k]   = FWarmOld[l,k]

        ;mapFColdYoung[dx+l,dy+k] = FColdYoung[l,k]
        ;mapFColdOld[dx+l,dy+k]   = FColdOld[l,k]

    endfor
endfor

;writefits,"test.fits",mapFyoung, hdr36

EXTAST, hdr36, astr
ad2xy, ac,dc,astr,xc,yc
DIST_ELLIPSE, ell, [xaxis36,yaxis36], xc, yc, ratio, pangle
region =  where( (pixelscale*ell) lt maxis )

size = SIZE(region, /N_ELEMENTS)
print, size

openw,lun, outfile,/get_lun, WIDTH=1024
printf, lun, "# Pixel    Fold    Fyoung    Lold    Lyoung    Ltot    Lstar"
for l=0L,size -1 ,1 do begin
   printf,lun, float(l), mapFold[region[l]], mapFyoung[region[l]], mapLold[region[l]], mapLyoung[region[l]], mapLtot[region[l]], mapLstar[region[l]];, mapFWarmOld[region[l]], mapFWarmYoung[region[l]], mapFColdOld[region[l]], mapFColdYoung[region[l]]
endfor

Free_lun, lun

end
 
