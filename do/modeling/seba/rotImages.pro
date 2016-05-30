pro rotImages
path = 'FinalRun/maps/'

im=readfits(path+"im500Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 313.67516
sxaddpar, rhdr, 'CDELT1', 0.01 # graden
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 342.23715
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim500Jy.fits',rim,rhdr

im=readfits(path+"im350Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 312.47332
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 341.65097
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim350Jy.fits',rim,rhdr

im=readfits(path+"im250Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 313.28115
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 342.04498
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim250Jy.fits',rim,rhdr

im=readfits(path+"im100Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim100Jy.fits',rim,rhdr

im=readfits(path+"im160Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 308.16025
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 351.54438
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim160Jy.fits',rim,rhdr

im=readfits(path+"im70Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 144.6502
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 179.34469
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim70Jy.fits',rim,rhdr

im=readfits(path+"im24Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 147.46725
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 159.91493
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim24Jy.fits',rim,rhdr

im=readfits(path+"imW4Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.70513
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 346.2763
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimW4Jy.fits',rim,rhdr

im=readfits(path+"imW3Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.70513
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 346.2763
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimW3Jy.fits',rim,rhdr

im=readfits(path+"imW2Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.70513
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 346.2763
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimW2Jy.fits',rim,rhdr

im=readfits(path+"imW1Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.70513
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 346.2763
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimW1Jy.fits',rim,rhdr

im=readfits(path+"im8Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 146.95841
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 157.95821
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim8Jy.fits',rim,rhdr

im=readfits(path+"im5.8Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 146.95841
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 157.95821
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim5.8Jy.fits',rim,rhdr

im=readfits(path+"im4.5Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 146.95841
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 157.95821
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim4.5Jy.fits',rim,rhdr

im=readfits(path+"im3.6Jy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 146.95841
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 157.95821
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotim3.6Jy.fits',rim,rhdr

im=readfits(path+"imzJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimzJy.fits',rim,rhdr

im=readfits(path+"imiJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimiJy.fits',rim,rhdr

im=readfits(path+"imrJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimrJy.fits',rim,rhdr

im=readfits(path+"imgJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimgJy.fits',rim,rhdr

im=readfits(path+"imuJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 311.31112
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 347.08413
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimuJy.fits',rim,rhdr

im=readfits(path+"imNUVJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.73724
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 255.41167
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimNUVJy.fits',rim,rhdr

im=readfits(path+"imFUVJy.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.73724
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 255.41167
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotimFUVJy.fits',rim,rhdr

im=readfits(path+"M31_Ldust_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_Ldust_filtered.fits',rim,rhdr

im=readfits(path+"M31_Mdust_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_Mdust_filtered.fits',rim,rhdr

im=readfits(path+"M31_T_C_ISM_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_T_C_ISM_filtered.fits',rim,rhdr

im=readfits(path+"M31_Mstar_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_Mstar_filtered.fits',rim,rhdr

im=readfits(path+"M31_sSFR_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_sSFR_filtered.fits',rim,rhdr

im=readfits(path+"M31_SFR_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_SFR_filtered.fits',rim,rhdr

im=readfits(path+"M31_tauv_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_tauv_filtered.fits',rim,rhdr

im=readfits(path+"M31_tvism_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_tvism_filtered.fits',rim,rhdr

im=readfits(path+"M31_T_W_BC_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_T_W_BC_filtered.fits',rim,rhdr

im=readfits(path+"M31_xi_PAH_tot_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_xi_PAH_tot_filtered.fits',rim,rhdr

im=readfits(path+"M31_fmu_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_fmu_filtered.fits',rim,rhdr

im=readfits(path+"M31_xi_C_tot_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_xi_C_tot_filtered.fits',rim,rhdr

im=readfits(path+"M31_L_C_tot_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_L_C_tot_filtered.fits',rim,rhdr

im=readfits(path+"M31_L_PAH_tot_filtered.fits",hdr)
hrot, im, hdr, rim, rhdr, 308, -1, -1, 0
sxaddpar, rhdr, 'CRPIX1', 254.54507
sxaddpar, rhdr, 'CDELT1', 0.01
sxaddpar, rhdr, 'CRVAL1', 0.0
sxaddpar, rhdr, 'CTYPE1', 'LINEAR'
sxaddpar, rhdr, 'CRPIX2', 254.80568
sxaddpar, rhdr, 'CDELT2', 0.01
sxaddpar, rhdr, 'CRVAL2', 0.0
sxaddpar, rhdr, 'CTYPE2', 'LINEAR'
sxaddpar, rhdr, 'CROTA1', 0.0
sxaddpar, rhdr, 'CROTA2', 0.0
writefits, path+'plotM31_L_PAH_tot_filtered.fits',rim,rhdr

end
