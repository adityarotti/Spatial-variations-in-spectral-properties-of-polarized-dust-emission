import healpy as h

nside=2048
lmax=2*nside
fwhm=(10./60.)*pi/180

blt=h.gauss_beam(fwhm,lmax=lmax)
ble=h.gauss_beam(fwhm*5,lmax=lmax)
blb=h.gauss_beam(fwhm*2,lmax=lmax)
bl=[blt,ble,blb]
h.write_cl("bl.fits",bl)

