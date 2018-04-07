import healpy as h

nside=512
boost=700.

mask=h.read_map("../datain/maps/mask/mask_gal60_nside512.fits")
mask=h.ud_grade(mask,nside)

bl=h.read_cl("../datain/spectra/353_TEB_beam.fits")
cldust=h.read_cl("../datain/spectra/dust_power_spectrum.fits")

n1=h.read_map("../datain/maps/simobsTQU/cmb_noise_TQU_map_year1_nside512_nrlz1.fits",(0,1,2))
n2=h.read_map("../datain/maps/simobsTQU/cmb_noise_TQU_map_year2_nside512_nrlz1.fits",(0,1,2))


map=zeros(h.nside2npix(nside),float)
alm=h.map2alm(map,lmax=2*nside)

l=3 ; m = 3
alm[h.Alm.getidx(2*nside,l,m)]=3.

l=6 ; m =1
alm[h.Alm.getidx(2*nside,l,m)]=2.

l=9 ; m =3
alm[h.Alm.getidx(2*nside,l,m)]=0.5

l=12 ; m =3
alm[h.Alm.getidx(2*nside,l,m)]=0.2

#l=17 ; m =5
#alm[h.Alm.getidx(2*nside,l,m)]=1.
#
#l=17 ; m =3
#alm[h.Alm.getidx(2*nside,l,m)]=2.
#

mockfrg=h.alm2map(alm,nside)
mockfrg=mockfrg/max(mockfrg)

tqufrg=[mockfrg[:]*0,mockfrg,mockfrg[:]*0.]
h.write_map("../datain/maps/simobsTQU/mockfrg_TQU.fits",tqufrg)

alm=h.map2alm(tqufrg,iter=8,lmax=2*nside)
for l in range(2*nside+1):
        for m in range(l+1):
                index=h.Alm.getidx(2*nside,l,m)
                for i in range(3):
                        alm[i][index]=alm[i][index]*bl[i][l]  #*sqrt(cldust[i][l])

frg=h.alm2map(alm,nside)
h.write_map("../datain/maps/simobsTQU/mockfrg_TQU_smoothed.fits",frg)

frg1=(frg[0]*boost+n1[0], frg[1]*boost+n1[1], frg[2]*boost+n1[2])
frg2=(frg[0]*boost+n2[0], frg[1]*boost+n2[1], frg[2]*boost+n2[2])

h.write_map("../datain/maps/simobsTQU/mockfrg_yr1_TQU.fits",frg1)
h.write_map("../datain/maps/simobsTQU/mockfrg_yr2_TQU.fits",frg2)

