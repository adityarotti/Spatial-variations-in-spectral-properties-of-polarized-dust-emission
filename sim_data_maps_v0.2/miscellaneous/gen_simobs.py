import healpy as h

nrlz=1
nside=512

pathout="../datain/maps/simobsTQU/"

cl=h.read_cl("../datain/spectra/dust_power_spectrum.fits")
cl[1]=cl[1]*100
cl[2]=cl[2]*100

bl=h.read_cl("../datain/spectra/353_TEB_beam.fits")

n1=h.read_map("../datain/maps/simobsTQU/cmb_noise_TQU_map_year1_nside512_nrlz1.fits",(0,1,2))
n2=h.read_map("../datain/maps/simobsTQU/cmb_noise_TQU_map_year2_nside512_nrlz1.fits",(0,1,2))

m=h.synfast(cl,nside,pixwin=True,pol=True,new=True,verbose=False)
filename="../datain/maps/simobsTQU/frg.fits"
h.write_map(filename,m)

# Smoothing alms
alm=h.map2alm(m,lmax=2*nside)
for l in range(2*nside+1):
	for m in range(l+1):
		index=h.Alm.getidx(2*nside,l,m)
		for i in range(3):
			alm[i][index]=alm[i][index]*bl[i][l]
		
m=h.alm2map(alm,nside)
filename="../datain/maps/simobsTQU/frg_smoothed.fits"
h.write_map(filename,m)

m1=(m[0]+n1[0],m[1]+n1[1],m[2]+n1[2])
filename="../datain/maps/simobsTQU/TQU_map_year1_nside512_nrlz1.fits"
h.write_map(filename,m1)

m2=(m[0]+n2[0],m[1]+n2[1],m[2]+n2[2])
filename="../datain/maps/simobsTQU/TQU_map_year2_nside512_nrlz1.fits"
h.write_map(filename,m2)

