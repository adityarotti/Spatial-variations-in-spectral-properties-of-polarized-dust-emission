import healpy as h

nrlz=50
nside=512
pathout="../datain/maps/dustmapTQU/"
cl=h.read_cl("../datain/spectra/dust_power_spectrum.fits")

for i in range(nrlz):
	print i
	m=h.synfast(cl,nside,pol=True,new=True,verbose=False)
	filename=pathout + "dustmap_nside" + str(nside) + "_nrlz" + str(i+1) + ".fits"
	h.write_map(filename,m)
