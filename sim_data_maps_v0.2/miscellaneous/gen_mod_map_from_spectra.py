import healpy as h

nrlz=1
nside=512
pathout="../datain/maps/simobsTQU/"
cl=h.read_cl("../datain/spectra/dust_power_spectrum.fits")
bl=h.read_cl("../datain/spectra/353_TEB_beam.fits")

cl[1]=cl[1]*(bl[1][:1025]**2.)
cl[2]=cl[2]*(bl[2][:1025]**2.)

dipmod=zeros(h.nside2npix(nside),float)
quadmod=zeros(h.nside2npix(nside),float)

for i in range(h.nside2npix(nside)):
	theta=h.pix2ang(nside,i)[0]
	dipmod[i]=1+cos(theta)
	quadmod[i]=1+cos(theta)**2.

for i in range(nrlz):
	print i
	ofrg=h.synfast(cl,nside,pol=True,new=True,verbose=False)
	filename=pathout + "dustmap_nside" + str(nside) + "_nrlz" + str(i+1) + ".fits"
	h.write_map(filename,ofrg)

	alm=h.map2alm(ofrg,lmax=2*nside)
	e=h.alm2map(alm[1],nside) ; e=e*dipmod ; alme=h.map2alm(e,lmax=2*nside)
	b=h.alm2map(alm[2],nside) ; b=b*dipmod ; almb=h.map2alm(b,lmax=2*nside)
	alm=(alm[0],alme,almb)
	frg=h.alm2map(alm,nside=nside)	
	filename=pathout + "dipmod_dustmap_nside" + str(nside) + "_nrlz" + str(i+1) + ".fits"
	h.write_map(filename,frg)

n1=h.read_map(pathout + "cmb_noise_TQU_map_year1_nside512_nrlz1.fits",(0,1,2))
n2=h.read_map(pathout + "cmb_noise_TQU_map_year2_nside512_nrlz1.fits",(0,1,2))

frg1=(frg[0]+n1[0], frg[1]+n1[1], frg[2]+n1[2])
frg2=(frg[0]+n2[0], frg[1]+n2[1], frg[2]+n2[2])

h.write_map(pathout + "dipmod_dustmap_yr1_TQU.fits",frg1)
h.write_map(pathout + "dipmod_dustmap_yr2_TQU.fits",frg2)

h.write_map(pathout + "dipmod.fits",dipmod)
