import healpy as h
import numpy as np
import matplotlib.pyplot as plt


def return_ea_theta(thetamin,thetamax,noea):
	thetamin=thetamin*np.pi/180.
	thetamax=thetamax*np.pi/180.
	A0=2.*np.pi*(np.cos(thetamin)-np.cos(thetamax))
	A=A0/noea

	thetabin=[]
	thetabin=np.append(thetabin,thetamin)
	for i in range(noea-1):
	    x=np.arccos(np.cos(thetabin[i]) - (A/(2.*np.pi)))
	    thetabin=np.append(thetabin,x)
	thetabin=np.append(thetabin,thetamax)

	thetaN=thetabin ; thetaS=(np.pi-thetaN)[::-1]
	thetaL=np.append(thetaN[:-1],thetaS[:-1]) ; thetaH=np.append(thetaN[1:],thetaS[1:])
	theta0=(thetaL+thetaH)/2
	return thetaL,thetaH,theta0
	
	

def return_ring_mask(nside,thetaL,thetaH,apow,mp=False,pathout="./"):
	thetaL=thetaL*np.pi/180.
	thetaH=thetaH*np.pi/180.
	apow=apow*np.pi/180.
	mask=np.zeros(h.nside2npix(nside),float)
	for i in range(h.nside2npix(nside)):
		theta=h.pix2ang(nside,i)[0]
		if thetaL>0 and thetaH<np.pi:
           		if theta>=thetaL+apow and theta<thetaH-apow:
                		mask[i]=1.
            		elif theta>=thetaL and theta<thetaL+apow:
                		mask[i]=np.cos(((theta-(thetaL+apow))/apow)*np.pi/2.)**2.
            		elif theta>thetaH-apow and theta<=thetaH:
                		mask[i]=np.cos(((theta-(thetaH-apow))/apow)*np.pi/2.)**2.
        	elif thetaL==0.:
            		if theta<=thetaH-apow:
                		mask[i]=1.
            		elif theta>thetaH-apow and theta<=thetaH:
                		mask[i]=np.cos(((theta-(thetaH-apow))/apow)*np.pi/2.)**2.
		elif thetaH==np.pi:
            		if theta>=thetaL+apow:
                		mask[i]=1.
            		elif theta>thetaL and theta<=thetaL+apow:
                		mask[i]=np.cos(((theta-(thetaL+apow))/apow)*np.pi/2.)**2.

	if mp:
		h.mollview(mask,title="")
		filename = pathout + "ring_mask_thetamin" + str(round(thetaL*180/np.pi,1)) + "_thetamax" + str(round(thetaH*180/np.pi,1)) 
		filename = filename + "_apow" + str(round(apow*180./np.pi,1)) +  "_nside" + str(nside) + ".pdf"
		plt.savefig(filename,dpi=150,bbox_inches='tight')
		plt.close()	
				
	return mask

def return_band_mask(wband,apow,nside,mp=False,pathout="./"):
	"""
	Take the band width(wband) and apodization width(apow) in degrees and the 
	NSIDE of the map as input and returns a apodized band mask
	"""
	mask=np.zeros(h.nside2npix(nside),float)
	for i in range(h.nside2npix(nside)):
		theta=h.pix2ang(nside,i)[0]*180/np.pi
	    	if theta<=90-wband-apow or theta>=90+wband+apow:
	        	mask[i]=1.
		elif theta>90-wband-apow and theta<=90-wband:
			mask[i]=np.cos(((theta-90+wband+apow)/apow)*np.pi/2.)
		elif theta<=90+wband+apow and theta>90+wband:
			mask[i]=np.cos(((theta-90-wband-apow)/apow)*np.pi/2.)

	if mp:
		h.mollview(mask,title="")
		filename=pathout + "band_mask_wband" + str(wband) + "_apow" + str(apow) + "_nside" + str(nside) + ".pdf"
		plt.savefig(filename,dpi=150,bbox_inches='tight')
		plt.close()
	return mask
		
	
	


def return_disc_mask(discsize,apow,nside,nsideout,pixel_index,mp=False,pathout="./"):
	mask=np.zeros(h.nside2npix(nside),float)
	theta0=h.pix2ang(nsideout,pixel_index)[0] ; phi0=h.pix2ang(nsideout,pixel_index)[1]
	mask_center=h.ang2pix(nside,theta0,phi0)
	# The coordinate of the pixel center updated to that in the mask resolution
	theta0=h.pix2ang(nside,mask_center)[0] ; phi0=h.pix2ang(nside,mask_center)[1]

	pixnum=np.arange(h.nside2npix(nside))
	theta=h.pix2ang(nside,pixnum)[0]
    	phi=h.pix2ang(nside,pixnum)[1]
	dotprod=np.sin(theta)*np.sin(theta0)*np.cos(phi-phi0)+np.cos(theta)*np.cos(theta0)
	dtheta=np.arccos(dotprod)*180./np.pi

	for i in range(h.nside2npix(nside)):
		if dtheta[i] <= (discsize-apow):
	            mask[i]=1.
	        elif dtheta[i] > (discsize-apow) and dtheta[i] <= discsize:
	            mask[i]=np.cos((dtheta[i]-(discsize-apow))*np.pi/(2.*apow))**2.

	if mp:
		h.mollview(mask,title="")
		filename = pathout + "disc_mask_size" + str(discsize) + "_apow" + str(apow) 
		filename = filename + "_nside" + str(nside) + "_mi" + str(pixel_index) + "_nsideout" + str(nsideout) + ".pdf"
		plt.savefig(filename,dpi=150,bbox_inches='tight')
		plt.close()

	return mask
		
	
	
	
	
