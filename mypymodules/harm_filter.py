import healpy as h
import numpy as np

class harmonic_filter(object):

	def __init__(self,lmin,lmax,maplmax,apow):
		self.lmin=lmin
		self.lmax=lmax
		self.maplmax=maplmax
		self.apow=apow

		self.ell=np.linspace(0,self.maplmax,self.maplmax+1)
		self.tfilter=np.zeros(self.ell.size,float)
		self.pfilter=np.zeros(self.ell.size,float)

		for i in range(self.ell.size):
			if i<self.lmin or i>self.lmax:
				self.tfilter[i]=0.
			elif i>=self.lmin+self.apow and i<=self.lmax-self.apow:
				self.tfilter[i]=1.
			elif i<self.lmax-self.apow:
				self.tfilter[i]=np.cos((self.ell[i]-self.lmin-self.apow)*np.pi/(2.*self.apow))**2.
			elif i>self.lmin+2.*self.apow:
				self.tfilter[i]=np.cos((self.ell[i]-self.lmax+self.apow)*np.pi/(2.*self.apow))**2.


		self.pfilter=self.tfilter*np.sqrt(1./((self.ell+2.)*(self.ell+1)*(self.ell-1)*self.ell))
		self.pfilter[0:2]=0.
			
		
	def return_filtered_map(self,xmap,pol):
		"""
		- Removes the prefactor applied to the polarization maps which evaluating the E/B maps.
                - Also does a filtering in harmonic space, using the apodization filters.
		"""
		temp_filter=self.tfilter
		if pol:
			temp_filter=self.pfilter
			
		alm=h.map2alm(xmap,lmax=self.lmax)
		for i in range(alm.size):
			ell=h.Alm.getlm(self.lmax,i)[0]
			alm[i]=alm[i]*temp_filter[ell]

		fmap=h.alm2map(alm,nside=h.get_nside(xmap),verbose=False)

		return fmap

	def return_filtered_alm(self,xmap,pol):
		"""
		- Removes the prefactor applied to the polarization maps which evaluating the E/B maps.
                - Also does a filtering in harmonic space, using the apodization filters.
		"""
		temp_filter=self.tfilter
		if pol:
			temp_filter=self.pfilter
			
		alm=h.map2alm(xmap,lmax=self.lmax)
		for i in range(alm.size):
			ell=h.Alm.getlm(self.lmax,i)[0]
			alm[i]=alm[i]*temp_filter[ell]

		return alm

	def return_filtered_auto_cl(self,xmap,pol):
		"""
		- Removes the prefactor applied to the polarization maps which evaluating the E/B maps.
                - Also does a filtering in harmonic space, using the apodization filters.
		"""
		temp_filter=self.tfilter
		if pol:
			temp_filter=self.pfilter
			
		alm=h.map2alm(xmap,lmax=self.lmax)
		for i in range(alm.size):
			ell=h.Alm.getlm(self.lmax,i)[0]
			alm[i]=alm[i]*temp_filter[ell]

		cl=h.alm2cl(alm)

		return cl


	def return_filtered_cross_cl(self,xmap1,xmap2,pol):
		"""
		- Removes the prefactor applied to the polarization maps which evaluating the E/B maps.
                - Also does a filtering in harmonic space, using the apodization filters.
		"""
		temp_filter=self.tfilter
		if pol:
			temp_filter=self.pfilter
			
		alm1=h.map2alm(xmap1,lmax=self.lmax)
		alm2=h.map2alm(xmap2,lmax=self.lmax)

		for i in range(alm1.size):
			ell=h.Alm.getlm(self.lmax,i)[0]
			alm1[i]=alm1[i]*temp_filter[ell]
			alm2[i]=alm2[i]*temp_filter[ell]

		cl=h.alm2cl(alm1,alm2)
		return cl
