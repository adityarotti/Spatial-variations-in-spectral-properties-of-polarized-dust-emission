import bmaster as ms
import numpy as np
import healpy as h
from matplotlib import pyplot as plt

# Assumes that the spectra provided have l=0,1 included. Excludes those explicitly by starting from the second element of the array.

class binned_master(object):
	
	def __init__(self,mask,lmin,lmax,masklmax,beam=np.ones(4096,float)):
		self.mask=mask
		self.lmin=max(lmin,2)
		self.lmax=lmax
		self.masklmax=masklmax
		self.beam=beam
		
		self.clmask=h.alm2cl(h.map2alm(self.mask,lmax=self.masklmax))
		self.mllp=ms.master.calc_kernel(self.clmask,self.lmin,self.lmax,self.masklmax)

		self.pbl=[]
		self.qlb=[]
		self.lbin=[]
		self.deltaell_bin=[]
		self.mbbp=[]

	# Defining the projection and the deprojection operators
	def setup_binning(self,deltaell):

		totell=self.lmax-self.lmin+1
		nbin=np.int(totell/deltaell)
		self.pbl=np.zeros((nbin,totell),float)
		self.qlb=np.zeros((totell,nbin),float)
		self.qlb_nobeam=np.zeros((totell,nbin),float)
		self.lbin=[]
		self.deltaell_bin=[]

		for i in range(nbin):
			bmin=i*deltaell
			bmax=min(bmin+deltaell-1,self.lmax-self.lmin)
			ell=np.linspace(bmin+self.lmin,bmax+self.lmin,bmax-bmin+1)
			ellmin=min(ell) ; ellmax=max(ell) 
			temp_bl=self.beam[int(ellmin):int(ellmax)+1]
			self.lbin=np.append(self.lbin,int(np.mean(ell)))
			norm=len(ell)
			f1=ell*(ell+1)/(2.*np.pi*norm)
			g1=2.*np.pi/(ell*(ell+1)) 
			self.pbl[i,bmin:bmax+1]=f1
			self.qlb[bmin:bmax+1,i]=g1*(temp_bl**2.)
			self.qlb_nobeam[bmin:bmax+1,i]=g1
			self.deltaell_bin=np.append(norm,self.deltaell_bin)
						
		self.mbbp=np.array(np.matrix(self.pbl)*np.matrix(self.mllp)*np.matrix(self.qlb))
		
	def return_mcs(self,cl):
		"""
		Returns the master corrected spectrum (mcs)
		"""
		mcl=ms.master.est_true_cl(cl[self.lmin:self.lmax+1],self.mllp,len(cl[self.lmin:self.lmax+1]))
		mcl=np.append(np.zeros(self.lmin,float),mcl)
		mcl=mcl/(self.beam[:self.lmax+1]**2.)
		return mcl


	def return_bmcs(self,cl):
		"""
		Returns the binned master corrected spectrum (bmcs)
		"""
		bcl=np.array(np.matrix(self.pbl)*np.transpose(np.matrix(cl[self.lmin:self.lmax+1])))[:,0]
		bcl=ms.master.est_true_cl(bcl,self.mbbp,len(bcl))

		#ubcl=np.array(np.matrix(self.qlb_nobeam)*np.transpose(np.matrix(bcl)))[:,0]
		#ubcl=np.append(np.zeros(self.lmin,float),ubcl)

		return self.lbin,bcl #,ubcl

	def return_bmcs_sub_thry(self,cl,clthry=np.zeros(2048,float)):
		"""
		Returns the binned master corrected spectrum (bmcs)
		"""
		bcl=np.array(np.matrix(self.pbl)*np.transpose(np.matrix(cl[self.lmin:self.lmax+1])))[:,0]
		bcl=ms.master.est_true_cl(bcl,self.mbbp,len(bcl))
		bclthry=np.array(np.matrix(self.pbl)*np.transpose(np.matrix(clthry[self.lmin:self.lmax+1])))[:,0]
		bcl=bcl-bclthry

		#ubcl=np.array(np.matrix(self.qlb_nobeam)*np.transpose(np.matrix(bcl)))[:,0]
		#ubcl=np.append(np.zeros(self.lmin,float),ubcl)

		return self.lbin,bcl #,ubcl


	def return_binned_spectra(self,cl):
		bcl=np.array(np.matrix(self.pbl)*np.transpose(np.matrix(cl[self.lmin:self.lmax+1])))[:,0]
		return self.lbin,bcl

	def plot_mask(self,maskindex,pathout):
		h.orthview(self.mask,title="",rot=(0,90))
		filename = pathout + "disc_mask" + "_mi" + str(maskindex) + ".pdf"
		plt.savefig(filename,dpi=150,bbox_inches='tight')
		plt.close()




	
		

