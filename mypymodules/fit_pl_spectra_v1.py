import healpy as h
import numpy as np
import binned_master as bm
import harm_filter as f
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.stats as stats
from matplotlib import pyplot as plt

from matplotlib import rcParams,rc
params = {'savefig.dpi': 200,
       	'axes.labelsize': 15,
	'font.size': 15,
	'axes.linewidth' : 2,
	'lines.linewidth' : 2,
	'legend.fontsize': 15,
	'xtick.labelsize': 15,
	'ytick.major.pad': 6,
	'xtick.major.pad': 6,
	'ytick.labelsize': 15,
	'text.usetex': True,
	'font.family':'sans-serif',
	'font.sans-serif':'FreeSans'}
rc('text.latex', preamble='\usepackage{sfmath}')
rcParams.update(params)


class fit_spectra(object):
	
	def __init__(self,lmin,lmax,mask,globalmask,masklmax,deltaell,blpwc=np.ones(2048,float),plot_kernel=False):
		self.lmin=lmin
		self.lmax=lmax

		# Here I ensure that the mask is multiplied by the global mask.
		self.mask=mask
		self.globalmask=globalmask
		if np.array_equal(self.mask,self.globalmask):
			self.globalmask[:]=1
		self.fsky=np.sum(self.mask*self.globalmask)/np.size(self.mask)

		self.masklmax=masklmax
		self.deltaell=deltaell
		self.blpwc=blpwc
	
		self.ell=np.arange(self.lmax+1)

		self.master=bm.binned_master(self.mask*self.globalmask,self.lmin,self.lmax,self.masklmax,self.blpwc)
		self.master.setup_binning(self.deltaell)

	def plfit(self,x,a,b):
		y=a*((x/80.)**(b+2))
		return y
	
	def plfit_sc(self,x,a):
		y=a*((x/80.)**(-2.43+2.))
		return y


	def fit_ee_spectra(self,map1,map2,clthry=np.zeros(2048,float),maskdata=True):
		"""
			Provide filtered maps to this routine
		"""

		self.fitee=np.zeros(6,float)

		if maskdata:
			alm1=h.map2alm(map1*self.mask,lmax=self.lmax)
			alm2=h.map2alm(map2*self.mask,lmax=self.lmax)
		else:
			alm1=h.map2alm(map1,lmax=self.lmax)
			alm2=h.map2alm(map2,lmax=self.lmax)	

		# Calculate the cross spectra.
		self.ubclee=self.master.return_mcs(h.alm2cl(alm1,alm2))
		self.lbin,self.clee=self.master.return_bmcs_sub_thry(h.alm2cl(alm1,alm2),clthry)

		# Evaluate the auto spectra to make the noise estimate.
		x,cl1=self.master.return_bmcs(h.alm2cl(alm1))
		x,cl2=self.master.return_bmcs(h.alm2cl(alm2))
		self.clee_err=np.sqrt(2./((2.*self.lbin+1)*self.fsky*self.deltaell))*(cl1+cl2)*0.5
		
		try:
			fit,fit_cov=curve_fit(self.plfit,self.lbin,self.clee,p0=[100,0],sigma=self.clee_err)
			self.fitee[0]=fit[0] ; self.fitee[1]=fit[1]
			self.fitee[2]=np.sqrt(np.diag(fit_cov))[0] ; self.fitee[3]=np.sqrt(np.diag(fit_cov))[1]
		except RuntimeError:
			self.fitee[0]=float('nan') ; self.fitee[1]=float('nan')
			self.fitee[2]=float('nan') ; self.fitee[3]=float('nan')

		try:
			fit_sc,fit_sc_cov=curve_fit(self.plfit_sc,self.lbin,self.clee,p0=[100],sigma=self.clee_err)
			self.fitee[4]=fit_sc[0] ; self.fitee[5]=np.sqrt(fit_sc_cov[0][0])
		except RuntimeError:
			self.fitee[5]=float('nan') ; self.fitee[5]=float('nan')		

	
		self.fitee_stat=np.zeros(6,float)

		if self.fitee[0]==self.fitee[0] and self.fitee[1]==self.fitee[1]:
			dof=len(self.lbin)-2.
			self.chi2n=stats.chi2(dof)
			self.fitee_stat[0]=np.sum(((self.clee-self.plfit(self.lbin,self.fitee[0],self.fitee[1]))**2.)/(self.clee_err**2.))
			self.fitee_stat[1]=dof
			self.fitee_stat[2]=quad(self.chi2n.pdf,0.,10.*dof)[0]-quad(self.chi2n.pdf,0,self.fitee_stat[0])[0]
		else:
			self.fitee_stat[0:3]=float("nan")

		if self.fitee[4]==self.fitee[4]:
			dof=len(self.lbin)-1.
			self.chi2n=stats.chi2(dof)
			self.fitee_stat[3]=np.sum(((self.clee-self.plfit_sc(self.lbin,self.fitee[4]))**2.)/(self.clee_err**2.))
			self.fitee_stat[4]=dof
			self.fitee_stat[5]=quad(self.chi2n.pdf,0,10.*dof)[0]-quad(self.chi2n.pdf,0,self.fitee_stat[3])[0]
		else:
			self.fitbb_stat[3:6]=float("nan")

	def fit_bb_spectra(self,map1,map2,clthry=np.zeros(2048,float),maskdata=True):
		"""
			Provide filtered maps to this routine
		"""

		self.fitbb=np.zeros(6,float)

		if maskdata:
			alm1=h.map2alm(map1*self.mask,lmax=self.lmax)
			alm2=h.map2alm(map2*self.mask,lmax=self.lmax)
		else:
			alm1=h.map2alm(map1,lmax=self.lmax)
			alm2=h.map2alm(map2,lmax=self.lmax)

		# Calculate the cross spectra.
		self.ubclbb=self.master.return_mcs(h.alm2cl(alm1,alm2))
		self.lbin,self.clbb=self.master.return_bmcs_sub_thry(h.alm2cl(alm1,alm2),clthry)

		# Evaluate the auto spectra to make the noise estimate.
		x,cl1=self.master.return_bmcs(h.alm2cl(alm1))
		x,cl2=self.master.return_bmcs(h.alm2cl(alm2))
		self.clbb_err=np.sqrt(2./((2.*self.lbin+1)*self.fsky*self.deltaell))*(cl1+cl2)*0.5
		
		try:
			fit,fit_cov=curve_fit(self.plfit,self.lbin,self.clbb,p0=[100,0],sigma=self.clbb_err)
			self.fitbb[0]=fit[0] ; self.fitbb[1]=fit[1]
			self.fitbb[2]=np.sqrt(np.diag(fit_cov))[0] ; self.fitbb[3]=np.sqrt(np.diag(fit_cov))[1]
		except RuntimeError:
			self.fitbb[0]=float('nan') ; self.fitbb[1]=float('nan')
			self.fitbb[2]=float('nan') ; self.fitbb[3]=float('nan')

		try:
			fit_sc,fit_sc_cov=curve_fit(self.plfit_sc,self.lbin,self.clbb,p0=[100],sigma=self.clbb_err)
			self.fitbb[4]=fit_sc[0] ; self.fitbb[5]=np.sqrt(fit_sc_cov[0][0])
		except RuntimeError:
			self.fitbb[5]=float('nan') ; self.fitbb[5]=float('nan')			

		self.fitbb_stat=np.zeros(6,float)

		if self.fitbb[0]==self.fitbb[0] and self.fitbb[1]==self.fitbb[1]:
			dof=len(self.lbin)-2.
			self.chi2n=stats.chi2(dof)
			self.fitbb_stat[0]=np.sum(((self.clbb-self.plfit(self.lbin,self.fitbb[0],self.fitbb[1]))**2.)/(self.clbb_err**2.))
			self.fitbb_stat[1]=dof
			self.fitbb_stat[2]=quad(self.chi2n.pdf,0,10.*dof)[0]-quad(self.chi2n.pdf,0,self.fitbb_stat[0])[0]
		else:
			self.fitbb_stat[0:3]=float("nan")

		if self.fitee[4]==self.fitee[4] :
			dof=len(self.lbin)-1.
			self.chi2n=stats.chi2(dof)
			self.fitbb_stat[3]=np.sum(((self.clbb-self.plfit_sc(self.lbin,self.fitbb[4]))**2.)/(self.clbb_err**2.))
			self.fitbb_stat[4]=dof
			self.fitbb_stat[5]=quad(self.chi2n.pdf,0,10.*dof)[0]-quad(self.chi2n.pdf,0,self.fitbb_stat[3])[0]
		else:
			self.fitbb_stat[3:6]=float("nan")


	def make_plot(self,path="./",figname="",i353=-1000,theta=-1000,clthry=np.zeros((2,2048),float),clcmb=np.zeros((2,2048),float),plot_unbinned=False,plotlog="linear",plotfit=True,plotthry=False,plotcmb=False):

		ell=np.arange(20)
		plt.figure()
		if len(figname)==0:
			figname=path + "plot.pdf"
		else:
			figname=path + figname

		
		plt.ttl=r"$\ell_{\rm min}=$" + str(self.lmin) + r" ; $\ell_{\rm max}=$" + str(self.lmax) + r" ; $\Delta \ell=$" + str(self.deltaell)
		if i353!=-1000:
			plt.ttl=plt.ttl + " ; i353=" + str(round(i353,3))
		
		if theta!=-1000:
			plt.ttl=plt.ttl + r" ; $\theta=$" + str(round(theta,1))
			
		f, (ax1, ax2) = plt.subplots(2, sharex=True) ; plt.subplots_adjust(hspace=0.2)

		if plotfit:
			ax1.plot(self.lbin,self.plfit(self.lbin,self.fitee[0],self.fitee[1]),"g-",lw=1,label=r"$A^{\rm EE}=$" + str(round(self.fitee[0],2)) + "$\pm$" + str(round(self.fitee[2],2)) + r"; $\alpha_E=$" + str(round(self.fitee[1],2)) + "$\pm$" + str(round(self.fitee[3],2)) + r"; $\chi^2=$" + str(round(self.fitee_stat[0],2)) + "; PTE=" + str(round(self.fitee_stat[2],3)))
		        ax1.plot(self.lbin,self.plfit_sc(self.lbin,self.fitee[4]),"r--",lw=1,label=r"$A^{\rm EE}=$" + str(round(self.fitee[4],2)) + "$\pm$" + str(round(self.fitee[5],2)) + r"; $\alpha^0_E=-2.42$" + r"; $\chi^2=$" + str(round(self.fitee_stat[3],2)) + "; PTE=" + str(round(self.fitee_stat[5],3)))

		if not(np.array_equal(clthry[0][:],np.zeros(2048,float))) and plotthry:
			ax1.plot(self.ell,self.ell*(self.ell+1)*clthry[0][:self.lmax+1]/(2.*np.pi),"k-",lw=1,alpha=0.5,label=r"$C_{\ell}^{\rm EE}$ dust $({\rm A}^{EE}=197.1$ ; ${\alpha}^{EE}=-2.42)$")

		if not(np.array_equal(clcmb[0][:],np.zeros(2048,float))) and plotcmb:
			ax1.plot(self.ell,self.ell*(self.ell+1)*clcmb[0][:self.lmax+1]/(2.*np.pi),"m-",lw=1,alpha=0.5,label=r"$C_{\ell}^{\rm EE}$ CMB")

	        ax1.errorbar(self.lbin,self.clee,yerr=self.clee_err,linestyle="none",marker='o',ms=3,mec="none",mfc="b",elinewidth=1,label=r"Binned $C_{\ell}^{\rm EE}$")
		if plot_unbinned:
			ax1.plot(self.ell,self.ell*(self.ell+1)*self.ubclee/(2.*np.pi),"c-",lw=1,alpha=0.3,label="")

	        ax1.hlines(0,0,max(self.lbin)+100,linestyle="dashed")

		#----------------------
		if plotfit:	
			ax2.plot(self.lbin,self.plfit(self.lbin,self.fitbb[0],self.fitbb[1]),"g-",lw=1,label=r"$A^{\rm BB}=$" + str(round(self.fitbb[0],2)) + "$\pm$" + str(round(self.fitbb[2],2)) + r"; $\alpha_B=$" + str(round(self.fitbb[1],2)) + "$\pm$" + str(round(self.fitee[3],2))  + r"; $\chi^2=$" + str(round(self.fitbb_stat[0],2)) + "; PTE=" + str(round(self.fitbb_stat[2],3)))
		        ax2.plot(self.lbin,self.plfit_sc(self.lbin,self.fitbb[4]),"r--",lw=1,label=r"$A^{\rm BB}=$" + str(round(self.fitbb[4],2)) + "$\pm$" + str(round(self.fitbb[5],2)) + r"; $\alpha^0_B=-2.42$" + r"; $\chi^2=$" + str(round(self.fitbb_stat[3],2)) + "; PTE=" + str(round(self.fitbb_stat[5],3)))

		if not(np.array_equal(clthry[1][:],np.zeros(2048,float))) and plotthry:
			ax2.plot(self.ell,self.ell*(self.ell+1)*clthry[1][:self.lmax+1]/(2.*np.pi),"k-",lw=1,alpha=0.5,label=r"$C_{\ell}^{\rm BB}$ dust $({\rm A}^{BB}=104.5$ ; ${\alpha}^{BB}=-2.44)$")

		if not(np.array_equal(clcmb[1][:],np.zeros(2048,float))) and plotcmb:
			ax2.plot(self.ell,self.ell*(self.ell+1)*clcmb[1][:self.lmax+1]/(2.*np.pi),"m-",lw=1,alpha=0.5,label=r"$C_{\ell}^{\rm BB}$ CMB")

	        ax2.errorbar(self.lbin,self.clbb,yerr=self.clbb_err,linestyle="none",marker='o',ms=3,mec="none",mfc="b",elinewidth=1,label=r"Binned $C_{\ell}^{\rm BB}$")
		if plot_unbinned:
			ax2.plot(self.ell,self.ell*(self.ell+1)*self.ubclbb/(2.*np.pi),"c-",lw=1,alpha=0.3,label="")
			
	        ax2.hlines(0,0,max(self.lbin)+100,linestyle="dashed")

		ax1.title.set_text(plt.ttl)
		ax1.legend(loc=0,fontsize=6)
	        ax2.legend(loc=0,fontsize=6)

		if plotlog=="logy":
			ax1.semilogy()
			ax2.semilogy()
		elif plotlog=="logx":
			ax1.semilogx()
			ax2.semilogx()
		elif plotlog=="loglog":
			ax1.loglog()
			ax2.loglog()

        	plt.xlabel("multipole $\ell$")
        	plt.ylabel(r"$\ell(\ell+1)C_{\ell}/(2 \pi)$")
		
		plt.savefig(figname,dpi=150,bbox_inches="tight")
		plt.close()


	def plot_kernel(self,pathout,maskindex,plot_binned=True):
		plt.imshow(np.log10(self.master.mllp))
		plt.colorbar()
		figname=pathout + "/kernel_mi" + str(maskindex) + "_lmin" + str(self.lmin) + "_lmax" + str(self.lmax) + "_dell" + str(self.deltaell) + ".pdf"
		plt.savefig(figname,dpi=150,bbox_inches="tight")
		plt.close()

		if plot_binned:
			plt.imshow(np.log10(self.master.mbbp))
			plt.colorbar()
			figname=pathout + "/kernel_binned_mi" + str(maskindex) + "_lmin" + str(self.lmin) + "_lmax" + str(self.lmax) + "_dell" + str(self.deltaell) + ".pdf"
			plt.savefig(figname,dpi=150,bbox_inches="tight")
			plt.close()

