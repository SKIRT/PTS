from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import random






# fit two blackbodies to the synthetic data



def funcinv(wa,D,f,T1):
    kv=k850*(850/wa)**2
    Md=f/(kv/D**2*blackbody_lam(wa, T1)*2.08*10**11)#*1.98*10**30/9.5/10**32/10**12*10**26	
    return np.log10(Md)







for i in range(len(data[:,0])):

	ydata = [data[i,30],data[i,34],data[i,36],data[i,38],data[i,40],data[i,42]]
	yerr = [data[i,31],data[i,35],data[i,37],data[i,39],data[i,41],data[i,43]]
	D = data[i,1]*3*10**5/67.30

	chi2 = float('inf')
	T1dist = np.arange(10,30,3)
	T2dist = np.arange(30,60,10)
	Mddist = np.arange(4,6,0.02)
	Mdprob = np.zeros(len(Mddist))
	ptot = 0
	for ii in range(len(Mddist)):
		for T2 in T2dist:
			for T1 in T1dist:
				popt = minimize(leastsq, [0.05], args=(Mddist[ii],wa,ydata,yerr,D,T1,T2),method='Nelder-Mead',options={'maxiter':200})
				Mdr_new = popt.x
				chi2_new = leastsq(Mdr_new,Mddist[ii],wa,ydata,yerr,D,T1,T2)
				if chi2 > chi2_new:
					chi2 = chi2_new
					Mdratio_sav = Mdr_new
					T1_sav = T1
					T2_sav = T2
					Md_sav = Mddist[ii]
					print Mddist[ii],Mdr_new,T1_sav,T2_sav
				prob=np.exp(-0.5*chi2_new)
				Mdprob[ii]+=prob
				ptot+=prob
	Mdprob = Mdprob / ptot
	print "tcold", Mddist[np.where(Mdprob==max(Mdprob))], percentiles(Mddist,Mdprob,16),percentiles(Mddist,Mdprob,50),percentiles(Mddist,Mdprob,84)
	plt.plot(Mddist, Mdprob)
	plt.title(dataID[i])
	plt.xlabel('T_cold (K)')
	plt.show()
	print Md_sav, T1_sav, T2_sav
	plt.errorbar(wa, ydata, yerr=yerr, fmt='bo')
	x2=np.arange(10.,600.,0.1)
	plt.loglog()
	plt.title(dataID[i])
	plt.ylim(0.001,1)
	plt.plot(x2, func2(x2,D,np.log10(Mdratio_sav)+Md_sav,T1_sav,np.log10(Mdratio_sav)+Md_sav,T2_sav), 'r', label="best fit")
	plt.xlabel('Wavelength (microns)')
	plt.ylabel('Flux (Jy)')	
	plt.plot(x2, func(x2,D,np.log10(Mdratio_sav)+Md_sav,T1_sav), ':', lw=2,label="cold dust: logMd = %s, Tc= %s K "%(np.log10(Mdratio_sav)+Md_sav,T1_sav) )
	plt.plot(x2, func(x2,D,np.log10(Mdratio_sav)+Md_sav,T2_sav), ':', lw=2,label="warm dust: logMd = %s, Tc= %s K "%(np.log10(Mdratio_sav)+Md_sav,T2_sav) )
	plt.legend(loc=4)
	plt.show()
	
plt.legend(frameon=False)
plt.savefig('fit_bb.png')
