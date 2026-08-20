import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.interpolate import interp1d

if __name__=="__main__":
	lambdas=np.logspace(-8,1,10000)
	result=np.zeros(len(lambdas))
	result2=np.zeros(len(lambdas))
	 
	qc=0.3e-3		# liquid water content
	mu=7			# shape parameter
	Dm=20e-6		# mean diameter
	lam=(mu+1)/Dm	# shape parameter
	n0=qc*lam**(mu+4)/gamma(mu+4)/(np.pi/6*1000)	# intercept
	mtun=0.5		# page 1315, Mitchell et al.
	
	import refractive_indices_s
	"""
	import refractive_indices as refractive_indices_s
	refractive_indices_s.lamr_h2o=refractive_indices_s.lam_h2o
	refractive_indices_s.lami_h2o=refractive_indices_s.lam_h2o
	refractive_indices_s.nr_h2o=refractive_indices_s.nr_h2o
	refractive_indices_s.ni_h2o=refractive_indices_s.ni_h2o
	"""
	scint_nrwbin=interp1d(refractive_indices_s.lamr_h2o,refractive_indices_s.nr_h2o)
	scint_niwbin=interp1d(refractive_indices_s.lami_h2o,refractive_indices_s.ni_h2o)
	for i in range(len(lambdas)):
		lambda1 = lambdas[i]
		if(lambda1 < np.min(refractive_indices_s.lamr_h2o)):
			nrwbin=refractive_indices_s.nr_h2o[0] 
		elif(lambda1 > np.max(refractive_indices_s.lamr_h2o)):
			nrwbin=refractive_indices_s.nr_h2o[-1] 
		else:
			nrwbin=scint_nrwbin(lambda1)
			
		if(lambda1 < np.min(refractive_indices_s.lami_h2o)):
			niwbin=refractive_indices_s.ni_h2o[0]
		elif(lambda1 > np.max(refractive_indices_s.lami_h2o)):
			niwbin=refractive_indices_s.ni_h2o[-1] 
		else:
			niwbin=scint_niwbin(lambda1)

		sigma=np.pi/4
		delta=2
		g=8*np.pi*niwbin/(3*lambda1)
		e0=0.25+0.6*(1.-np.exp(-8.*np.pi*niwbin/3.))**2 # eq 8
		e=e0/lambda1 # equation 28
		ra=0.7393*nrwbin-0.6069 # equation 6
		rext=ra*0.5 # equation 11
		a1=0.25+0.25*np.exp(-1167.*niwbin) # equation5
		kmax=mtun/e0 # equation 9
		a2=ra/(kmax**mtun*np.exp(-mtun)*lambda1**mtun)
		a3=rext/(kmax**mtun*np.exp(-mtun)*lambda1**mtun)
		a4=0.06*np.pi/lambda1
		a5=(np.pi/lambda1)**(-2./3.)
		a6=1.0
	
		gam_a_1=gamma(mu+3)
		gam_a_2=gamma(3+mtun+mu)
	
		# absorption
		babs= \
			sigma*n0 * \
				gam_a_1 / \
				(lam**(3.+mu)) - \
			sigma*n0 * \
				gam_a_1 / \
				((lam+g)**(3.+mu)) + \
			a1*sigma*n0 * \
				gam_a_1 / \
				((lam+g)**(3.+mu)) - \
			a1*sigma*n0 * \
				gam_a_1 / \
				((lam+2.*g)**(3.+mu)) + \
			a2*sigma*n0 * \
				gam_a_2 / \
				((lam+e)**(3.+mu+mtun)) - \
			a2*sigma*n0 * \
				gam_a_2 / \
				((lam+e+g)**(3.+mu+mtun))
		m2=sigma*n0*gamma(mu+3)/lam**(mu+3)
		gam_3=gamma(mu+7/3)
		gam_4=gamma(mu+2)
		gam_5=gamma(mu+1)
		gam_6=gamma(mtun+mu+2)
		gam_7=gamma(mtun+mu+1)
		imag1=complex(0,1)
		ncom=complex(nrwbin,-niwbin)
		q=imag1*(ncom-complex(1,0))*complex(2*np.pi/lambda1,0)
		
		# extinction
		test= \
			np.pi*n0*gam_a_1 / (2.*lam**(3.+mu)) + \
			a3*np.pi*n0*gam_a_2 / (2.*(lam+e)**(3.+mu+mtun)) + \
			a6*sigma*a5*n0*gam_3* \
				(lam**(-(mu+7./3.))-(lam+a4)**(-(mu+7./3.)))+ \
			np.pi*n0*np.real( \
			 complex(gam_4,0.) / (q*(complex(lam,0.)+q)**(complex(mu+2.,0.))) + \
			 complex(gam_5,0.) / (q*q)*((complex(lam,0.)+q)**(-(complex(mu+1.,0.))) - \
				complex(lam,0.)**(-complex(mu+1.,0.))))+ \
			a3*np.pi*n0*np.real( \
			 complex(gam_6,0.) / (q*(complex(lam+e,0.)+q)** complex(mu+mtun+2.,0.)) + \
			 complex(gam_7,0.) / (q*q)*((complex(lam+e,0.)+q)**(-complex(mu+mtun+1.,0.)) - \
			   complex(lam+e,0.)**(-complex(mu+mtun+1.,0.))))
	
		result[i]=babs  #/m2;
		result2[i]=test  #/m2;
		result2[i]=np.maximum(result2[i],result[i])
	# create a plot
	plt.ion()
	plt.figure()
	plt.subplot(211)
	# total extinction
	plt.plot(3e8/lambdas,1-np.exp(-((result2)*200)))
	# absorption
	plt.plot(3e8/lambdas,1-np.exp(-((result)*200)))
	
	plt.xlabel("frequency, $\\nu$ (s$^{-1}$)")
	plt.ylabel('Total amount of light extinction or absorption')
	plt.legend(['Total extinction','Absorption'])
	plt.xscale('log')
	plt.title('Cloud of thickness 200m (qc=' + str(qc*1000) + ' g kg$^{-1}$, ' + \
		str(Dm*1e6) + ' $\\mu$m diameter drops)')
	
	plt.subplot(212)
	# real refractive index
	plt.plot(3e8/np.array(refractive_indices_s.lamr_h2o),refractive_indices_s.nr_h2o)
	# absorption
	plt.plot(3e8/np.array(refractive_indices_s.lami_h2o),refractive_indices_s.ni_h2o)
	
	plt.xlabel("frequency, $\\nu$ (s$^{-1}$)")
	plt.ylabel('nr, ni')
	plt.legend(['nr','ni'])
	plt.xscale('log')
	plt.yscale('log')
	
	
	