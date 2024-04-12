import numpy as np
from scipy.integrate import odeint
from scipy.special import gamma
import distutils

"""
	this function calculates the autoconversion of ice to snow
"""
def derivs(y, t, lambdaimin,ci,Dis, di, tsaut, mui, cs,ds):
	qi = y[2]
	ni = y[3]
	qs = y[0]
	ns = y[1]
	
	lami = (ci*gamma(1.0 +mui+ di)/gamma(1.0 + mui)* \
		ni/np.maximum(qi,1.e-30))**(1.0/di)
	
	
	# equation A7 in QJRMS, Field et al. (2023)
	if (lami <= lambdaimin):
		psaut = qi/tsaut * ((lambdaimin/lami)**di-1.0)
		nsaut = psaut/(ci*Dis**di)	
		niaut = psaut/(ci*Dis**di) #-nsaut #-psaut*ni/np.maximum(qi,1e-30)
	else:
		psaut=0.0
		nsaut=0.0
		niaut=0.0
			
	if((ni < 1.0) | (qi < 1.e-20)):
		psaut=0.0
		nsaut=0.0
		niaut=0.0

	# qs, ns, qi, ni
	dydt = [psaut, nsaut, -psaut, niaut]
	
	return dydt




if __name__=="__main__":
	import matplotlib.pyplot as plt
	import sys
	
	old_um_flag = bool(distutils.util.strtobool(sys.argv[1]))
# 	old_um_flag = True
	N=10e3		# initial ice number concentration
	Q=1.0e-3		# initial ice mixing ratio
	
	mui = 0.0		# shape parameter of ice
	mus = 0.0		# shape parameter of snow
	ci=200.0*np.pi/6.0	# mass dimension prefactor
	di=3.0				# mass-dimension exponent
	cs=0.026
	ds=2.0
	
	Dis = 100.0e-6		# this size that newly formed snow particles are
	Dimax = 100.0e-6	# max average size of ice before converting to snow
	if old_um_flag:
		Dimax = 50.0e-6	# max average size of ice before converting to snow
		Dis = 50.0e-6	# max average size of ice before converting to snow
		lambdaimin = (1.0 + mui)/ Dimax # corresponding lambda
	else:
		lambdaimin = (1.0 + mui +di)/ Dimax # corresponding lambda
	
	tsaut = 60.0	# autoconversion time-scale
	
	y0 = [0.0,0.0, Q, N]
	t = np.linspace(0,600.,1001)
	sol = odeint(derivs, y0, t, args=(lambdaimin, ci, Dis, di, tsaut, mui, cs, ds))


	plt.ion()
	plt.figure()
	qs=sol[:,0]
	ns=sol[:,1]
	qi=sol[:,2]
	ni=sol[:,3]
	lami = (ci*gamma(1.0 +mui+ di)/gamma(1.0 + mui)* \
		ni/np.maximum(qi,1.e-30))**(1.0/di)
	lams = (cs*gamma(1.0 +mus+ ds)/gamma(1.0 + mus)* \
		ns/np.maximum(qs,1.e-30))**(1.0/ds)
	
	if old_um_flag:
		diam_ice=(1.0 + mui)/lami
		diam_snow=(1.0 + mus)/lams
	else:
		diam_ice=(1.0 + mui + di)/lami
		diam_snow=(1.0 + mus + ds)/lams
	plt.subplot(211)
	plt.plot(t, Dimax*np.ones(len(t))*1e6)
	plt.plot(t, diam_ice*1e6)
	plt.ylabel('$D_{ice}$ ($\mu$m)')
	plt.xlabel('time (s)')
	plt.legend(['Threshold to start auto','Actual ice diameter'])
	if old_um_flag:
		plt.title('Old UM')
	else:
		plt.title('New UM')
	
	
	plt.subplot(212)
	plt.plot(t, Dimax*np.ones(len(t))*1e6)
	plt.plot(t, diam_snow*1e6)
	plt.ylabel('$D_{snow}$ ($\mu$m)')
	plt.xlabel('time (s)')
	plt.legend(['Threshold to start auto','Actual snow diameter'])
	
	n0s=ns*lams**(1.0+mus)/gamma(1.0+mus)
	zfactor = (6.0*cs/np.pi/1000.)**2*n0s*gamma(1.0+mus+2.0*ds) / \
		lams**(1.0+mus+2.0*ds) / (1e-3**6)
	if old_um_flag:
		fig=plt.figure()
		
		plt.plot(t,np.log10( zfactor ))
		#plt.plot(t,ns)
	else:
		plt.figure(fig)
		plt.plot(t,np.log10( zfactor ),'--')
		#plt.plot(t,ns,'--')
		plt.legend(['old UM','new UM'])
		plt.xlabel('t (s)')
		plt.ylabel('Z')	