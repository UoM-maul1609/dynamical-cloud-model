from scipy.special import gamma
import numpy as np
import matplotlib.pyplot as plt


if __name__=="__main__":
	mu1=20
	q=0.3e-3
	n=100.0e6
	c=np.pi/6.0*1000.0
	d=3
	
	lambda1=(c*n/q*gamma(mu1+1+d)/gamma(mu1+1))**(1.0/d)
	n0=n*lambda1**(mu1+1.0)/gamma(mu1+1)
	# calculate m1,m2 and m3, for real data you 
	# would need to do sum(Ni*di), sum(Ni*di**2), sum(Ni*di**6)
	M1=n0*gamma(mu1+2) / lambda1**(mu1+2)	
	M2=n0*gamma(mu1+3) / lambda1**(mu1+3)
	M6=n0*gamma(mu1+7) / lambda1**(mu1+7)
	
	F=M2**5/(M6*M1**4)
	
	# see https://www.aoml.noaa.gov/ftp/hrd/rblack/ForGreg/Heymsfield_et_al_2002.pdf
	# equation 5
	p=[1.0-F, 8.0-18.0*F,24.0-119.0*F,\
		32.0-342.0*F,16.0-360.0*F]
	roots=np.roots(p)
	
	mu2=roots
	lambda2=M1/M2*gamma(mu2+3)/gamma(mu2+2)	
	n02=M1/gamma(mu2+2)*lambda2**(mu2+2)
	ind,=np.where(lambda2>0)
	
	lambda2=np.real(lambda2[ind[0]])
	n02=np.real(n02[ind[0]])
	mu2=np.real(mu2[ind[0]])
	
	d=np.linspace(0,100e-6,1000)
	plt.ion()
	plt.plot(d,n0*d**(mu1)*np.exp(-lambda1*d))
	