import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

dt=10.
dz=100.
n=120
zpeak1=5000.
zstdev1=500.
peak1=1e-3
zpeak2=5000.
zstdev2=500.
peak2=1e3

# parameters for rain
a1=4854.0
a2=-446.0
b1=1.0
b2=0.782
f=0.5
g1=0.0
g2=4085.35
c=522.0
d=3.0
mu=2.5
initialise=True
def initialisation1():
	z=np.mgrid[0:dz*n:dz]
	
	q=peak1*np.exp(-0.5*((z-zpeak1)/zstdev1)**2)
	n1=peak2*np.exp(-0.5*((z-zpeak2)/zstdev2)**2)
	
		
	return (z,q,n1)

def calcparameters(q,n1):
	lambda0 = (c*gamma(mu+d+1.0)*n1/(gamma(mu+1.0)*q))**(1.0/d)
	n0 = n1*lambda0**(mu+1.0) / gamma(mu+1.0)
	
	return (lambda0,n0)


def calcfallspeeds(lambda0,n0):
	vq = a1*gamma(b1+d+mu+1.0) / lambda0**b1 + a2*gamma(b2+d+mu+1.0) / lambda0**b2
	vq = vq / gamma(d+mu+1.0)
	# calculate the number weighted fall-speed
	vn = a1*gamma(b1+mu+1.0) / lambda0**b1 + a2*gamma(b2+mu+1.0) / lambda0**b2
	vn = vn / gamma(mu+1.0)
	
	return (vq,vn)

def sedimentation1(dt,dz,z,q,n1):

	# calculate parameters
	(lambda0,n0)=calcparameters(q,n1)
	
	# calculate the mass weighted fall-speed
	(vq,vn)=calcfallspeeds(lambda0,n0)
	
	# rotstayn solution
	qold=q.copy()
	nold=n1.copy()
	for i in range(n-1):
		# flux top - q
		rf = qold[i+1]*vq[i+1]
		a=vq[i]*dt/dz
		q[i] = qold[i]*np.exp(-a) + rf/vq[i]*(1.0-np.exp(-a))
		# flux top - n
		a=vn[i]*dt/dz
		rf = nold[i+1]*vn[i+1]
		a=vn[i]*dt/dz
		n1[i] = nold[i]*np.exp(-a) + rf/vn[i]*(1.0-np.exp(-a))
	return (z,q,n1,vq,vn)


def sedimentation2(dt,dz,z,q,n1):
	# this is for upwind method
	# calculate parameters
	(lambda0,n0)=calcparameters(q,n1)
	
	# calculate the mass weighted fall-speed
	(vq,vn)=calcfallspeeds(lambda0,n0)
	ind,=np.where(q<1e-8)
	ind2,=np.where(n1<1e-3)

	vmax=np.max(np.maximum(np.max(np.abs(vq[ind])),np.max(np.abs(vn[ind2]))))
	dt1 = 0.5*dz/vmax
	nsubsteps=int(np.maximum(1.0,np.ceil(dt/dt1)))
	dt1=dt/nsubsteps
	print(nsubsteps)

	# upwind solution
	fz_r=np.zeros(n)
	fz_l=np.zeros(n)
	fz_r1=np.zeros(n)
	fz_l1=np.zeros(n)
	for j in range(nsubsteps):
		for i in range(n-1):
			fz_r[i]= ( (-vq[i+1]+np.abs(-vq[i+1]))*q[i] + \
					   (-vq[i+1]-np.abs(-vq[i+1]))*q[i+1] )*dt1/(2.0*dz)
			fz_l[i+1]= ( (-vq[i+1]+np.abs(-vq[i+1]))*q[i] + \
					   (-vq[i+1]-np.abs(-vq[i+1]))*q[i+1] )*dt1/(2.0*dz)
					   
			fz_r1[i]= ( (-vn[i+1]+np.abs(-vn[i+1]))*n1[i] + \
					   (-vn[i+1]-np.abs(-vn[i+1]))*n1[i+1] )*dt1/(2.0*dz)
			fz_l1[i+1]= ( (-vn[i+1]+np.abs(-vn[i+1]))*n1[i] + \
					   (-vn[i+1]-np.abs(-vn[i+1]))*n1[i+1] )*dt1/(2.0*dz)
		fz_l[0]=fz_r[0]
		fz_l1[0]=fz_r1[0]
		fz_r[-1]=fz_l[-1]
		fz_r1[-1]=fz_l1[-1]
		q=q-(fz_r-fz_l)
		n1=n1-(fz_r1-fz_l1)

	return (z,q,n1,vq,vn)


if __name__=="__main__":
	plt.ion()

	if initialise:
		(z,q,n1)=initialisation1()
		print('n,q: ' + str(np.sum(n1)) + ', ' + str(np.sum(q)))
		
	plt.subplot(121)
	plt.plot(q,z)	
	plt.xlabel('q$_{rain}$')
	plt.ylabel('z (m)')
	plt.subplot(122)
	plt.plot(n1,z)
	plt.xlabel('n$_{rain}$')
	plt.ylabel('z (m)')
	
	for i in range(60):
		(z,q,n1,vq,vn)=sedimentation1(dt,dz,z,q,n1)
	plt.subplot(121)
	plt.plot(q,z)	
	plt.subplot(122)
	plt.plot(n1,z)
	
	(z,q,n1)=initialisation1()
	for i in range(60):
		(z,q,n1,vq,vn)=sedimentation2(dt,dz,z,q,n1)
	plt.subplot(121)
	plt.plot(q,z)	
	plt.subplot(122)
	plt.plot(n1,z)
	
	
# 	plt.subplot(133)
# 	plt.plot(vq,z)
# 	plt.plot(vn,z)
	