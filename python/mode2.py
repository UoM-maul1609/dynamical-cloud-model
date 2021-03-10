import numpy as np
from scipy import integrate

"""
some constants 
"""
TTR=273.15
Lf=333550.
cw=4200.
gamma_liq=0.072 # surface tension of water
DEcrit=0.2 # see page 3039 of Phillips et al.
phi=0.5
T=265.
f=-cw*(T-TTR)/Lf # fraction frozen at the end of stage 1 of freezing
alphar=2.5
alphai=2.5


def test(y,x):
    return x*y**2

"""
Collisions between precipitating particles of different species - see page 27 of LEM
-need to calculate the number of rain drops colliding onto ice - a sink of rain
-need to calculate the mass of rain drops colliding onto ice - a sink of rain and source
- of ice mass
"""
def collisions_precip(diamr,diami,n0r,lambda0r,n0i,lambda0i,ar,br,ai,bi,cr,dr,ci,di,f):
    # number of collisions between rain and ice
    integrand1=np.pi/4.*(diamr+diami)**2*np.abs(ar*diamr**br-ai*diami**bi)* \
        n0r*diamr**alphar*np.exp(-lambda0r*diamr)*n0i*diami**alphai*np.exp(-lambda0i*diami)
        
    # calculate number of fragments for each of these collisions
    mr=cr*diamr**dr  
    mi=ci*diami**di  
    (nfrag,nfrag_freeze1,nfrag_freeze2,nfrag_liq)= \
        mode2_phillips(f,mr,mi,ar,br,ai,bi,cr,dr,ci,di)
        
    # number of new ice from mode 2
    #integrand1=integrand1*nfrag_freeze2

    # mass of rain collected by ice  - not sure equation 112 is correct way around?
#    integrand2=np.pi/4.*(diamr+diami)**2*np.abs(ar*diamr**br+ai*diami**bi)* \
#        cr*diamr**dr*n0r*np.exp(-lambda0r*diamr)*n0i*np.exp(-lambda0i*diami)
        
    #return (integrand1,integrand2)
    return integrand1
    
"""
Collisions between precipitating particles of different species - see page 27 of LEM
-need to calculate the number of rain drops colliding onto ice - a sink of rain
-need to calculate the mass of rain drops colliding onto ice - a sink of rain and source
- of ice mass
"""
def collisions_precip_mass(mr,mi,n0r,lambda0r,n0i,lambda0i,ar,br,ai,bi,cr,dr,ci,di,f):
    # number of collisions between rain and ice
    diamr=(mr/cr)**(1/dr)
    diami=(mi/ci)**(1/di)
    integrand1=np.pi/4.*(diamr+diami)**2*np.abs(ar*diamr**br-ai*diami**bi)* \
        n0r*diamr**alphar*np.exp(-lambda0r*diamr)*n0i*diami**alphai*np.exp(-lambda0i*diami)* \
        (diamr**(1.-dr)) /(cr*dr)*(diami**(1.-di)) /(ci*di)
        
    # calculate number of fragments for each of these collisions
    (nfrag,nfrag_freeze1,nfrag_freeze2,nfrag_liq)= \
        mode2_phillips_mass(f,mr,mi,ar,br,ai,bi,cr,dr,ci,di)
    
    # number of new ice from mode 2
    integrand1=integrand1*nfrag_freeze2

    # mass of rain collected by ice  - not sure equation 112 is correct way around?
#    integrand2=np.pi/4.*(diamr+diami)**2*np.abs(ar*diamr**br+ai*diami**bi)* \
#        cr*diamr**dr*n0r*np.exp(-lambda0r*diamr)*n0i*np.exp(-lambda0i*diami)
        
    #return (integrand1,integrand2)
    return integrand1

"""
calculates the parameters of distribution for given moments: N, q
"""
def calculate_parameters(N,q,c,d,alpha):
    # N=n0*gamma(1+alpha)/lambda0**(1+alpha)
    # q = c*n0*gamma(1+d+alpha)/lambda0**(1+d+alpha)
    # N/q=lambda**d*gamma(1+alpha)/(c*gamma(1+d+alpha))
    # lambda0=((c*gamma(1+d+alpha)/gamma(1+alpha))*N/q)**(1/d)
    lambda0=((c*np.math.gamma(1+d+alpha)/np.math.gamma(1+alpha))*N/q)**(1.0/d)
    n0=N*lambda0**(1+alpha)/np.math.gamma(1.0+alpha)
    return (n0,lambda0)

"""
calculates the fragments due to different modes - Phillips et al - check
"""
def mode2_phillips(f,mr,mi,ar,br,ai,bi,cr,dr,ci,di):
    if mr < mi:
        # diameters from mass-dimension relation (inverse)
        diamr=(mr/cr)**(1/dr)
        if diamr>0.15e-3:
            diami=(mi/ci)**(1/di)
    
            # fall speeds
            vr=ar*diamr**br
            vi=ai*diami**bi
            # CKE from equation 6
            K0=0.5*(mr*mi/(mr+mi))*(vr-vi)**2
            # DE parameter
            DE=K0/(gamma_liq*np.pi*diamr**2)
            # number of fragments in splash
            nfrag=3.*max(DE-DEcrit,0)
            # number of fragments in splash that freeze due to mode1
            nfrag_freeze1=nfrag*f
            # number of fragments in splash that freeze due to mode2
            nfrag_freeze2=nfrag*(1-f)*phi
        else:
            nfrag=0.
            nfrag_freeze1=0.
            nfrag_freeze2=0.
    else:
        nfrag=0.
        nfrag_freeze1=0.
        nfrag_freeze2=0.   

    # number of liquid drops left in splash
    nfrag_liq=nfrag-nfrag_freeze1-nfrag_freeze2

    return (nfrag,nfrag_freeze1,nfrag_freeze2,nfrag_liq)
    
    
"""
calculates the fragments due to different modes - Phillips et al - check
"""
def mode2_phillips_mass(f,mr,mi,ar,br,ai,bi,cr,dr,ci,di):
    # diameters from mass-dimension relation (inverse)
    diamr=(mr/cr)**(1/dr)
    diami=(mi/ci)**(1/di)

    # fall speeds
    vr=ar*diamr**br
    vi=ai*diami**bi
    # CKE from equation 6
    K0=0.5*(mr*mi/(mr+mi))*(vr-vi)**2
    # DE parameter
    DE=K0/(gamma_liq*np.pi*diamr**2)
    # number of fragments in splash
    nfrag=3.*max(DE-DEcrit,0)
    # number of fragments in splash that freeze due to mode1
    nfrag_freeze1=nfrag*f
    # number of fragments in splash that freeze due to mode2
    nfrag_freeze2=nfrag*(1-f)*phi

    # number of liquid drops left in splash
    nfrag_liq=nfrag-nfrag_freeze1-nfrag_freeze2

    return (nfrag,nfrag_freeze1,nfrag_freeze2,nfrag_liq)
    
    
Ns=np.logspace(-2,4,5)*1000.
qs=np.logspace(-6,-3,6)
Ni1=np.logspace(-1,2,6)*1000.
eval2=np.zeros((5,6))
for i in range(len(Ns)):
    for j in range(len(Ni1)):
        # rain water conc / m-3 
        Nr=1.0e3 
        # rain water content kg / m-3 
        qr=3.e-4
        Nr=Ns[i]
        #qr=qs[j]
        
        # mass-dimension relations for rain
        cr=np.pi/6*1000.
        dr=3.0
        # fallspeed-dimension relation for rain - table 17
        ar=836.
        br=0.8
        # parameters for rain
        (n0r,lambda0r)=calculate_parameters(Nr,qr,cr,dr,alphar)

 
        # ice water conc / m-3 
        Ni=0.01e3 
        # ice water content kg / m-3 
        qi=.5e-3
        Ni=Ni1[j]
        #qi=qs[j]
        # mass-dimension relations for ice - see table A, page 46
        ci=0.04
        di=2.0
        # fallspeed-dimension relation for ice, Ferrier94 - table 17
        ai=8.97
        bi=0.42
        # parameters for ice
        (n0i,lambda0i)=calculate_parameters(Ni,qi,ci,di,alphai)


        mrthresh=cr*150e-6**dr
        mrupper=np.max([cr*(10e-3)**dr,mrthresh])
        miupper=np.max([cr*(10e-3)**dr,10.e-3])
        eval1=integrate.dblquad(collisions_precip_mass,mrthresh,mrupper,\
            lambda x: x, lambda x: miupper,\
             args=(n0r,lambda0r,n0i,lambda0i,ar,br,ai,bi,cr,dr,ci,di,f))
             
        eval2[i][j]=eval1[0]
        print(str(eval2[i][j]) + ' ' + str(lambda0r)+ ' ' + str(lambda0i))

import matplotlib.pyplot as plt
from numpy import ma
from matplotlib import ticker, cm

plt.ion()
plt.figure()
#cs=plt.contourf(Ni1/1000.,Ns/1000.,eval2,locator=ticker.LogLocator(), cmap=cm.PuBu_r)
cs=plt.contourf(Ni1/1000.,Ns/1000.,eval2, cmap=cm.PuBu_r)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Nr (L^{-1})')
plt.xlabel('Ni (L^{-1})')
cbar = plt.colorbar(cs)    
