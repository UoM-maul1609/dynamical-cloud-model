import numpy as np

fs=0.2 # fraction of bare soil
fv=0.8 # fraction of vegetated soil

Fs=10. # net downward minus upward solar irradiance at the top of canopy plus over bare soil (W /m2)
Fi=2.  #  downward thermal ir
fsFs=fs*Fs # net downward solar irradiance over bare soil

epsilon=0.95
sigma_boltz=5.67e-23
Tg=280.

Fng=fsFs+Fi-epsilon*sigma_boltz*Tg**4


# Equation 8.103
ustar=2. 
uaf=0.83*ustar
cf=0.01-0.003/uaf
Rf=1./(cf*uaf) # aerodynamic resistance in the foliage (sm-1)


# Equation 8.114 reduction in transpiration due to drying up of the soil towards wilting point
w_wilt=0.114
wgs=0.435
wg_avg=0.2 # average liquid water content in the root layers of soil
wcr=0.75*wgs
Fc=np.maximum(np.minimum((wg_avg-w_wilt)/(wcr-w_wilt),1.0),1.0e-12)


# Leaf stomata resistance
Tf=280.
Rmin=100. # page 749
Rst=Rmin/Fc*(1.0+(200./(Fs+0.1))**2)*400./((Tf-273.15)*(313.15-Tf))


print(Rst)