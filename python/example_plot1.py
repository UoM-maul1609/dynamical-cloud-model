import matplotlib
# matplotlib.use('Agg') # display not needed (if run over server)

from matplotlib import pyplot as plt

from netCDF4 import Dataset 
import numpy as np


plt.ion()
# open file
nc1 = Dataset('/tmp/output0.nc')
nc2 = Dataset('/tmp/output1.nc')

# extract needed variables
y = nc1['y'][:]
z = nc1['z'][:]
ii=130
itop=43
qnc1 = np.squeeze(nc1['q'][ii,0,:,:,14])/1e6
qni1 = np.squeeze(nc1['q'][ii,0,:,:,30])/1e3
th1 = np.squeeze(nc1['th'][ii,0,:,0:itop])
tref1 = nc1['trefn'][0:itop]
precip1 = np.mean(np.squeeze(nc1['precip'][:,0,:,0,0]),axis=1)
precip1a = np.squeeze(nc1['precip'][ii,0,:,:,0])

qnc2 = np.squeeze(nc2['q'][ii,0,:,:,14])/1e6
qni2 = np.squeeze(nc2['q'][ii,0,:,:,30])/1e3
th2 = np.squeeze(nc2['th'][ii,0,:,0:itop])
precip2 = np.mean(np.squeeze(nc2['precip'][:,0,:,0,0]),axis=1)
precip2a = np.squeeze(nc2['precip'][ii,0,:,:,0])
tref2 = nc2['trefn'][0:itop]

time=nc1['time'][:]

fig=plt.figure()
# plot out
ax1=plt.subplot(221)
plt.pcolormesh((y-50.0)/1000,z/1000,qni1.transpose(),\
	shading='gouraud', norm=matplotlib.colors.LogNorm())
ax1.set_xlim((-7.5,7.5))
plt.title('Simulation without mode-2')
cb1=plt.colorbar()
cb1.set_label('ice conc (L$^{-1}$)')
plt.clim((0.01,50))
# plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)
plt.contour((y-50.0)/1000,z/1000,qnc1.transpose(),colors='c',levels=[0.01,50,100])
# plt.contour((y-50.0)/1000,z/1000,precip1a.transpose(),colors='r',levels=[0.05,10])

CS1=plt.contour((y-50.0)/1000,z[0:itop]/1000,(th1+tref1-273.15).transpose(),colors='m',levels=[-7.,0])
ax1.clabel(CS1, CS1.levels, inline=True,fontsize=10)
plt.ylabel('z (km)')
plt.xlabel('x (km)')
plt.text(0.1,0.9,'(a) no mode-2',transform=ax1.transAxes)


# plot out
ax2=plt.subplot(222)
plt.pcolormesh((y-50.0)/1000,z/1000,qni2.transpose(),\
	shading='gouraud', norm=matplotlib.colors.LogNorm())
ax2.set_xlim((-7.5,7.5))
plt.title('Simulation with mode-2')
cb2=plt.colorbar()
cb2.set_label('ice conc (L$^{-1}$)')

plt.clim((0.01,50))
# plt.clim((0,40))
# plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)


plt.contour((y-50.0)/1000,z/1000,qnc2.transpose(),colors='c',levels=[0.01,50,100])
# plt.contour((y-50.0)/1000,z/1000,precip2a.transpose(),colors='r',levels=[0.05,10])

CS2=plt.contour((y-50.0)/1000,z[0:itop]/1000,(th2+tref2-273.15).transpose(),colors='m',levels=[-7.,0])

ax2.clabel(CS2, CS2.levels, inline=True, fontsize=10)
plt.ylabel('z (km)')
plt.xlabel('x (km)')
plt.text(0.1,0.9,'(b) mode-2',transform=ax2.transAxes)



ax1.set_xlim((-6,6))
ax2.set_xlim((-6,6))

ax3=plt.subplot(212)
plt.plot(time/60,precip1)
plt.plot(time/60,precip2)
plt.yscale('log')
plt.ylim((0.001,5))
plt.ylabel('domain av. precip.\nrate (mm hr$^{-1}$)')
plt.xlabel('time (mins)')
plt.legend(['without mode-2','with mode-2'],loc=3)
plt.text(0.1,0.9,'(c) Precipitation rate',transform=ax3.transAxes)

# print to a file
nc1.close()
nc2.close()
fig.tight_layout()

plt.savefig('/tmp/foo.png', bbox_inches='tight')
