import matplotlib
# matplotlib.use('Agg') # display not needed (if run over server)

from matplotlib import pyplot as plt

from netCDF4 import Dataset 
import numpy as np

# open file
nc = Dataset('/tmp/output7.nc')

# extract needed variables
y = nc['y'][:]
print(y-50.0)
z = nc['z'][:]
ii=-1
th = np.squeeze(nc['q'][ii,0,:,:,14])/1e6
# th = np.squeeze(nc['q'][ii,0,:,:,15])+np.squeeze(nc['q'][ii,0,:,:,23])
# th = np.squeeze(nc['q'][ii,0,:,:,15]) /1e6
# th = np.squeeze(nc['q'][ii,0,:,:,0]) 
# th = np.squeeze(nc['q'][ii,0,:,:,15])*1000.
# th = np.squeeze(nc['q'][ii,0,:,:,30]) /1000.
	#+np.squeeze(nc['q'][ii,0,:,:,23])+np.squeeze(nc['q'][ii,0,:,:,31])

print(np.shape(th))

time=nc['time'][:]


# plot out
# plt.subplot(211)
plt.pcolormesh(y-50.0,z,th.transpose(),cmap='RdBu_r',shading='gouraud')
plt.title('time: ' + str(time[ii]) + 's')
# plt.clim((0,0.1))
# plt.clim((0,1e-5))
plt.colorbar()
# plt.clim((0,10))
cl=plt.clim()

plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)


th = np.squeeze(nc['precip'][ii,0,:,:,0])
plt.contour(y-50.0,z,th.transpose(),colors='c',levels=[0.01,1,10])

th = np.squeeze(nc['q'][ii,0,:,:,30])/1000.0
plt.contour(y-50.0,z,th.transpose(),colors='m',levels=[0.1,1,10,50])

plt.clim(cl)
plt.xlim((-10000,10000))
# plt.subplot(212)
# print(plt.xlim())
# plt.plot(y-50.0,th[:,0])
# plt.grid()

# print to a file
plt.savefig('/tmp/foo.png', bbox_inches='tight')
nc.close()
