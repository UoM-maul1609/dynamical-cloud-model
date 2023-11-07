import matplotlib
matplotlib.use('Agg') # display not needed (if run over server)

from matplotlib import pyplot as plt

from netCDF4 import Dataset 
import numpy as np

# open file
nc = Dataset('/tmp/output.nc')

# extract needed variables
y = nc['y'][:]
z = nc['z'][:]
th = np.squeeze(nc['q'][1,0,:,:,30])

print(np.shape(th))

time=nc['time'][:]


# plot out
plt.pcolormesh(y,z,th.transpose(),cmap='RdBu_r',shading='gouraud')
plt.title('time: ' + str(time[-1]) + 's')

# plt.clim((-0.05,0.05))
plt.colorbar()

plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)


# print to a file
plt.savefig('/tmp/foo.png', bbox_inches='tight')
