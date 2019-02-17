import matplotlib
matplotlib.use('Agg') # display not needed (if run over server)

from matplotlib import pyplot as plt

from netCDF4 import Dataset 

# open file
nc = Dataset('/tmp/output.nc')

# extract needed variables
y = nc['y'][:]
z = nc['z'][:]
th = nc['th'][-1,:,:,:]


# plot out
plt.pcolor(y,z,th[0].transpose())
plt.clim((-0.05,0.05))

plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)


# print to a file
plt.savefig('/tmp/foo.png', bbox_inches='tight')