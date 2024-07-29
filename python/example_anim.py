import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from netCDF4 import Dataset 

fileNames=['/tmp/output0.nc','/tmp/output1.nc','/tmp/output3.nc','/tmp/output7.nc']
titles=['All SIP','No SIP','Only HM','All SIP 1000 cm$^{-3}$']

plt.ion()
fig=plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(4, 1)
ax1=plt.subplot(gs[0,0])
ax2=plt.subplot(gs[1,0])
ax3=plt.subplot(gs[2,0])
ax4=plt.subplot(gs[3,0])
axes=[ax1,ax2,ax3,ax4]
gs.tight_layout(fig)
i=0;ii=-1
nc = Dataset(fileNames[i])
y = nc['y'][:]/1000
z = nc['z'][:]/1000
th = np.squeeze(nc['q'][ii,0,:,:,14])/1e6
axes[0].clear()
quad1=axes[i].pcolormesh(y-0.050,z,th.transpose(),cmap='gist_ncar',shading='gouraud', \
	norm=matplotlib.colors.LogNorm())
quad1.set_clim((10,1000))
cb=fig.colorbar(quad1,ax=axes)
cb.set_label('CDNC (cm$^{-3}$)')


# cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
def draw(ii):
	for i in range(len(axes)):
		axes[i].clear()
			
		# open file
		nc = Dataset(fileNames[i])
		# extract needed variables
		y = nc['y'][:]/1000
		z = nc['z'][:]/1000
		th = np.squeeze(nc['q'][ii,0,:,:,14])/1e6
		time=nc['time'][:]
		
		quad1=axes[i].pcolormesh(y-0.050,z,th.transpose(),cmap='gist_ncar',shading='gouraud', \
			norm=matplotlib.colors.LogNorm())
# 		axes[i].set_title(titles[i] + ' and time: ' + str(time[ii]) + 's')
		time_text = axes[i].text(0.02, 0.9, \
			titles[i] + ' and time: ' + str(round(time[ii]/60.0,2)) + 'mins', transform=axes[i].transAxes)
		quad1.set_clim((10,1000))
# 		if((ii==0) and (i==0)):
# 			cb=fig.colorbar(quad1,ax=axes)
# 			cb.set_label('CDNC (cm$^{-3}$)')
		axes[i].set_aspect('equal')
# 		axes[i].autoscale(tight=True)
		axes[i].set_xlim((-7.500,7.500))
		axes[i].set_ylabel('z (km)')
		axes[i].set_xlabel('y (km)')
		
		th = np.squeeze(nc['precip'][ii,0,:,:,0])
		quad2=axes[i].contour(y-0.05,z,th.transpose(),colors='c',levels=[0.01,1,10])
		
		th = np.squeeze(nc['q'][ii,0,:,:,30])/1000.0
		quad3=axes[i].contour(y-0.05,z,th.transpose(),colors='m',levels=[0.1,1,10,50])
		nc.close()
		
	
def init():
    pass

def update_img(iter):
	draw(iter)

# gs.tight_layout(fig)

anim = animation.FuncAnimation(fig,update_img,init_func=init,frames=360,interval=50,blit=False,repeat=False)
plt.show()

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# anim.save('/tmp/animation.mp4',writer=writer)
print('Finished!!')

"""
	now plot out the precip
"""
fileNames=['/tmp/output0.nc','/tmp/output1.nc','/tmp/output2.nc','/tmp/output3.nc','/tmp/output4.nc','/tmp/output7.nc']
titles=['All SIP','No SIP','No HM','Only HM','No HM or M2','All SIP 1000 cm$^{-3}$']

plt.figure()
dat=[]
str1=['-','-','-','--','--','--']
ipeak=[90,132,132,90,132,180]
for i in range(len(fileNames)):
	# open file
	nc = Dataset(fileNames[i])
	# extract needed variables
	time=nc['time'][:]
	precip=np.mean(nc['precip'][:,0,:,0,0],axis=1)
	ice=np.max(np.max(nc['q'][:,0,:,:,30]/1000.0,axis=2),axis=1)
	
# 	plt.plot(time/60.0,np.cumsum(precip)*10.0/3600)
	dat1=np.cumsum(precip)*10.0/3600
	dat.append(dat1[-1])
	
# 	plt.legend(titles)
	ii=ipeak[i]
	rain=np.squeeze(np.max(nc['q'][ii,0,:,:,22]/1000.0,axis=0)).data[:]
	plt.subplot(222)
	plt.plot(rain[:40],nc['trefn'][:40]-273.15,ls=str1[i])
	plt.xlabel('Inst. Rain Conc. (L$^{-1}$)')
	plt.ylabel('T ($^\circ$C)')
	if (i==0):
		plt.gca().invert_yaxis()
	plt.legend(titles)
	
# 	plt.legend(titles)
	plt.subplot(223)
	plt.plot(time/60.0,ice,ls=str1[i])
	plt.ylim((-4.63426437377929,150))
	plt.ylabel('Maximum Ice Conc. (L$^{-1}$)')
	plt.xlabel('time (mins)')
	plt.legend(titles,loc=2)
	
# 	plt.legend(titles)
	ice=np.squeeze(np.max(nc['q'][ii,0,:,:,30]/1000.0,axis=0)).data[:]
	plt.subplot(224)
	plt.plot(ice[:40],nc['trefn'][:40]-273.15,ls=str1[i])
	plt.xlabel('Inst. Ice Conc. (L$^{-1}$)')
	plt.ylabel('T ($^\circ$C)')
	if (i==0):
		plt.gca().invert_yaxis()
	plt.legend(titles)
	
	nc.close()
plt.subplot(221)
plt.bar(titles,dat)
plt.ylabel('Accum. precipitation (mm)')
plt.savefig('/tmp/summary.png')
