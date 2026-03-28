#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import sys
import time
import multiprocessing
import gc

fontsize=25
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

beta=1.
D=0.3
epsilon=2.5
rho0=5.
LX=512
LY=512
init=3
ran=0
rhomax=30.
tmax=600000
NCPU=4
multi=True
movie=False


for arg in sys.argv[1:]:
	if "-beta=" in arg:
		beta=float(arg[6:])
	elif "-D=" in arg:
		D=float(arg[3:])
	elif "-epsilon=" in arg:
		epsilon=float(arg[9:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-gamma=" in arg:
		gamma=float(arg[7:])
	elif "-rhomax=" in arg:
		rhomax=float(arg[8:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU>1:
	multi=True
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

DT=150
DeltaT=1./(4*D+np.exp(4*beta))
rhomid=0.5*rhomax

colors_map1=[(1,0,0),(0,1,0),(0,0,1),(0,0,0)]
cmap1=LinearSegmentedColormap.from_list('my_list1', colors_map1, N=256)

colors_map2=[(0,0,0),(1,0.5,0),(1,1,0),(1,1,1)]
cmap2=matplotlib.colors.LinearSegmentedColormap.from_list('my_list2', colors_map2, N=256)

def Snapshot(i):
	if not os.path.isfile('snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,D,epsilon,rho0,LX,LY,init,ran,i)):
		t=i*DT
		
		
		x=np.linspace(0,LX,LX)
		y=np.linspace(0,LY,LY)

		fig=plt.figure(figsize=(14,6))
		gs=matplotlib.gridspec.GridSpec(1,2,width_ratios=[1,1],height_ratios=[1],left=0.06,right=1.0,bottom=0.08,top=0.92,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		
		STATE=np.loadtxt('data_APM_dynamics2d/APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,LX,LY,init,ran,i*DT))
		#RHO=np.loadtxt('data_APM_dynamics2d/APM_density_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,LX,LY,init,ran,i*DT))
		STATE=np.ma.masked_where(STATE!=STATE,STATE)
		
		#cmap=plt.get_cmap('bwr')
		plt.pcolormesh(x,y,STATE+1,vmin=1,vmax=4,rasterized=True,cmap=cmap1)
		cb=plt.colorbar(ticks=[1,2,3,4])
		cb.solids.set_rasterized(True)
		cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.03*LY,'$t=%d$'%(i*DT*DeltaT),ha='left',va='center',fontsize=20)
		plt.text(LX,1.03*LY,'$\\beta=%.8g$, $D=%.8g$, $\epsilon=%.8g$, $\\rho_0=%.8g$'%(beta,D,epsilon,rho0),ha='right',va='center',fontsize=20)
		
		
		ax=plt.subplot(gs[0,1])
		
		RHO=np.loadtxt('data_APM_dynamics2d/APM_density_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,LX,LY,init,ran,i*DT))
		
		#cmap=plt.get_cmap('Greys')
		plt.pcolormesh(x,y,RHO,vmin=0,vmax=rhomax,rasterized=True,cmap=cmap2)
		cb=plt.colorbar(ticks=[0,rhomax/3.,2*rhomax/3.,rhomax])
		cb.solids.set_rasterized(True)
		cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.03*LY,'$t=%d$'%(i*DT*DeltaT),ha='left',va='center',fontsize=20)
		plt.text(LX,1.03*LY,'$\\beta=%.8g$, $D=%.8g$, $\epsilon=%.8g$, $\\rho_0=%.8g$'%(beta,D,epsilon,rho0),ha='right',va='center',fontsize=20)
		
		plt.savefig('snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,D,epsilon,rho0,LX,LY,init,ran,i),dpi=120)
		plt.close()
		
		print('-snap=%d/%d -t=%d -tcpu=%d'%(i+1,Nsnap,i*DT,time.time()-clock))
		del fig,RHO,STATE


os.system('mkdir -p snapshots/')

i=0
ARG=[]
while os.path.isfile('data_APM_dynamics2d/APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,LX,LY,init,ran,i*DT)) and i*DT<=tmax:
	ARG.append(i)
	i+=1
	
Nsnap=len(ARG)
	
print('%d Snapshots'%len(ARG))
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

if movie:
	os.system('mkdir -p movies')	
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s 1680x720 movies/movie_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(beta,D,epsilon,rho0,LX,LY,init,ran,beta,D,epsilon,rho0,LX,LY,init,ran))
	if Nsnap>tmax/DT:
		os.system('rm snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_*.png'%(beta,D,epsilon,rho0,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
