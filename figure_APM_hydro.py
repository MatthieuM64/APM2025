#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import time
import multiprocessing

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
#matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

epsilon=2.4
rho0=3
beta=0.75
LX=200
LY=200
Rd=10
rhod=15
init=0
tmax=10000
NCPU=4
multi=True
movie=False

for arg in sys.argv[1:]:
	if "-epsilon=" in arg:
		epsilon=float(arg[9:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-beta=" in arg:
		beta=float(arg[6:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-Rd=" in arg:
		Rd=float(arg[4:])
	elif "-rhod=" in arg:
		rhod=float(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

def delta_snap(t):
	if t<1000:
		return 2
	else:
		return 20

dpi=160

if init==0:
	colors_map=[(0,0,1),(1,1,1),(1,0,0)]
	cmap=matplotlib.colors.LinearSegmentedColormap.from_list('my_list', colors_map, N=256)
else:
	colors_map=[(0,0,1),(1,1,1),(0,1,0)]
	cmap=matplotlib.colors.LinearSegmentedColormap.from_list('my_list', colors_map, N=256)

def Snapshot(i):
	t=TIME[i]
	path='snapshots/figure_APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_%d.txt'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,i)
	if not os.path.isfile(path):
		fig=plt.figure(figsize=(8,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.07,right=0.97,bottom=0.05, top=0.95,hspace=0.25,wspace=0.25)
		
		ax=plt.subplot(gs[0,0])
		
		MAG=np.loadtxt('data_APM4_dynamics2d/APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_t=%d.txt'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,t))
		N=len(MAG)
		x=np.linspace(0,LX,N+1)
		y=np.linspace(0,LY,N+1)
		
		### Plot the densities with different colormaps (when is the highest state)
		plt.pcolormesh(LX-x,y,MAG,rasterized=True,vmin=-1,vmax=1,cmap=cmap)
		cb=plt.colorbar(ticks=[-1,0,1])
		cb.solids.set_rasterized(True)
		if init==0:
			cb.ax.set_ylabel('$m^{\sigma=1}/\\rho$',rotation=90,labelpad=-10)
		elif init==1:
			cb.ax.set_ylabel('$m^{\sigma=2}/\\rho$',rotation=90,labelpad=-10)
		
		plt.axis('equal')	
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))

		plt.text(0.01*LX,0.04*LY,'$t=%d$'%(t),ha='left',va='center',fontsize=17)
		plt.text(LX,1.04*LY,'$\\beta=%.8g$, $\epsilon=%.8g$, $\\rho_0=%.8g$'%(beta,epsilon,rho0),ha='right',va='center',fontsize=17)
		plt.text(0.99*LX,0.95*LY,'$r_d=%.8g$, $\\rho_0^d=%.8g$'%(Rd,rhod),ha='right',va='center',fontsize=17)

		plt.savefig(path,dpi=dpi)
		plt.close()		
		
		print('-snap=%d/%d -t=%.8g -tcpu=%d'%(i+1,Nsnap,t,time.time()-clock))
		del MAG,fig


os.system('mkdir -p snapshots')

i=0
t=0
ARG=[]
TIME=[]
while os.path.isfile('data_APM4_dynamics2d/APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_t=%d.txt'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,t)) and t<=tmax:
	ARG.append(i)
	TIME.append(t)
	i+=1
	t+=delta_snap(t)

Nsnap=len(ARG)	
print('%d Snapshots'%Nsnap)
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
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_%%01d.png -c:v h264 -r 25 -crf 25 -s %dx%d movies/movie_APM4_hydro_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d.mp4'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,8*dpi,6*dpi,beta,epsilon,rho0,LX,LY,Rd,rhod,init))

print('OK - time=%d sec'%(time.time()-clock))
