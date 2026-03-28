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

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

beta=1.
D=1.
epsilon=2.7
rho0=10.
Rd=10
rhod=12
LX=200
LY=200
init=0
ran=0
rhomax=3.
t0=6000
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
	elif "-Rd=" in arg:
		Rd=int(arg[4:])
	elif "-rhod=" in arg:
		rhod=float(arg[6:])
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
	elif "-ran=" in arg:
		ran=int(arg[5:])
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
	
def delta_snap(t):
	if t<100000-1:
		return 100
	else:
		return 1000

DeltaT=1./(4*D+np.exp(4*beta))
rhomid=0.5*rhomax

colors=[(1,0,0),(0,1,0),(0,0,1),(0,0,0)]
cmap_name='my_list'
cmap=LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

def Snapshot(i):
	if not os.path.isfile('snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran,i)):
		t=TIME[i]
		STATE=np.loadtxt('data_APM_dynamics2d/APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran,t))
		STATE=np.ma.masked_where(STATE!=STATE,STATE)
		
		x=np.linspace(0,LX,LX)
		y=np.linspace(0,LY,LY)

		fig=plt.figure(figsize=(7,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.1,right=1.04,bottom=0.06,top=0.94,hspace=0.1,wspace=0.1)
		
		ax=plt.subplot(gs[0,0])
		
		
		
		#cmap=plt.get_cmap('bwr')
		plt.pcolormesh(x,y,STATE,vmin=0,vmax=3,rasterized=True,cmap=cmap)
		cb=plt.colorbar(ticks=[0,1,2,3])
		cb.solids.set_rasterized(True)
		cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.03*LY,'$t=%d$'%((t-t0)*DeltaT),ha='left',va='center',fontsize=15)
		plt.text(LX,1.03*LY,'$\\beta=%.8g$, $D=%.8g$, $\\varepsilon=%.8g$, $\\rho_0=%.8g$, $R_d=%d$, $\\rho_d=%.8g$'%(beta,D,epsilon,rho0,Rd,rhod),ha='right',va='center',fontsize=15)
		
		plt.savefig('snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran,i),dpi=180)
		plt.close()
		
		print('-snap=%d/%d -t=%d -tcpu=%d'%(i+1,Nsnap,t-t0,time.time()-clock))
		del fig,RHO,STATE


os.system('mkdir -p snapshots/')

i=0
t=t0
ARG=[]
TIME=[]
while os.path.isfile('data_APM_dynamics2d/APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran,t)) and t<=tmax:
	ARG.append(i)
	TIME.append(t)
	i+=1
	t+=delta_snap(t)
	
Nsnap=len(ARG)
	
print('%d Snapshots'%len(ARG))
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG[::-1])
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

if movie:
	os.system('mkdir -p movies')	
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s 1260x1080 movies/movie_APM_state_beta=%.8g_D=%.8g_epsilon=%.8g_rho0=%.8g_Rd=%d_rhod=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran,beta,D,epsilon,rho0,Rd,rhod,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
