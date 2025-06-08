from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from matplotlib import cm
from matplotlib.ticker import StrMethodFormatter
import matplotlib.colors as col
import matplotlib.pyplot as plt
import numpy as np
import sys,json
import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font',family='Arial') ## might be commented out if now work

title='dbh-unsub-pop' # title of this figure

#Modified by Leticia Gomes on Oct 27,2022
#Modification: read .dat generated from trajectory3 control file
data = np.loadtxt("average-dbh-unsub.dat", skiprows=1, dtype=float) 

#Format .dat file timestep, Ekin, Epot0-3, Etot0-3,Pop0-3
#Population info will start at index  10
p1=data[:,10]
p2=data[:,11]
p3=data[:,12]
p4=data[:,13]


### Figure 1
fig=plt.figure()
ax = fig.add_subplot()
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(True)
plt.subplots_adjust(wspace=0.3,hspace=0.3,bottom=0.2)
ax.axes.tick_params(axis='both',direction='in')
#### Format
xlabel=np.arange(0.1,2.01,0.3)
ylabel=np.arange(0,1.01,0.2)
#fig.suptitle('%s' % (title),fontsize=18)

ax.set_xlabel(r'Simulation timestep',fontsize=14,labelpad=8)
ax.set_ylabel(r'State population',fontsize=14,labelpad=8)
ax.set_xlim(0,2)
ax.set_ylim(0,1)
ax.set_xticks(xlabel)
ax.set_xticklabels(xlabel,fontsize=14,rotation=0)
ax.set_yticks(ylabel)
ax.set_yticklabels(ylabel,fontsize=14)

ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))

ax.set_aspect(1/ax.get_data_ratio())

#### Plot data

x=np.arange(0,2,0.001)
y1=np.mean(p1,axis=0)
y2=np.mean(p2,axis=0)
y3=np.mean(p3,axis=0)
y4=np.mean(p4,axis=0)
#y5=np.mean(p5,axis=0)
#y6=np.mean(p6,axis=0)
#y7=np.mean(p7,axis=0)
#y8=np.mean(p8,axis=0)

#cutoff=2000
ax.plot(x,p1,marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='black')
ax.plot(x,p2,marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='red')
#ax.plot(x,p3,marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='blue')

#ax.plot(x[0:cutoff],y2[0:cutoff],marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='red')
#ax.plot(x[0:cutoff],y3[0:cutoff],marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='blue')

#ax.plot(x[0:cutoff],y4[0:cutoff],marker='o',markersize=0.,linewidth=2.0,alpha=1.0,color='green')

#Print the half lives for the different states
for n, p in enumerate(p1): #print the half life for the ground state
    if p>= 0.5:
        print(n) 
        break
for n, p in enumerate(p2): #print the half life for the first excited state
    if p<= 0.5:
        print(n)
        break
#plt.show()
fig.savefig("dbh-unsub-pop.png" ,bbox_inches="tight",dpi=400)
s0=np.max(y1)
s1=np.max(y2)
s2=np.max(y3)
s3=np.max(y4)


