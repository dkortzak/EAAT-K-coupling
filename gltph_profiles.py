#!/usr/bin/env python
from scipy import stats, integrate
import matplotlib.pyplot as plt
import numpy as np

def inverse_line(x,a,b):
    return (1/a)*x-b/a
def bootstrap(bsprofsx,bsprofs,param):
  pbs0=np.zeros(1000)
  Keqs1=np.zeros(1000)
  for i in range(1000):
    x=inverse_line(bsprofsx[i,:],*param)-9.
    N=np.argmin(x**2)
    sam=-bsprofs[i,:]/integrate.simps(bsprofs[i,:],x)   
    Keqs1[i]=(1./(sam[:N].sum()/sam[N:].sum()))
    pbs0[i]=-integrate.simps(sam[N:],x[N:])
  print "Mean, SD Pclose from bootstraping"
  print np.mean(pbs0),np.std(pbs0)
  print "95CI for Pclose from bootstraping"
  print np.sort(pbs0)[25],np.sort(pbs0)[975]
  
  
#linear fit of eigenvalue vs HP distance
param_gltph_wt_apo=[-0.30243795,  4.1552556 ]
param_gltph_wt_k=[-0.22719732,  3.61298254]


plt.figure(figsize=(3,3))
inbs1tt=np.loadtxt("data/umbrella_sampling_GLTPH/apo/bsProfs.xvg",comments=["@","#","&"])
bsprofs1tt=np.reshape(inbs1tt[:,1],newshape=(1000,200))
bsprofs1xtt=np.reshape(inbs1tt[:,0],newshape=(1000,200))
x=bsprofs1xtt[0,:]
y1=np.mean(bsprofs1tt,axis=0)
yerr=np.std(bsprofs1tt,axis=0)
scale=-integrate.simps(y1,inverse_line(x,*param_gltph_wt_apo))
y1=y1/scale
yerr=yerr/scale
bootstrap(bsprofs1xtt,bsprofs1tt,param_gltph_wt_k)
plt.fill_between(inverse_line(x,*param_gltph_wt_apo),y1-yerr,y1+yerr,color="k",alpha=0.25)
plt.errorbar(inverse_line(x,*param_gltph_wt_apo),y1,fmt="-",color="k")
inbs1tt=np.loadtxt("data/umbrella_sampling_GLTPH/K1bound/bsProfs.xvg",comments=["@","#","&"])
bsprofs1tt=np.reshape(inbs1tt[:,1],newshape=(1000,200))
bsprofs1xtt=np.reshape(inbs1tt[:,0],newshape=(1000,200))
y=np.mean(bsprofs1tt,axis=0)
yerr=np.std(bsprofs1tt,axis=0)
x=bsprofs1xtt[0,:]
scale=-integrate.simps(y,inverse_line(x,*param_gltph_wt_k))
y=y/scale
yerr=yerr/scale
plt.fill_between(inverse_line(x,*param_gltph_wt_k),y-yerr,y+yerr,color="cornflowerblue",alpha=0.25)
plt.errorbar(inverse_line(x,*param_gltph_wt_k),y,fmt="-",color="b")
plt.xlim(4,23)
plt.show()
bootstrap(bsprofs1xtt,bsprofs1tt,param_gltph_wt_k)