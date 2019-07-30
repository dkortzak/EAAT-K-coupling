#!/usr/bin/env python
from scipy import stats, integrate
import matplotlib.pyplot as plt
import numpy as np

def line(x,a,b):
    return a*x+b
  
def bootstrap(bsprofsx,bsprofs,param):
  pbs0=np.zeros(1000)
  Keqs1=np.zeros(1000)
  b1=np.zeros(1000)
  b2=np.zeros(1000)
  for i in range(1000):
    x=line(bsprofsx[i,:],*param)-9.
    N=np.argmin(x**2)
    sam=-bsprofs[i,:]/integrate.simps(bsprofs[i,:],x)
    Keqs1[i]=(1./(sam[:N].sum()/sam[N:].sum()))
    pbs0[i]=-integrate.simps(sam[N:],x[N:])
  print "Mean, SD Pclose from bootstraping"
  print np.mean(pbs0),np.std(pbs0)
  print "95CI for Pclose from bootstraping"
  print np.sort(pbs0)[25],np.sort(pbs0)[975]
  return np.mean(pbs0),np.std(pbs0)

# linear fit of eigenvector projection vs HP distances
param_eaat1_wt_apo=[-2.67204882,  9.44401941]
param_eaat1_wt_k=[-3.06904798,  9.35961961]
param_eaat1_la_apo=[-2.65044501,  9.40658272]
param_eaat1_la_k=[ -2.9759215 ,  10.10027717]
param_eaat1_lt_apo=[ -2.79933929,  10.42532599]
param_eaat1_lt_k=[-2.98847417,  9.61308078]
param_eaat1_ra_apo=[-2.70166366,  9.2092695]
param_eaat1_ra_k=[-3.04853192,  9.65655686]


def plot_profile(name,fname1,fname2,params1,params2):
  print name
  plt.figure(figsize=(3,3))
  inbs1tt=np.loadtxt(fname1,comments=["@","#","&"])
  bsprofs1tt=np.reshape(inbs1tt[:,1],newshape=(1000,200))
  bsprofs1xtt=np.reshape(inbs1tt[:,0],newshape=(1000,200))
  y=np.mean(bsprofs1tt,axis=0)
  x=bsprofs1xtt[0,:]
  yerr=np.std(bsprofs1tt,axis=0)
  scale=-integrate.simps(y,line(x,*params1))
  y=y/scale
  yerr=yerr/scale
  plt.fill_between(line(x,*params1),y-yerr,y+yerr,color="dimgray",alpha=0.25)
  plt.errorbar(line(x,*params1),y,color="dimgray")
  m_wt_apo,sd_wt_apo=bootstrap(bsprofs1xtt,bsprofs1tt,params1)
  inbs1tt=np.loadtxt(fname2,comments=["@","#","&"])
  bsprofs1tt=np.reshape(inbs1tt[:,1],newshape=(1000,200))
  bsprofs1xtt=np.reshape(inbs1tt[:,0],newshape=(1000,200))
  y=np.mean(bsprofs1tt,axis=0)
  yerr=np.std(bsprofs1tt,axis=0)
  x=bsprofs1xtt[0,:]
  scale=-integrate.simps(y,line(x,*params2))
  y=y/scale
  yerr=yerr/scale
  plt.fill_between(line(x,*params2),y-yerr,y+yerr,color="blue",alpha=0.25)
  plt.errorbar(line(x,*params2),y,color="blue")
  plt.xlim(4,23)
  plt.ylim(0,0.35)
  plt.plot([9,9],[0.,0.2],color="r")
  plt.show()
  m_wt_k,sd_wt_k=bootstrap(bsprofs1xtt,bsprofs1tt,params2)


plot_profile("EAAT1 WT","data/umbrella_sampling_EAAT1/WT/apo/bsProfs.xvg","data/umbrella_sampling_EAAT1/WT/K1bound/bsProfs.xvg",param_eaat1_wt_apo,param_eaat1_wt_k)
plot_profile("EAAT1 R479A","data/umbrella_sampling_EAAT1/R479A/apo/bsProfs.xvg","data/umbrella_sampling_EAAT1/R479A/K1bound/bsProfs.xvg",param_eaat1_ra_apo,param_eaat1_ra_k)
plot_profile("EAAT1 L448T","data/umbrella_sampling_EAAT1/L448T/apo/bsProfs.xvg","data/umbrella_sampling_EAAT1/L448T/K1bound/bsProfs.xvg",param_eaat1_lt_apo,param_eaat1_lt_k)
plot_profile("EAAT1 L448A","data/umbrella_sampling_EAAT1/L448A/apo/bsProfs.xvg","data/umbrella_sampling_EAAT1/L448A/K1bound/bsProfs.xvg",param_eaat1_la_apo,param_eaat1_la_k)
