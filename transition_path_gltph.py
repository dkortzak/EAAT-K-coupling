#!/usr/bin/env python
from msmbuilder.msm import MarkovStateModel
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import msmexplorer as msme
trajs=[]
for i in range(1,76):
    for chain in ["M","N","O"]:
        trajs.append(np.loadtxt("data/freeMD_GLTPH_IFC/run"+str(i)+"/states_"+chain+".xvg",comments=["@","#"])[:,1])
msm = MarkovStateModel(lag_time=1)   
msm.fit(trajs)
if hasattr(msm, 'all_populations_'):
        tmat = msm.all_transmats_.mean(0)
        pop = msm.all_populations_.mean(0)
elif hasattr(msm, 'populations_'):
        tmat = msm.transmat_
        pop = msm.populations_
graph = nx.Graph(tmat)
pos = nx.spring_layout(graph)
f=dict({})
statesm=[0,1,10,11,100,101,110,111]
statelabel=["Apo","K3","K2","K2K3","K1","K1K3","K1K2","K1\nK2K3"]
for i in range(len(msm.mapping_.items())):
    f[int(msm.mapping_.items()[i][1])]=statelabel[statesm.index(int(msm.mapping_.items()[i][0]))]
ax=plt.figure(figsize=(6,3))
msme.plot_msm_network(msm,  pos=pos,node_color='pomegranate',edge_color='gray',labels=f,width=.5)
msme.plot_tpaths(msm,pos=pos,sources=[0],sinks=[4],num_paths=2,edge_color='r',labels=f,arrows=True)
plt.show()

trajs=[]
for i in range(1,80):
    for chain in ["M","N","O"]:
        trajs.append(np.loadtxt("data/freeMD_GLTPH_OFC/run"+str(i)+"/states_"+chain+".xvg",comments=["@","#"])[:,1])
msm = MarkovStateModel(lag_time=1,ergodic_cutoff="off")    
msm.fit(trajs)
if hasattr(msm, 'all_populations_'):
        tmat = msm.all_transmats_.mean(0)
        pop = msm.all_populations_.mean(0)
elif hasattr(msm, 'populations_'):
        tmat = msm.transmat_
        pop = msm.populations_
graph = nx.Graph(tmat)
pos = nx.spring_layout(graph)
f=dict({})
statesm=[0,1,10,11,100,101,110,111]
statelabel=["Apo","K3","K2","K2K3","K1","K1K3","K1K2","K1\nK2K3"]
for i in range(len(msm.mapping_.items())):
    f[int(msm.mapping_.items()[i][1])]=statelabel[statesm.index(int(msm.mapping_.items()[i][0]))]    
ax=plt.figure(figsize=(6,3))
msme.plot_msm_network(msm,  pos=pos,node_color='pomegranate',edge_color='gray',labels=f,width=.5)
msme.plot_tpaths(msm,pos=pos,sources=[0],sinks=[4],num_paths=2,edge_color='r',labels=f,arrows=False)
plt.show()
