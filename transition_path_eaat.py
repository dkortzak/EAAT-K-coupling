#!/usr/bin/env python
from msmbuilder.msm import MarkovStateModel
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import msmexplorer as msme
trajs=[]
for i in range(1,104):
    for chain in ["M","N","O"]:
        trajs.append(np.loadtxt("data/freeMD_EAAT1/run"+str(i)+"/states_"+chain+".xvg",comments=["@","#"])[2:,1])
msm = MarkovStateModel(lag_time=1,ergodic_cutoff="off",reversible_type=None)    
msm.fit(trajs)
if hasattr(msm, 'all_populations_'):
        tmat = msm.all_transmats_.mean(0)
        pop = msm.all_populations_.mean(0)
elif hasattr(msm, 'populations_'):
        tmat = msm.transmat_
        pop = msm.populations_
graph = nx.Graph(tmat)
pos = nx.spring_layout(graph)
statesm=[0,10,100,110,1000,1010,1100,1110,1,11,101,111,1001,1011,1101,1111]
statelabel=["Apo","K3","K2","K2K3","K1","K1K3","K1K2","K1\nK2K3","K4","K3K4","K2K4","K2\nK3K4","K1K4","K1\nK3K4","K1\nK2K4","K1K2\nK3K4"]
f=dict({})
print msm.mapping_.items()
for i in range(16):
    f[int(msm.mapping_.items()[i][1])]=statelabel[statesm.index(int(msm.mapping_.items()[i][0]))]
ax=plt.figure(figsize=(6,6))
msme.plot_msm_network(msm,  pos=pos,node_color='pomegranate',edge_color='gray',labels=f,width=.5)
msme.plot_tpaths(msm,pos=pos,sources=[0],sinks=[8],num_paths=2,edge_color='r',labels=f,arrows=True)
plt.show()