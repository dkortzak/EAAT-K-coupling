#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.distances as mdad
def get_distances3():
    sitepos=np.asarray([K1site.center_of_mass(),K2site.center_of_mass(),K3site.center_of_mass(),K4site.center_of_mass()],dtype="float32")

    distances1 = mdad.distance_array(reference=Kions.positions, configuration=K1site.positions)
    distances2 = mdad.distance_array(reference=Kions.positions, configuration=K2site.positions)
    distances3 = mdad.distance_array(reference=Kions.positions, configuration=K3site.positions)
    distances4 = mdad.distance_array(reference=Kions.positions, configuration=K4site.positions)
    
    distances1=np.sqrt(.25*np.sum(distances1**2,axis=1))
    distances2=np.sqrt(1./3.*np.sum(distances2**2,axis=1))
    distances3=np.sqrt(.25*np.sum(distances3**2,axis=1))
    distances4=np.sqrt(1./3.*np.sum(distances4**2,axis=1))
    distances=np.transpose([distances1,distances2,distances3,distances4])

    return np.array(distances)
  
def assign_states3(laststate,distances):

        current_state=np.zeros(959)
        current_state_tmp=np.zeros(959)

        clearly_bound=np.zeros(shape=(4,959),dtype=bool)
        clearly_unbound=np.zeros(shape=(4,959),dtype=bool)
        clearly_bound[0,:] = distances[0,:]<3.55
        clearly_bound[1,:] = distances[1,:]<3.3
        clearly_bound[2,:] = distances[2,:]<3.55
        clearly_bound[3,:] = distances[3,:]<3.55
        
        clearly_unbound[0,:] = distances[0,:]>11.4
        clearly_unbound[1,:] = distances[1,:]>5.0
        clearly_unbound[2,:] = distances[2,:]>7.0
        clearly_unbound[3,:] = distances[3,:]>4.5
        ambigously_bound = np.invert(clearly_bound) * np.invert(clearly_unbound)
        current_state[clearly_unbound.all(axis=0)] = 0
        (sitescb, Kcb) = np.where(clearly_bound)
        kcbx=[]
        current_state[Kcb] = sitescb+1
        (sitesuab, Kuab) = np.where(ambigously_bound)
            
        if ambigously_bound.any():
                distances[np.invert(ambigously_bound)]=np.inf
                current_state[ambigously_bound.any(axis=0)] = np.argmin(distances[:,ambigously_bound.any(axis=0)],axis=0)+1

        zzz=Kuab
        zzz1=np.where(laststate==0)
        zzz3=np.setdiff1d(zzz,zzz1)
        zzz4=np.intersect1d(zzz,zzz1)
        current_state[zzz4]=0
        current_state[zzz3]=laststate[zzz3]
        current_state[Kcb] = sitescb+1
        return current_state
      
def full2red(state):
    r=""
    for i in range(1,5):
        if i in state:
            r=r+"1"
        else:
            r=r+"0"
    return str(r)   

N=0
#loop over simulations
for j in range(1,2):
  #loop over chains
  for chain in ["M","N","O"]:
    
    traj="data/freeMD_EAAT1/run"+str(j)+"/trj_fit_M"+str(k)+".xtc"
    pdb="data/freeMD_EAAT1/chM.pdb"
    uni = mda.Universe(pdb,traj) 
    out=open("data/freeMD_EAAT1/run"+str(j)+"/states_"+chain+".xvg","w")
    
    #select ligating oxygens for each site
    K1sel="segid M and ((resid 374 and name O) or \
    (resid 378 and name O) or (resid 463 and name O) or (resid 467 and name OD1 ))"
    K2sel="segid M and ((resid 376 and name O) or \
    (resid 417 and name O) or (resid 420 and name O))"
    K3sel="segid M and ((resid 343 and name O) or \
    (resid 456 and (name OD1 or name OD2)) or (resid 460 and name OG1))"
    K4sel="segid M and ((resid 386 and (name OE1 or name OE2)) or (resid 416 and name O))"
    Kions=uni.selectAtoms("name K")
    K1site=uni.select_atoms(K1sel)
    K2site=uni.select_atoms(K2sel)
    K3site=uni.select_atoms(K3sel)
    K4site=uni.select_atoms(K4sel)
    
    N=len(Kions) 
    current_state1=np.zeros(N)
    #loop over frames in trajectory
    for frame in uni.trajectory[0:]:
        #calculate mean squared distance for each K+ and assign states accordingly
        current_state1=assign_states2(current_state1,np.transpose(get_distances3()))
        out.write(str(50*frame)+"\t"+full2red(current_state1[np.where(current_state1!=0)])+"\n")
    out.close()
