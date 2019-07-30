#!/usr/bin/env python
import numpy as np
from random import randint

anadir="data/freeMD_GLTPH_OFC"
chains=["M","N","O"]
def dwell_GLTPHOFC():
    states=["000","001","010","011","100","101","110","111"]
    transitions=np.zeros(shape=(len(states),len(states)),dtype=int)
    total_times=np.zeros(len(states))
    dwelltimes=[]
    for i in states:
        dwelltimes.append([])
    for i in range(len(states)):
        for j in states:
            dwelltimes[i].append([])
    dwelltime=0
    rates=np.zeros(shape=(len(states),len(states)))
    Keq=np.zeros(shape=(len(states),len(states)))
    for chain in chains:
     for run in range(1,80):
        in1=tuple(open(anadir+"/run"+str(run)+"/states_"+chain+".xvg"))
        for i in range(0,len(in1)-1):
            old=in1[i].split()[1]
            new=in1[i+1].split()[1]
            transitions[states.index(old),states.index(new)]=int(transitions[states.index(old),states.index(new)])+1
            total_times[states.index(old)]=total_times[states.index(old)]+50
            dwelltime=dwelltime+50
    for i in states:
        for j in states:
            rates[states.index(i),states.index(j)]=float(transitions[states.index(i),states.index(j)])/(total_times[states.index(i)]*1e-12)
    for i in states:
        for j in states:
            if transitions[states.index(i),states.index(j)]>1 and transitions[states.index(j),states.index(i)]>1:
                Keq[states.index(i),states.index(j)]=rates[states.index(j),states.index(i)]/rates[states.index(i),states.index(j)]
    return transitions,total_times*1e-12,rates,Keq

def dwell_GLTPHOFC_bs():
    Nbs=1000
    states=["000","001","010","011","100","101","110","111"]
    bschains=["M","N","O"]
    bsruns=range(1,81)
    bsres=[]
    for bss in range(0,Nbs):
        transitions=np.zeros(shape=(len(states),len(states)),dtype=int)
        total_times=np.zeros(len(states))
        rates=np.zeros(shape=(len(states),len(states)))
        Keq=np.zeros(shape=(len(states),len(states)))
        for count in range(240):
            k1=randint(0,2)
            k2=randint(1,79)

            in1=tuple(open(anadir+"/run"+str(k2)+"/states_"+bschains[k1]+".xvg"))
            for i in range(0,len(in1)-1):
                old=in1[i].split()[1]
                new=in1[i+1].split()[1]
                transitions[states.index(old),states.index(new)]=int(transitions[states.index(old),states.index(new)])+1
                total_times[states.index(old)]=total_times[states.index(old)]+50
        for i in states:
            for j in states:
                rates[states.index(i),states.index(j)]=float(transitions[states.index(i),states.index(j)])/(total_times[states.index(i)]*1e-12)
        for i in states:
            for j in states:
                if transitions[states.index(i),states.index(j)]>0 and transitions[states.index(j),states.index(i)]>0:
                    Keq[states.index(i),states.index(j)]=rates[states.index(j),states.index(i)]/rates[states.index(i),states.index(j)]
        bsres.append([transitions,total_times*1e-12,rates,Keq])
    return bsres


print dwell_GLTPHOFC()

res=dwell_GLTPHOFC_bs()
print "---------Keq-----------"
z=np.zeros((1000,8,8))
for i in range(1,100): 
    z[i]=res[i][3]
print "means"
print np.mean(z,axis=0)
print "sds"
print np.std(z,axis=0)
print "---------rates-----------"
b=np.zeros((1000,8,8))
for i in range(1,100): 
    b[i]=res[i][2]
print "means"
print np.mean(b,axis=0)
print "sds"
print np.std(b,axis=0)