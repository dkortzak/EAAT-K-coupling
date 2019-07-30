#!/usr/bin/env python
import numpy as np
from random import randint

def dwell_EAAT1():
    states=["0000","0010","0100","0110","1000","1010","1100","1110","0001","0011","0101","0111","1001","1011","1101","1111"]
    transitions=np.zeros(shape=(len(states),len(states)),dtype=int)
    total_times=np.zeros(len(states))
    rates=np.zeros(shape=(len(states),len(states)))
    Keq=np.zeros(shape=(len(states),len(states)))
    for chain in ["M","N","O"]:
     for run in range(1,104):
        in1=tuple(open("data/freeMD_EAAT1/run"+str(run)+"/states_"+chain+".xvg"))
        for i in range(2,len(in1)-1):
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
    return transitions,total_times*1e-12,rates,Keq

def dwell_EAAT1_bs():
    Nbs=1000
    states=["0000","0010","0100","0110","1000","1010","1100","1110","0001","0011","0101","0111","1001","1011","1101","1111"]
    transitions=np.zeros(shape=(len(states),len(states)),dtype=int)
    total_times=np.zeros(len(states))
    dwelltimes=[]
    for i in states:
        dwelltimes.append([])
    for i in range(len(states)):
        for j in states:
            dwelltimes[i].append([])  
    bschains=["M","N","O"]
    bsruns=range(1,104)
    bsres=[]
    for bss in range(0,Nbs):
        transitions=np.zeros(shape=(len(states),len(states)),dtype=int)
        total_times=np.zeros(len(states))
        rates=np.zeros(shape=(len(states),len(states)))
        Keq=np.zeros(shape=(len(states),len(states)))
        for count in range(len(bsruns)*3):
            k1=randint(0,2)
            k2=randint(1,len(bsruns))
            in1=tuple(open("data/freeMD_EAAT1/run"+str(bsruns[k2-1])+"/states_"+bschains[k1]+".xvg"))
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
                if transitions[states.index(i),states.index(j)]>1 and transitions[states.index(j),states.index(i)]>1:
                    Keq[states.index(i),states.index(j)]=rates[states.index(j),states.index(i)]/rates[states.index(i),states.index(j)]
        bsres.append([transitions,total_times*1e-12,rates,Keq])
    return bsres

print dwell_EAAT1()
res=dwell_EAAT1_bs()
print "---------Keq-----------"
b=np.zeros((1000,16,16))
for i in range(1,1000): 
    b[i]=res[i][3]
print "means"
print np.mean(b,axis=0)
print "sds"
print np.std(b,axis=0)
print "---------rates-----------"
b=np.zeros((1000,16,16))
for i in range(1,1000): 
    b[i]=res[i][2]
print "means"
print np.mean(b,axis=0)
print "sds"
print np.std(b,axis=0)