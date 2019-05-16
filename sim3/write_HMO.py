"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

r0=0
v0=1
k0=1
m0=1
dt0=0.05

#step 0: variables
ndt=40
nsteps=1000
nrecord=1000
''' turfv '''
t               = array('f',[0])
u               = array('f',[0])
r               = array('f',[r0])
f               = array('f',[-k0*r0])
v               = array('f',[v0])
''' system '''
r0              = array('f',[r0])
v0              = array('f',[v0])
m               = array('f',[m0])
k               = array('f',[k0])
dt              = array('f',[0])
omega_dt        = array('f',[0])
n               = array('f',[0])
id              = array('i',[0])
''' energy vs time '''
ke_v            = array('f',[0])
ke_u            = array('f',[0])
u               = array('f',[0])

''' trees '''
file=TFile("data/HMO_dat_r%.0f_v%.0f.root" % (r0[0],v0[0]), "recreate" );
tree=TTree("time","data_storage")
tree_freq=TTree("freq","data_storage")
tree_energy=TTree("energy","data_storage")
tree.Branch('t',t,'t/F') # turfv-branches
tree.Branch('u',u,'u/F')
tree.Branch('r',r,'r/F')
tree.Branch('f',f,'f/F')
tree.Branch('v',v,'v/F')
tree.Branch('r0',r0,'r0/F') # system-branches
tree.Branch('v0',v0,'v0/F')
tree.Branch('m',m,'m/F')
tree.Branch('k',k,'k/F')
tree.Branch('dt',dt,'dt/F')
tree.Branch('Wdt',omega_dt,'omega_dt/F')
tree.Branch('nstep',n,'n/F')
tree.Branch('id',id,'n/I')
tree.Branch('kev',ke_v,'ke_v/F') # energy vs time-branches
tree.Branch('keu',ke_u,'ke_u/F')
tree.Branch('u',u,'u/F')

''' freq vs omega_dt '''
t1              = array('f',[0])
t2              = array('f',[0])
t1_i            = array('f',[0])
t2_i            = array('f',[0])
freq            = array('f',[0])
freq_intpl      = array('f',[0])
freq_sum        = array('f',[0])
freq_sum_i      = array('f',[0])

ncycle          = array('f',[0])
totalF_sum      = array('f',[0])
totalF_avg      = array('f',[0])
totalF_sum_i    = array('f',[0])
totalF_avg_i    = array('f',[0])

tree_freq.Branch('f_sum',totalF_sum,'totalF_sum/F')
tree_freq.Branch('f_avg',totalF_avg,'totalF_avg/F')
tree_freq.Branch('f_i_sum',totalF_sum_i,'totalF_sum_i/F')
tree_freq.Branch('f_i_avg',totalF_avg_i,'totalF_avg_i/F')
tree_freq.Branch('ncycle',ncycle,'ncycle/F')
tree_freq.Branch('dt',dt,'dt/F')
tree_freq.Branch('Wdt',omega_dt,'omega_dt/F')

''' energy vs dt-branches '''
totalE_sum_v    = array('f',[0])
totalE_sum_u    = array('f',[0])
totalE_avg_v    = array('f',[0])
totalE_avg_u    = array('f',[0])

tree_energy.Branch('totalE_sum_v',totalE_sum_v,'totalE_sum_v/F')
tree_energy.Branch('totalE_sum_u',totalE_sum_u,'totalE_sum_u/F')
tree_energy.Branch('totalE_avg_v',totalE_avg_v,'totalE_avg_v/F')
tree_energy.Branch('totalE_avg_u',totalE_avg_u,'totalE_avg_u/F')
tree_energy.Branch('dt',dt,'dt/F')
tree_energy.Branch('Wdt',omega_dt,'omega_dt/F')

tree.Fill() #fill initial conditions

#step 1: langevin algorithm
for n_dt in range(ndt):
    r[0]        = r0[0]
    v[0]        = v0[0]
    t[0]        = 0

    ''' avg freq vs omega_dt '''
    t1[0]       =0
    t2[0]       =0
    t1_i[0]     =0
    t2_i[0]     =0
    ncycle[0]   =0
    
    totalF_sum[0]=0
    totalF_avg[0]=0
    totalF_sum_i[0]=0
    totalF_avg_i[0]=0
    
    ''' avg energy vs omega_dt '''
    totalE_sum_v[0]=0
    totalE_sum_u[0]=0
    totalE_avg_v[0]=0
    totalE_avg_u[0]=0
    
    for i in range(0,nsteps):
        rprev=r[0]  #r(i-1)
        tprev=t[0]  #t(i-1)

        dt[0]=dt0*(1+n_dt)    #*
        omega_dt[0]=dt[0]*(k[0]/m[0])**(0.5)
        
        t[0]=dt[0]*(i+1)
        u[0]=v[0]+dt[0]*f[0]/(2*m[0])
        r[0]=r[0]+v[0]*dt[0]+(dt[0]**2)*f[0]/(2*m[0])
        f[0]=-k[0]*r[0]
        v[0]=u[0]+dt[0]*f[0]/(2*m[0])

        rnew=r[0]     #r(i)
        tnew=t[0]     #t(i)

        ke_v[0]=0.5*m[0]*v[0]*v[0]
        ke_u[0]=0.5*m[0]*u[0]*u[0]
        u[0]=0.5*m[0]*r[0]*r[0]
        #print("KE_V: {} KE_U: {} PE_R: {}".format(ke_v[0],ke_u[0],u[0]))

        if rnew<=0 and rprev>0:
            if t1[0]==0 and t2[0]==0:
                t1[0]=i*dt[0]
                t1_i[0]=(0-rnew)*(tnew-tprev)/(rnew-rprev)+tprev
            elif t1[0]!=0 and t2[0]==0:
                t2[0]=i*dt[0]
                t2_i[0]=(0-rnew)*(tnew-tprev)/(rnew-rprev)+tprev
                ncycle[0]+=1
            elif t1[0]!=0 and t2[0]!=0:
                ncycle[0]+=1
                freq[0]=(2*3.14)/(t2[0]-t1[0])
                freq_intpl[0]=(2*3.14)/(t2_i[0]-t1_i[0])
                
                totalF_sum[0]+=freq[0]
                #print("sum of frequencies={}, ncycle={}".format(totalF_sum[0],ncycle[0]))
                totalF_sum_i[0]+=freq_intpl[0]
                
                t1[0]=t2[0] #reset time points
                t2[0]=i*dt[0]

        if i<nrecord: # triggering condition
            if i==990: print("ndt: {} dt: {} id: {}".format(n_dt,dt[0],id[0]))
            n[0]=i
            id[0]=n_dt
            tree.Fill()
            
            totalE_sum_v[0]+=(ke_v[0]+u[0])
            totalE_sum_u[0]+=(ke_u[0]+u[0])
        
    totalE_avg_v[0]=totalE_sum_v[0]/nrecord
    totalE_avg_u[0]=totalE_sum_u[0]/nrecord
    tree_energy.Fill()

    dummy1=totalF_sum[0]/float(ncycle[0])
    dummy2=totalF_sum_i[0]/float(ncycle[0])
    totalF_avg[0]=dummy1
    totalF_avg_i[0]=dummy2
    tree_freq.Fill()


print("DONE!")
#step3: write to disk
tree.Write()
tree_freq.Write()
tree_energy.Write()
file.Close()
