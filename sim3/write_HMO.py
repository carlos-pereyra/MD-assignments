"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
from array import array
import random
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
import time
start_time=time.time()

r0=1
v0=0
k0=1
m0=1
dt0=0.05

alpha_list=[0,0.1,1,10]

Kb0=1.38E-23
Tp0=1

#step 0: variables
ndt=40
nsteps=int(10E3)
nrecord=int(1E3)
nal=len(alpha_list)

''' turfv '''
t               = array('f',[0])
u               = array('f',[0])
r               = array('f',[r0])
f               = array('f',[-k0*r0])
v               = array('f',[v0])
''' control variables '''
r0              = array('f',[r0])
v0              = array('f',[v0])
m               = array('f',[m0])
k               = array('f',[k0])
alpha           = array('f',[0])
alpha_id        = array('i',[0])
xi_1            = array('f',[0])
xi_2            = array('f',[0])
eta_1           = array('f',[0])
beta            = array('f',[0])
Tp              = array('f',[Tp0])
''' system '''
n               = array('f',[0])
dt              = array('f',[0])
id              = array('i',[0])
Kb              = array('f',[Kb0])
omega_dt        = array('f',[0])
omega_dt_id     = array('i',[0])
size            = array('i',[0])
''' energy vs time '''
ke_v            = array('f',[0])
ke_u            = array('f',[0])
u               = array('f',[0])

''' trees '''
file=TFile("data/HMO_wNoise_dat_r%.0f_v%.0f_n%.0f.root" % (r0[0],v0[0],nsteps), "recreate" );
tree=TTree("time","data_storage")
tree_freq=TTree("freq","data_storage")
tree_energy=TTree("energy","data_storage")
#tree_alpha=TTree("alpha","data_storage")

# turfv-branches
tree.Branch('t',t,'t/F')
tree.Branch('u',u,'u/F')
tree.Branch('r',r,'r/F')
tree.Branch('f',f,'f/F')
tree.Branch('v',v,'v/F')

# control variable-branches
tree.Branch('r0',r0,'r0/F')
tree.Branch('v0',v0,'v0/F')
tree.Branch('alpha',alpha,'alpha/F')
tree.Branch('alpha_id',alpha_id,'alpha_id/I')
tree.Branch('xi1',xi_1,'xi_1/F')
tree.Branch('xi2',xi_2,'xi_2/F')
tree.Branch('eta1',eta_1,'eta_1/F')
tree.Branch('beta',beta,'beta/F')

# system-branches
tree.Branch('n',n,'n/I')
tree.Branch('m',m,'m/F')
tree.Branch('k',k,'k/F')
tree.Branch('dt',dt,'dt/F')
tree.Branch('id',id,'id/I')
tree.Branch('Wdt',omega_dt,'omega_dt/F')
tree.Branch('WdtId',omega_dt_id,'omega_dt_id/I')
tree.Branch('size',size,'size/I')

# energy vs time-branches
tree.Branch('kev',ke_v,'ke_v/F')
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
tree_energy.Branch('alpha',alpha,'alpha/F')
tree_energy.Branch('alpha_id',alpha_id,'alpha_id/I')

#tree.Fill() #fill initial conditions

# step 1: langevin algorithm - studies
for n_al in range(0,nal):
    alpha[0]=alpha_list[n_al]
    alpha_id[0] =n_al

    for n_dt in range(ndt):
        omega_dt_id[0]=n_dt

        r[0]        =r0[0]
        v[0]        =v0[0]
        t[0]        =0

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

            a=(1 - alpha[0]*dt[0]/(2*m[0]) )/( 1 + alpha[0]*dt[0]/(2*m[0]) )
            b=1/(1 + alpha[0]*dt[0]/(2*m[0]) )

            '''=============== random number process ============='''
            xi_1[0]=random.uniform(0, 1) # xi1 calculation
            if xi_1[0]==0:
                print("change xi-1: {}".format(xi_1[0]))
                xi_1[0]=1-xi_1[0]

            xi_2[0]=random.uniform(0, 1) # xi2 calculation
            if xi_2[0]==0:
                print("change xi-2: {}".format(xi_2[0]))
                xi_2[0]=1-xi2[0]

            eta_1[0]=np.sqrt(-2*np.log(xi_1[0]))*np.cos(2*np.pi*xi_2[0])
            beta[0]=np.sqrt(2*alpha[0]*Kb[0]*Tp[0])*eta_1[0]
            '''==================================================='''

            t[0]=dt[0]*(i+1)
            u[0]=np.sqrt(b)*( v[0] + (dt[0]*f[0]+beta[0])/(2*m[0]) )
            r[0]=r[0]+np.sqrt(b)*dt[0]*u[0]
            f[0]=-k[0]*r[0]
            v[0]=(a/np.sqrt(b))*u[0]+(dt[0]*f[0]+beta[0])/(2*m[0])

            rnew=r[0]     #r(i)
            tnew=t[0]     #t(i)

            ke_v[0]=0.5*m[0]*v[0]*v[0]
            ke_u[0]=0.5*m[0]*u[0]*u[0]
            u[0]=0.5*m[0]*r[0]*r[0]

            totalE_sum_v[0]+=(ke_v[0]+u[0])
            totalE_sum_u[0]+=(ke_u[0]+u[0])

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
                n[0]=i          #nstep
                id[0]=n_dt      #id number
                size[0]=nsteps
                tree.Fill()

        totalE_avg_v[0]=totalE_sum_v[0]/float(nsteps)
        totalE_avg_u[0]=totalE_sum_u[0]/float(nsteps)
        tree_energy.Fill()

        dummy1=totalF_sum[0]/float(ncycle[0])
        dummy2=totalF_sum_i[0]/float(ncycle[0])
        totalF_avg[0]=dummy1
        totalF_avg_i[0]=dummy2
        tree_freq.Fill()

    print("alpha %s done --- %s seconds --- changing" % (alpha[0],time.time()-start_time))

print("DONE!")
#step3: write to disk
tree.Write()
tree_freq.Write()
tree_energy.Write()
file.Close()
