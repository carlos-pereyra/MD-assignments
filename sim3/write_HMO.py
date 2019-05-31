"""
*  author: carlos p
*  purpose: solve langevin equation with random noise and damping
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
dt0=0.01

alpha_list=[0,0.1,1,10]
temp_list=[0.1,1,10]
Kb0=1
#Kb0=1.38E-23

#step 0: variables
ndt=40
nsteps=int(1E3)
nrecord=int(1E3)
nal=int(len(alpha_list))
ntemp=int(len(temp_list))

'''==============================
    variables
    =============================='''
t               = array('f',[0])
u               = array('f',[0])
r               = array('f',[r0])
f               = array('f',[-k0*r0])
v               = array('f',[v0])
'''
    control variables
    '''
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
temp            = array('f',[0])
temp_id         = array('i',[0])
'''
    system
    '''
n               = array('f',[0])
dt              = array('f',[0])
dt_id           = array('i',[0])
Kb              = array('f',[Kb0])
Wdt             = array('f',[0])
size            = array('i',[0])
'''
    energy
    '''
ke_v            = array('f',[0])
ke_u            = array('f',[0])
pe              = array('f',[0])
'''
    frequency
    '''
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
'''
    average energy
    '''
totalE_sum_v    = array('f',[0])
totalE_sum_u    = array('f',[0])
totalE_avg_v    = array('f',[0])
totalE_avg_u    = array('f',[0])

'''==============================
    trees
    =============================='''
file=TFile("data/hmo/HMO_dat_r{:.0f}_v{:.0f}_nstep{}_ndt{}_nalpha{}.root".format(r0[0],v0[0],nsteps,ndt,nal), "recreate" );
tree=TTree("timeDat","data_storage")
tree2=TTree("freqDat","frequency and dt study")    # frequency vs dt
tree3=TTree("energyDat","energy and dt study")     # energy vs dt
#tree_alpha=TTree("alpha","data_storage")

'''==============================
    tree 1: u,v,f,ke,pe and time
    =============================='''
tree.Branch('t',t,'t/F')
tree.Branch('u',u,'u/F')
tree.Branch('r',r,'r/F')
tree.Branch('f',f,'f/F')
tree.Branch('v',v,'v/F')
# control variable-branches
tree.Branch('r0',r0,'r0/F')
tree.Branch('v0',v0,'v0/F')
tree.Branch('xi1',xi_1,'xi_1/F')
tree.Branch('xi2',xi_2,'xi_2/F')
tree.Branch('eta1',eta_1,'eta_1/F')
tree.Branch('beta',beta,'beta/F')
# system-branches
tree.Branch('n',n,'n/I')
tree.Branch('m',m,'m/F')
tree.Branch('k',k,'k/F')
tree.Branch('Wdt',Wdt,'Wdt/F')
tree.Branch('size',size,'size/I')
# energy vs time-branches
tree.Branch('kev',ke_v,'ke_v/F')
tree.Branch('keu',ke_u,'ke_u/F')
tree.Branch('pe',pe,'pe/F')
# auxiliary
tree.Branch('dt',dt,'dt/F')
tree.Branch('dt_id',dt_id,'dt_id/I')
tree.Branch('alpha',alpha,'alpha/F')
tree.Branch('alpha_id',alpha_id,'alpha_id/I')
tree.Branch('temp',temp,'temp/F')
tree.Branch('temp_id',temp_id,'temp_id/I')

'''==============================
    tree 2: frequency and dt
    =============================='''
tree2.Branch('f_sum',totalF_sum,'totalF_sum/F')
tree2.Branch('f_avg',totalF_avg,'totalF_avg/F')
tree2.Branch('f_i_sum',totalF_sum_i,'totalF_sum_i/F')
tree2.Branch('f_i_avg',totalF_avg_i,'totalF_avg_i/F')
tree2.Branch('ncycle',ncycle,'ncycle/F')
tree2.Branch('dt',dt,'dt/F')
tree2.Branch('Wdt',Wdt,'Wdt/F')
tree2.Branch('alpha_id',alpha_id,'alpha_id/I')

'''==============================
    tree 3: energy and dt and alpha and T
    =============================='''
tree3.Branch('totalE_sum_v',totalE_sum_v,'totalE_sum_v/F')
tree3.Branch('totalE_sum_u',totalE_sum_u,'totalE_sum_u/F')
tree3.Branch('totalE_avg_v',totalE_avg_v,'totalE_avg_v/F')
tree3.Branch('totalE_avg_u',totalE_avg_u,'totalE_avg_u/F')
tree3.Branch('dt',dt,'dt/F')
tree3.Branch('Wdt',Wdt,'Wdt/F')
tree3.Branch('alpha',alpha,'alpha/F')
tree3.Branch('alpha_id',alpha_id,'alpha_id/I')
tree3.Branch('temp',temp,'temp/F')
tree3.Branch('temp_id',temp_id,'temp_id/I')

# step 1: langevin algorithm - studies
for n_temp in range(0,ntemp):
    temp[0]=temp_list[n_temp]
    temp_id[0]=n_temp
    
    for n_al in range(0,nal):
        alpha[0]=alpha_list[n_al]
        alpha_id[0] =n_al

        for n_dt in range(ndt):
            dt_id[0]    = n_dt
            r[0]        = r0[0]
            v[0]        = v0[0]
            t[0]        = 0
            '''===============
                avg freq
                ==============='''
            t1[0]       =0
            t2[0]       =0
            t1_i[0]     =0
            t2_i[0]     =0
            ncycle[0]   =0
            totalF_sum[0]=0
            totalF_avg[0]=0
            totalF_sum_i[0]=0
            totalF_avg_i[0]=0

            '''===============
                avg energy
                ==============='''
            totalE_sum_v[0]=0
            totalE_sum_u[0]=0
            totalE_avg_v[0]=0
            totalE_avg_u[0]=0

            for n_step in range(0,nsteps):
                rprev=r[0]  #r(i-1)
                tprev=t[0]  #t(i-1)

                dt[0]=dt0*(1+n_dt)    #*
                Wdt[0]=dt[0]*(k[0]/m[0])**(0.5)

                a=(1-alpha[0]*dt[0]/(2*m[0]) )/( 1 + alpha[0]*dt[0]/(2*m[0]) )
                b=1/(1+alpha[0]*dt[0]/(2*m[0]) )

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
                beta[0]=np.sqrt(2*alpha[0]*Kb[0]*temp[0]*dt[0])*eta_1[0]

                '''===============
                    langevin eqn
                    ==============='''
                t[0]=dt[0]*(n_step+1)
                u[0]=np.sqrt(b)*( v[0] + (dt[0]*f[0]+beta[0])/(2*m[0]) )
                r[0]=r[0]+np.sqrt(b)*dt[0]*u[0]
                f[0]=-k[0]*r[0]
                v[0]=(a/np.sqrt(b))*u[0]+(dt[0]*f[0]+beta[0])/(2*m[0])
                ke_v[0]=0.5*m[0]*v[0]*v[0]
                ke_u[0]=0.5*m[0]*u[0]*u[0]
                pe[0]=0.5*m[0]*r[0]*r[0]
                totalE_sum_v[0]+=(ke_v[0]+pe[0])
                totalE_sum_u[0]+=(ke_u[0]+pe[0])
                #print("v {} ke {} pe {} total {}".format(v[0],ke_v[0],pe[0],totalE_sum_v[0]))
                
                '''===============
                   frequency calc
                   ==============='''
                rnew=r[0]     #r(i)
                tnew=t[0]     #t(i)
                if rnew<=0 and rprev>0:
                    if t1[0]==0 and t2[0]==0:
                        t1[0]=n_step*dt[0]
                        t1_i[0]=(0-rnew)*(tnew-tprev)/(rnew-rprev)+tprev
                    elif t1[0]!=0 and t2[0]==0:
                        t2[0]=n_step*dt[0]
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
                        t2[0]=n_step*dt[0]

                '''===============
                    record data
                    ==============='''
                if n_step<nrecord: # triggering condition
                    n[0]=n_step         #nstep
                    #dt_id[0]=n_dt      #id number
                    size[0]=nsteps
                    tree.Fill()

            totalE_avg_v[0]=totalE_sum_v[0]/float(nsteps)
            totalE_avg_u[0]=totalE_sum_u[0]/float(nsteps)
            tree3.Fill() #compare E vs alpha vs temp

            #totalF_avg[0]=totalF_sum[0]/float(ncycle[0])
            #totalF_avg_i[0]=totalF_sum_i[0]/float(ncycle[0])
            #tree2.Fill()

        print("temp %s alpha %s done --- %s seconds --- changing" % (temp[0],alpha[0],time.time()-start_time))




print("DONE!")
#step3: write to disk
tree.Write()
tree2.Write() #
tree3.Write() #
file.Close()
