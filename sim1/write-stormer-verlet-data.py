"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

gROOT.ProcessLine(
"struct position {\
   std::vector<float>   x;\
   std::vector<float>   y;\
};\
struct velocity {\
   std::vector<float>   x;\
   std::vector<float>   y;\
};\
struct half {\
   std::vector<float>   x;\
   std::vector<float>   y;\
};\
struct force {\
   std::vector<float>   x;\
   std::vector<float>   y;\
};\
struct Particle {\
    std::vector<float> t;\
    position r;\
    velocity v;\
    half u;\
    force f;\
    Float_t k;\
    Float_t m;\
    Float_t dt;\
    Int_t id;\
    std::vector<float> KE_V;\
    std::vector<float> KE_U;\
    std::vector<float> PE_R;\
    Float_t KE_PE_V_SUM;\
    Float_t KE_PE_U_SUM;\
    Float_t KE_PE_V_AVG;\
    Float_t KE_PE_U_AVG;\
    Int_t ncycle;\
    Float_t freq;\
    Float_t freq_intpl;\
    Float_t omega_dt;\
    Int_t r0;\
    Int_t v0;\
    Float_t t1;\
    Float_t t2;\
    Int_t size;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

from ROOT import Box
from ROOT import Particle

box=Box()
r0=1
v0=0

#step0: allocate memory
natoms=60
m=[0]*natoms    #*
k=[0]*natoms    #*
dt=[0]*natoms   #*
t=[0]*natoms
u=[0]*natoms
r=[0]*natoms
f=[0]*natoms
v=[0]*natoms
#total average energy
t1=[0]*natoms
t1_intpl=[0]*natoms
t2=[0]*natoms
t2_intpl=[0]*natoms
KE_V=[0]*natoms
KE_U=[0]*natoms
PE_R=[0]*natoms
KE_PE_V_SUM=[0]*natoms
KE_PE_U_SUM=[0]*natoms
KE_PE_V_AVG=[0]*natoms
KE_PE_U_AVG=[0]*natoms
ncycle=[0]*natoms
freq=[0]*natoms
freq_intpl=[0]*natoms
omega_dt=[0]*natoms
box.atoms.resize(natoms)

file=TFile("data/stormer-verlet-dat-r{}-v{}.root".format(r0,v0), "recreate" );
tree=TTree("time-sequence", "particle object storage")
tree.Branch('box',box)

#step1: initial conditions
'''
for p in range(0,natoms):
    m[p]=1          #*
    k[p]=1          #*
    dt[p]=0.05*(1+p)    #*
    f0=-k[p]*r0
    #
    t[p]=0
    u[p]=v0+dt[p]*f0/(2*m[p])
    r[p]=r0
    f[p]=-k[p]*r0
    v[p]=v0
    #
    KE_V[p]=0.5*m[p]*v[p]*v[p]
    KE_U[p]=0.5*m[p]*v[p]*v[p]
    PE_R[p]=0.5*m[p]*r[p]*r[p]
'''

#step2: stormer-verlet
nsteps=1000
for n in range(0,nsteps):
    for p in range(0,natoms):
        rprev=r[p]     #r(i-1)
        tprev=t[p]  #t(i-1)

        if n==0:
            m[p]=1          #*
            k[p]=1          #*
            dt[p]=0.05*(1+p)    #*
            f0=-k[p]*r0
            #
            t[p]=0
            u[p]=v0+dt[p]*f0/(2*m[p])
            r[p]=r0
            f[p]=-k[p]*r0
            v[p]=v0
            #
            KE_V[p]=0.5*m[p]*v[p]*v[p]
            KE_U[p]=0.5*m[p]*v[p]*v[p]
            PE_R[p]=0.5*m[p]*r[p]*r[p]
        elif n>0:
            t[p]=dt[p]*(n+1)
            u[p]=v[p]+dt[p]*f[p]/(2*m[p])
            r[p]=r[p]+v[p]*dt[p]+(dt[p]**2)*f[p]/(2*m[p])
            f[p]=-k[p]*r[p]
            v[p]=u[p]+dt[p]*f[p]/(2*m[p])

        rnew=r[p]     #r(i)
        tnew=t[p]     #t(i)

        KE_V[p]=0.5*m[p]*v[p]*v[p]
        KE_U[p]=0.5*m[p]*u[p]*u[p]
        PE_R[p]=0.5*m[p]*r[p]*r[p]
        #print("KE_V: {} KE_U: {} PE_R: {}, NET_EN: {}".format(KE_V[p],KE_U[p],PE_R[p],KE_V[p]+PE_R[p]))
        omega_dt[p]=((k[p]/m[p])**(0.5))*(dt[p])

        if rnew<=0 and rprev>0:
            if t1[p]==0 and t2[p]==0:
                t1[p]=n*dt[p]
                #t1 approx
                t1_intpl[p]=(0-rnew)*(tnew-tprev)/(rnew-rprev)+tprev
                box.atoms[p].t1=t1[p]
            elif t1[p]!=0 and t2[p]==0:
                t2[p]=n*dt[p]
                #t2 approx
                t2_intpl[p]=(0-rnew)*(tnew-tprev)/(rnew-rprev)+tprev
                box.atoms[p].t2=t2[p]
                ncycle[p]=ncycle[p]+1
            elif t1[p]!=0 and t2[p]!=0:
                ncycle[p]=ncycle[p]+1
                freq[p]=(2*3.14)/(t2[p]-t1[p])
                freq_intpl[p]=(2*3.14)/(t2_intpl[p]-t1_intpl[p])
                #reset time points
                t1[p]=t2[p]
                t2[p]=n*dt[p]

        #if ncycle[p]>0:
        #    KE_PE_V_SUM[p]=KE_PE_V_SUM[p]+KE_V[p]+PE_R[p]
        #    KE_PE_U_SUM[p]=KE_PE_U_SUM[p]+KE_U[p]+PE_R[p]

        #if rnew<=0 and rprev>0 and t1[p]!=0 and t2[p]!=0:
        #    KE_PE_V_AVG[p]=KE_PE_V_SUM[p]/(ncycle[p])
        #    KE_PE_U_AVG[p]=KE_PE_U_SUM[p]/(ncycle[p])

        if n<1000: # triggering condition
            KE_PE_V_SUM[p]=KE_PE_V_SUM[p]+KE_V[p]+PE_R[p]
            KE_PE_U_SUM[p]=KE_PE_U_SUM[p]+KE_U[p]+PE_R[p]
            box.atoms[p].t.push_back(t[p])
            box.atoms[p].u.x.push_back(u[p])
            box.atoms[p].r.x.push_back(r[p])
            box.atoms[p].f.x.push_back(f[p])
            box.atoms[p].v.x.push_back(v[p])
            box.atoms[p].k=k[p]
            box.atoms[p].m=m[p]
            box.atoms[p].dt=dt[p]
            box.atoms[p].id=p
            box.atoms[p].KE_V.push_back(KE_V[p])
            box.atoms[p].KE_U.push_back(KE_U[p])
            box.atoms[p].PE_R.push_back(PE_R[p])
            box.atoms[p].KE_PE_V_SUM=KE_PE_V_SUM[p]
            box.atoms[p].KE_PE_U_SUM=KE_PE_U_SUM[p]
            #box.atoms[p].KE_PE_V_AVG=KE_PE_V_AVG[p]
            #box.atoms[p].KE_PE_U_AVG=KE_PE_U_AVG[p]
            box.atoms[p].ncycle=ncycle[p]
            box.atoms[p].freq=freq[p]
            box.atoms[p].freq_intpl=freq_intpl[p]
            box.atoms[p].omega_dt=omega_dt[p]
            box.atoms[p].r0=r0
            box.atoms[p].v0=v0
            box.atoms[p].size=n

#step3: write to disk
tree.Fill()
tree.Write()
file.Close()
