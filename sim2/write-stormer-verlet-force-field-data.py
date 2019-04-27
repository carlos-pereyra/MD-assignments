"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
import random
from datetime import datetime

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
    Float_t rx0;\
    Float_t ry0;\
    Float_t vx0;\
    Float_t vy0;\
    position r;\
    velocity v;\
    half u;\
    force f;\
    Float_t k;\
    Float_t m;\
    Float_t dt;\
    std::vector<float> pe;\
    Float_t KE_V;\
    Float_t KE_U;\
    Float_t KE_PE_V_SUM;\
    Float_t KE_PE_U_SUM;\
    Float_t KE_PE_V_AVG;\
    Float_t KE_PE_U_AVG;\
    Int_t ncycle;\
    Float_t freq;\
    Float_t omega_dt;\
    Int_t id;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

from ROOT import Box
from ROOT import Particle

box=Box()
random.seed(datetime.now())

l=40
natoms=3
nsteps=1000
#rx0=5
#ry0=5
rx0=[2*i+x for i,x in enumerate(sorted(random.sample(range(1,l-1), natoms)))]
ry0=[5]*natoms
#vx0=-1
#vy0=1
vx0=[-1]*natoms
vy0=[1]*natoms


def potential_field(xi,nu):
    e0=1
    r0=1

    u_left_boundary=((1/xi)**12-2*(1/xi)**6)
    u_right_boundary=((1/(l-xi))**12-2*(1/(l-xi))**6)
    u_bottom_boundary=((1/nu)**12-2*(1/nu)**6)
    u_top_boundary=((1/(l-nu))**12-2*(1/(l-nu))**6)
    u=u_left_boundary+u_right_boundary+u_bottom_boundary+u_top_boundary

    return u


def force_field(xi,nu):
    e0=1
    r0=1

    fx_left_boundary=12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary=-12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx=fx_left_boundary+fx_right_boundary

    fy_bottom_boundary=12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary=-12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy=fy_bottom_boundary+fy_top_boundary

    return fx,fy


#step0: allocate memory
m=[0]*natoms    #*
k=[0]*natoms    #*
dt=[0]*natoms   #*
t=[0]*natoms
ux=[0]*natoms
uy=[0]*natoms
rx=[0]*natoms
ry=[0]*natoms
fx=[0]*natoms
fy=[0]*natoms
vx=[0]*natoms
vy=[0]*natoms
#total average energy
t1=[0]*natoms
t1_intpl=[0]*natoms
t2=[0]*natoms
t2_intpl=[0]*natoms
#energy
pe=[0]*natoms
ke_x=[0]*natoms
#KE_U=[0]*natoms
#KE_PE_V_SUM=[0]*natoms
#KE_PE_U_SUM=[0]*natoms
#KE_PE_V_AVG=[0]*natoms
#KE_PE_U_AVG=[0]*natoms
#ncycle=[0]*natoms
#freq=[0]*natoms
#freq_intpl=[0]*natoms
#omega_dt=[0]*natoms
box.atoms.resize(natoms)

file=TFile("data/stormer-verlet-dat-rx{}-vx{}.root".format(rx0[0],vx0[0]), "recreate" );
tree=TTree("time-sequence", "particle object storage")
tree.Branch('box',box)


#step1: stormer-verlet
for n in range(0,nsteps):
    for p in range(0,natoms):
        rprev=rx[p]     #r(i-1)
        tprev=t[p]      #t(i-1)

        #initial conditions
        if n==0:
            t[p]=0
            m[p]=1              #*
            k[p]=1              #*
            dt[p]=0.05*(1+p)    #*

            rx[p]=rx0[p]
            ry[p]=ry0[p]

            fx0, fy0=force_field(rx[p],ry[p])

            ux[p]=vx0[p]+dt[p]*fx0/(2*m[p])
            uy[p]=vy0[p]+dt[p]*fy0/(2*m[p])

            pe[p]=potential_field(rx[p],ry[p]) #phi is a scalar!
            fx[p],fy[p]=force_field(rx[p],ry[p])

            vx[p]=vx0[p]
            vy[p]=vy0[p]

        elif n>0:
            t[p]=dt[p]*(n+1)
            ux[p]=vx[p]+dt[p]*fx[p]/(2*m[p])
            uy[p]=vy[p]+dt[p]*fy[p]/(2*m[p])

            rx[p]=rx[p]+vx[p]*dt[p]+(dt[p]**2)*fx[p]/(2*m[p])
            ry[p]=ry[p]+vy[p]*dt[p]+(dt[p]**2)*fy[p]/(2*m[p])

            pe[p]=potential_field(rx[p],ry[p]) #phi is a scalar!
            fx[p],fy[p]=force_field(rx[p],ry[p])

            vx[p]=ux[p]+dt[p]*fx[p]/(2*m[p])
            vy[p]=uy[p]+dt[p]*fy[p]/(2*m[p])

        #print("n: {}, p: {}, rx: {}, ry: {}".format(n,p,rx[p],ry[p]))
        if n<1000: #and ncycle[p]>0: # triggering condition
            box.atoms[p].rx0=rx0[p]

            box.atoms[p].vx0=vx0[p]
            box.atoms[p].t.push_back(t[p])
            box.atoms[p].u.x.push_back(ux[p])
            box.atoms[p].u.y.push_back(uy[p])
            box.atoms[p].r.x.push_back(rx[p])
            box.atoms[p].r.y.push_back(ry[p])
            box.atoms[p].f.x.push_back(fx[p])
            box.atoms[p].f.y.push_back(fy[p])
            box.atoms[p].v.x.push_back(vx[p])
            box.atoms[p].v.y.push_back(vy[p])
            box.atoms[p].pe.push_back(pe[p])
            box.atoms[p].k=k[p]
            box.atoms[p].m=m[p]
            box.atoms[p].dt=dt[p]
            box.atoms[p].id=p
            #box.atoms[p].KE_PE_V_SUM=KE_PE_V_SUM[p]
            #box.atoms[p].KE_PE_U_SUM=KE_PE_U_SUM[p]
            #box.atoms[p].KE_PE_V_AVG=KE_PE_V_AVG[p]
            #box.atoms[p].KE_PE_U_AVG=KE_PE_U_AVG[p]
            #box.atoms[p].ncycle=ncycle[p]
            #box.atoms[p].freq=freq[p]
            #box.atoms[p].freqIntpl=freq_intpl[p]
            #box.atoms[p].omega_dt=omega_dt[p]


#step3: write to disk
tree.Fill()
tree.Write()
file.Close()
