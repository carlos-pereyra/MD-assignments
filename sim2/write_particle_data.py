"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
from numpy import random
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
from datetime import datetime

gROOT.ProcessLine("struct Particle {\
    std::vector<float> t;\
    std::vector<float> ux;\
    std::vector<float> uy;\
    std::vector<float> x;\
    std::vector<float> y;\
    std::vector<float> fx;\
    std::vector<float> fy;\
    std::vector<float> vx;\
    std::vector<float> vy;\
    Float_t k;\
    Float_t m;\
    Float_t dt;\
    std::vector<float> ke_avg;\
    std::vector<float> u_wall_avg;\
    std::vector<float> u_lj_avg;\
    std::vector<float> total_energy_avg;\
    Int_t id;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

from ROOT import Box
from ROOT import Particle

box=Box()
random.seed(1) #datetime.now()

l=20
natoms=10
nsteps=10000

#   STEP0:
#        * Initial Conditions
#        * Global Variables
m=[1]*natoms
k=[1]*natoms
dt=[0.05]*natoms
t=[0]*natoms
ux=[0]*natoms
uy=[0]*natoms
rx=[0]*natoms
ry=[0]*natoms
i=0
while rx[natoms-1]==0 and ry[natoms-1]==0:
    rx[i]=random.randint(1,l)
    ry[i]=random.randint(1,l)
    for j in range(0,len(rx)):
        if rx[j]==rx[i] and ry[j]==ry[i] and i!=j:
            break
        else:
            i+=1
            break

fx=[0]*natoms
fy=[0]*natoms
vx=[2]*natoms
vy=[3]*natoms
#energy
r_mag=np.zeros((natoms,natoms))
u_lj=np.zeros((natoms,natoms))
ke_sum=0
u_wall_sum=0
u_lj_sum=0

#   STEP1:
#       * Simple Functions
#           * Wall Potential
#           * Lennard Jones Potential
#           * Wall Force
#           * Particle Interaction Force^* still in progress.
def potential_wall(xi,nu):
    u_left_boundary=((1/xi)**12-2*(1/xi)**6)
    u_right_boundary=((1/(l-xi))**12-2*(1/(l-xi))**6)
    u_bottom_boundary=((1/nu)**12-2*(1/nu)**6)
    u_top_boundary=((1/(l-nu))**12-2*(1/(l-nu))**6)
    u=u_left_boundary+u_right_boundary+u_bottom_boundary+u_top_boundary
    return u

def potential_lennard_jones(theta):
    uLJ=(1/theta)**12-2*(1/theta)**6
    return uLJ

def force_field(xi,nu):
    fx_left_boundary=12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary=-12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx=fx_left_boundary+fx_right_boundary

    fy_bottom_boundary=12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary=-12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy=fy_bottom_boundary+fy_top_boundary

    return fx,fy

box.atoms.resize(natoms)
#file=TFile("data/particle_box_rx{}_ry{}_vx{}_vy{}_multiparticle.root".format(rx[0],ry[0],vx[0],vy[0]), "recreate" );
file=TFile("data/particle_box_rx0_ry0_vx0_vy0_multiparticle.root", "recreate" );
tree=TTree("time-sequence", "particle object storage")
tree.Branch('box',box)

#   STEP2:
#       * Verlet Algorithm + Energy
#
for n in range(0,nsteps):
    ke_sum=0
    u_wall_sum=0
    u_lj_sum=0
    for p in range(0,natoms):
        #verlet
        t[p]=dt[p]*n
        ux[p]=vx[p]+dt[p]*fx[p]/(2*m[p])
        uy[p]=vy[p]+dt[p]*fy[p]/(2*m[p])
        rx[p]=rx[p]+vx[p]*dt[p]+(dt[p]**2)*fx[p]/(2*m[p])
        ry[p]=ry[p]+vy[p]*dt[p]+(dt[p]**2)*fy[p]/(2*m[p])
        #kinetic energy
        ke_sum+=0.5*m[p]*((vx[p]**2+vy[p]**2)**0.5)**2
        #wall potential
        u_wall=potential_wall(rx[p],ry[p])
        u_wall_sum+=u_wall #energy of all particles (at specific time step)
        #interaction potential
        for j in range(p,natoms):
            r_mag[p][j]=((rx[p]-rx[j])**2+(ry[p]-ry[j])**2)**1/2.
            r_mag[j][p]=-r_mag[p][j]
            if p==j:
                r_mag[p][j]=0
                u_lj[p][j]=0 #phi_lj is a scalar!
            else:
                u_lj[p][j]=potential_lennard_jones(r_mag[p][j])

        u_lj_sum+=np.sum(u_lj[p])
        fx[p],fy[p]=force_field(rx[p],ry[p])
        vx[p]=ux[p]+dt[p]*fx[p]/(2*m[p])
        vy[p]=uy[p]+dt[p]*fy[p]/(2*m[p])

        #triggering condition
        if n<1000:
            box.atoms[p].t.push_back(t[p])
            box.atoms[p].ux.push_back(ux[p])
            box.atoms[p].uy.push_back(uy[p])
            box.atoms[p].x.push_back(rx[p])
            box.atoms[p].y.push_back(ry[p])
            box.atoms[p].fx.push_back(fx[p])
            box.atoms[p].fy.push_back(fy[p])
            box.atoms[p].vx.push_back(vx[p])
            box.atoms[p].vy.push_back(vy[p])
            box.atoms[p].k=k[p]
            box.atoms[p].m=m[p]
            box.atoms[p].dt=dt[p]
            box.atoms[p].id=p
            box.atoms[p].ke_avg.push_back(ke_sum/natoms)
            box.atoms[p].u_wall_avg.push_back(u_wall_sum/natoms) #potential energy per particle
            box.atoms[p].u_lj_avg.push_back(u_lj_sum/natoms)
            total_energy_avg=u_wall_sum/natoms+u_lj_sum/natoms+ke_sum/natoms
            box.atoms[p].total_energy_avg.push_back(total_energy_avg)


#   STEP3:
#       * write to disk
tree.Fill()
tree.Write()
file.Close()
