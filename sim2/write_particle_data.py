"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
from numpy import random
from datetime import datetime
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

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
std::vector<float> l_mv;\
std::vector<float> u_wl;\
std::vector<float> u_lj;\
std::vector<float> total_energy;\
std::vector<float> l_mv_avg;\
std::vector<float> u_wl_avg;\
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
time=datetime.now()
time=time.strftime("%s")
random.seed(int(time))

l=20
natoms=10
nsteps=100000

#   STEP0:
#        * Initial Conditions
#        * Global Variables
m=[1]*natoms
k=[1]*natoms
dt=[0.01]*natoms
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
r=np.zeros((natoms,natoms))
u_lj_mat=np.zeros((natoms,natoms))
l_mv_sum=np.zeros(natoms)
u_wl_sum=np.zeros(natoms)
u_lj_sum=np.zeros(natoms)

#file setup
box.atoms.resize(natoms)
file=TFile("data/particle_box_n{:d}.root".format(natoms), "recreate" );
tree=TTree("time_tree", "store_particle_data_in_time")
tree.Branch('box',box)

#   STEP1:
#       * Simple Functions
#           * Wall Potential
#           * Lennard Jones Potential
#           * Wall Force
#           * Particle Interaction Force^* still in progress.
def wall_potential(xi,nu):
    u_left_boundary     = ((1/xi)**12-2*(1/xi)**6)
    u_right_boundary    = ((1/(l-xi))**12-2*(1/(l-xi))**6)
    u_bottom_boundary   = ((1/nu)**12-2*(1/nu)**6)
    u_top_boundary      = ((1/(l-nu))**12-2*(1/(l-nu))**6)
    
    u=u_left_boundary+u_right_boundary+u_bottom_boundary+u_top_boundary
    return u

def lennard_jones_potential(r):
    u_lj=(1/r)**12-2*(1/r)**6
    return u_lj

def wall_force(xi,nu):
    fx_left_boundary    = 12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary   = -12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx                  = fx_left_boundary+fx_right_boundary

    fy_bottom_boundary  = 12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary     = -12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy                  = fy_bottom_boundary+fy_top_boundary

    return fx,fy


#   STEP2:
#       * Verlet Algorithm + Energy
#
for n in range(0,nsteps):
    for p in range(0,natoms):
        ''' ENERGY-CALCS'''
        #interaction potential
        for j in range(p,natoms):
            r[p][j]=((rx[p]-rx[j])**2+(ry[p]-ry[j])**2)**1/2.
            r[j][p]=r[p][j]
            if p==j:
                r[p][j]=0
                u_lj_mat[p][j]=0
            else:
                u_lj_mat[p][j]=lennard_jones_potential(r[p][j])
                u_lj_mat[j][p]=u_lj_mat[p][j]
    
        l_mv=0.5*m[p]*(vx[p]**2+vy[p]**2)
        u_wl=wall_potential(rx[p],ry[p])
        u_lj=np.sum(u_lj_mat[p])
        
        #energy sums
        u_lj_sum[p]+=np.sum(u_lj_mat[p])
        l_mv_sum[p]+=l_mv #kinetic energy
        u_wl_sum[p]+=u_wl
        
        ''' VERLET METHOD '''
        t[p]=dt[p]*n
        ux[p]=vx[p]+dt[p]*fx[p]/(2*m[p])
        uy[p]=vy[p]+dt[p]*fy[p]/(2*m[p])
        #rx[p]=rx[p]+vx[p]*dt[p]+(dt[p]**2)*fx[p]/(2*m[p])
        #ry[p]=ry[p]+vy[p]*dt[p]+(dt[p]**2)*fy[p]/(2*m[p])
        rx[p]=rx[p]+ux[p]*dt[p]
        ry[p]=ry[p]+uy[p]*dt[p]
        
        fx[p],fy[p]=wall_force(rx[p],ry[p])
        vx[p]=ux[p]+dt[p]*fx[p]/(2*m[p])
        vy[p]=uy[p]+dt[p]*fy[p]/(2*m[p])

        #triggering condition
        if n<10000:
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
            box.atoms[p].l_mv.push_back(l_mv) #kinetic energy
            box.atoms[p].u_wl.push_back(u_wl) #wall energy
            box.atoms[p].u_lj.push_back(u_lj) #lennard jones energy
            box.atoms[p].total_energy.push_back(l_mv+u_wl+u_lj)


#   STEP3:
#       * write to disk
tree.Fill()
tree.Write()
file.Close()
