"""
*  author: carlos p
*  purpose: solve particle in a box
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
Int_t id;\
std::vector<float> l_mv;\
std::vector<float> u_wl;\
std::vector<float> u_lj;\
std::vector<float> total_u_per_atom;\
std::vector<float> total_u;\
std::vector<float> dtt;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

from ROOT import Box
box=Box()

random.seed(int(datetime.now().strftime("%s")))
l=20
natoms=20
nsteps=10000
nstore=8000
dt=np.ones(natoms)*0.01

#   STEP0:
#        * Initial Conditions
#        * Global Variables
t=np.zeros(natoms)
ux=np.zeros(natoms)
uy=np.zeros(natoms)
rx=random.rand(natoms)*l
ry=random.rand(natoms)*l
for i in range(0,len(rx)):
    for j in range(i,len(rx)):
        while (np.sqrt((rx[i]-rx[j])**2+(ry[i]-ry[j])**2)<1 and i!=j) or (np.sqrt((rx[i]-l)**2+(ry[i]-l)**2)<1):
            print("too close together i {} j {}".format(i,j))
            rx[i]=random.random()*l
            ry[i]=random.random()*l

fx=np.zeros(natoms)
fy=np.zeros(natoms)
vx=np.ones(natoms)*2
vy=np.ones(natoms)*3
#energy
r=np.zeros((natoms,natoms))
u_lj_mat=np.zeros((natoms,natoms))
l_mv_n=np.zeros(natoms)
u_wl_n=np.zeros(natoms)
u_lj_n=np.zeros(natoms)
u_per_atom=0
u_per_time=0
#system constants
m=np.ones(natoms)
k=np.zeros(natoms)
#file out
box.atoms.resize(natoms)
file=TFile("data/particle_box_n{:d}_dt1_{:d}.root".format(natoms,int(1/dt[0])), "recreate" );
tree=TTree("time_tree", "store_particle_data")
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
    l_mv_n=np.zeros(natoms)
    u_wl_n=np.zeros(natoms)
    u_lj_n=np.zeros(natoms)
    potential_energy_over_time=0
    for p in range(0,natoms):
        ''' Lennard-Jones '''
        for j in range(p,natoms): #interaction potential
            r[p][j]=((rx[p]-rx[j])**2+(ry[p]-ry[j])**2)**1/2.
            r[j][p]=r[p][j]
            if p==j:
                r[p][j]=0
                u_lj_mat[p][j]=0
            else:
                u_lj_mat[p][j]=lennard_jones_potential(r[p][j])
                u_lj_mat[j][p]=u_lj_mat[p][j]
    
        ''' Each Atom's Instantaneous Energy '''
        u_lj=np.sum(u_lj_mat[p])
        l_mv=0.5*m[p]*(vx[p]**2+vy[p]**2)
        u_wl=wall_potential(rx[p],ry[p])
        ''' Sum Energy for each Atom '''
        u_lj_n[p]+=np.sum(u_lj_mat[p])
        l_mv_n[p]+=l_mv
        u_wl_n[p]+=u_wl
        ''' Verlet-Method '''
        t[p]=dt[p]*n
        ux[p]=vx[p]+dt[p]*fx[p]/(2*m[p])
        uy[p]=vy[p]+dt[p]*fy[p]/(2*m[p])
        rx[p]=rx[p]+ux[p]*dt[p]
        ry[p]=ry[p]+uy[p]*dt[p]
        fx[p],fy[p]=wall_force(rx[p],ry[p])
        vx[p]=ux[p]+dt[p]*fx[p]/(2*m[p])
        vy[p]=uy[p]+dt[p]*fy[p]/(2*m[p])
        
        #triggering condition
        if n<nstore:
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
            box.atoms[p].total_u_per_atom.push_back(l_mv+u_wl+u_lj)
            potential_energy_over_time+=l_mv+u_wl+u_lj

    if n==nstore:
        box.atoms[p].total_u.push_back(potential_energy_over_time/(n+1))
        box.atoms[p].dtt.push_back(dt[0])

#   STEP3:
#       * write to disk
tree.Fill()
tree.Write()
file.Close()

#print("rx: {}, ry: {}, vx: {}, vy: {}, fx: {}, fy: {}, ux: {}, uy: {} id: {} dt: {} m: {}".format(rx[p],ry[p],vx[p],vy[p],fx[p],fy[p],ux[p],uy[p],p,dt[p],m[p]))
