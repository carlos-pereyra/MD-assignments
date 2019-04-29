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
import numpy

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
random.seed(datetime.now())

l=20
natoms=5
nsteps=1000

rx0=[i+x for i,x in enumerate(sorted(random.sample(range(1,l-1), natoms)))]
ry0=[i+x for i,x in enumerate(sorted(random.sample(range(1,l-1), natoms)))]
#rx0=[5]*natoms
#ry0=[5]*natoms
vx0=[2]*natoms
vy0=[3]*natoms


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
#energy
r_mag=np.zeros((natoms,natoms))
u_lj=np.zeros((natoms,natoms))
ke_sum=0
u_wall_sum=0
u_lj_sum=0

box.atoms.resize(natoms)
file=TFile("test/particle_box_rx{}_ry{}_vx{}_vy{}_multiparticle.root".format(rx0[0],ry0[0],vx0[0],vy0[0]), "recreate" );
tree=TTree("time-sequence", "particle object storage")
tree.Branch('box',box)

#step1: stormer-verlet
for n in range(0,nsteps):
    ke_sum=0
    u_wall_sum=0
    u_lj_sum=0
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
            fx0,fy0=force_field(rx[p],ry[p])
            ux[p]=vx0[p]+dt[p]*fx0/(2*m[p])
            uy[p]=vy0[p]+dt[p]*fy0/(2*m[p])
            fx[p],fy[p]=force_field(rx[p],ry[p])
            vx[p]=vx0[p]
            vy[p]=vy0[p]

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
                    u_lj[p][j]=0 #phi_lj is a scalar!
                else:
                    u_lj[p][j]=potential_lennard_jones(r_mag[p][j])
                    u_lj[j][p]=u_lj[p][j]
            u_lj_sum+=np.sum(u_lj[p])

        #all other time
        elif n>0:
            t[p]=dt[p]*(n+1)
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
        if n<1000: #and ncycle[p]>0: # triggering condition
            box.atoms[p].rx0=rx0[p]
            box.atoms[p].ry0=ry0[p]
            box.atoms[p].vx0=vx0[p]
            box.atoms[p].vy0=vy0[p]
            box.atoms[p].t.push_back(t[p])
            box.atoms[p].u.x.push_back(ux[p])
            box.atoms[p].u.y.push_back(uy[p])
            box.atoms[p].r.x.push_back(rx[p])
            box.atoms[p].r.y.push_back(ry[p])
            box.atoms[p].f.x.push_back(fx[p])
            box.atoms[p].f.y.push_back(fy[p])
            box.atoms[p].v.x.push_back(vx[p])
            box.atoms[p].v.y.push_back(vy[p])
            box.atoms[p].k=k[p]
            box.atoms[p].m=m[p]
            box.atoms[p].dt=dt[p]
            box.atoms[p].id=p
            box.atoms[p].ke_avg.push_back(ke_sum/natoms)
            box.atoms[p].u_wall_avg.push_back(u_wall_sum/natoms) #potential energy per particle
            box.atoms[p].u_lj_avg.push_back(u_lj_sum/natoms)
            total_energy_avg=u_wall_sum/natoms+u_lj_sum/natoms+ke_sum/natoms
            box.atoms[p].total_energy_avg.push_back(total_energy_avg)
            '''
            for p in range(0,natoms):
                for j in range(p,natoms):
                    box.atoms[p].u_lj.push_back(u_lj[p][j])
                    box.atoms[j].u_lj.push_back(u_lj[j][p])
            '''

#print("u_lj: {}, u_lj[0][0]: {}".format(u_lj[0],u_lj[0][0]))

#step3: write to disk
tree.Fill()
tree.Write()
file.Close()

'''
test=np.zeros((5,5))
for i in range(0,5):
    for j in range(i,5):
        test[i][j]=i+0.5
        test[j][i]=-test[i][j]
        if i==j:
            test[i][j]=0
'''
