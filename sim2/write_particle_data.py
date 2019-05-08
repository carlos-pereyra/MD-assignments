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
from array import array

random.seed(int(datetime.now().strftime("%s")))
l=20
natoms=10
nsteps=10000


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
#       * Variables
file=TFile("data/test.root", "recreate" );
tree=TTree("time-sequence", "particle object storage")

t_dum=array('f',[0]*natoms)
ux_dum=array('f',[0]*natoms)
uy_dum=array('f',[0]*natoms)
x_dum=array('f',[0]*natoms)
y_dum=array('f',[0]*natoms)
fx_dum=array('f',[0]*natoms)
fy_dum=array('f',[0]*natoms)
vx_dum=array('f',[2]*natoms)
vy_dum=array('f',[3]*natoms)
dt_dum=array('f',[0.01]*natoms)
m_dum=array('f',[1]*natoms)
ke_dum=array('f',[0]*natoms)
pe_wall_dum=array('f',[0]*natoms)
pe_jone_dum=array('f',[0]*natoms)

t=array('f',[0])
ux=array('f',[0])
uy=array('f',[0])
x=array('f',[0])
y=array('f',[0])
fx=array('f',[0])
fy=array('f',[0])
vx=array('f',[0])
vy=array('f',[0])
id=array('i',[0])
dt=array('f',[0])
m=array('f',[1])
ke=array('f',[0])
pe_wall=array('f',[0])
pe_jone=array('f',[0])
tree.Branch('t',t,'t/F')
tree.Branch('ux',ux,'ux/F')
tree.Branch('uy',uy,'uy/F')
tree.Branch('x',x,'x/F')
tree.Branch('y',y,'y/F')
tree.Branch('fx',fx,'fx/F')
tree.Branch('fy',fy,'fy/F')
tree.Branch('vx',vx,'vx/F')
tree.Branch('vy',vy,'vy/F')

tree.Branch('id',id,'id/I')
tree.Branch('dt',dt,'dt/F')
tree.Branch('m',m,'m/F')
tree.Branch('ke',ke,'ke/F')
tree.Branch('pe_wall',pe_wall,'pe_wall/F')
tree.Branch('pe_jone',pe_jone,'pe_jone/F')

for i in range(0,len(x_dum)):
    for j in range(i,len(x_dum)):
        while (np.sqrt((x_dum[i]-x_dum[j])**2+(y_dum[i]-y_dum[j])**2)<1) or (np.sqrt((x_dum[i]-l)**2+(y_dum[i]-l)**2)<1):
            print("too close together i {} j {}".format(i,j))
            x_dum[i]=random.random()*l
            y_dum[i]=random.random()*l
            if(i==j):
                break

for i in range(0,len(x_dum)):
    for j in range(i,len(x_dum)):
        while (np.sqrt((x_dum[i]-x_dum[j])**2+(y_dum[i]-y_dum[j])**2)<1) or (np.sqrt((x_dum[i]-l)**2+(y_dum[i]-l)**2)<1):
            print("too close together i {} j {}".format(i,j))
            x_dum[i]=random.random()*l
            y_dum[i]=random.random()*l
            if(i==j):
                break

r=np.zeros((natoms,natoms))
u_lj_mat=np.zeros((natoms,natoms))


#   STEP3:
#       * Verlet Algorithm + Energy
for inc in range(0,1):
    dt_dum=array('f',[0.01*(inc+1)]*natoms)
    for n in range(0,nsteps):
        for p in range(0,natoms):
            ''' Lennard-Jones '''
            for j in range(p,natoms): #interaction potential
                r[p][j]=((x_dum[p]-x_dum[j])**2+(y_dum[p]-y_dum[j])**2)**1/2.
                r[j][p]=r[p][j]
                if p==j:
                    r[p][j]=0
                    u_lj_mat[p][j]=0
                else:
                    u_lj_mat[p][j]=lennard_jones_potential(r[p][j])
                    u_lj_mat[j][p]=u_lj_mat[p][j]
        
            ''' Each Atom's Instantaneous Energy '''
            pe_jone_dum[p]=np.sum(u_lj_mat[p])
            ke_dum[p]=0.5*m_dum[p]*(vx_dum[p]**2+vy_dum[p]**2)
            pe_wall_dum[p]=wall_potential(x_dum[p],y_dum[p])
            ''' Sum Energy for each Atom '''
            #u_lj_n[p]+=np.sum(u_lj_mat[p])
            #l_mv_n[p]+=l_mv
            #u_wl_n[p]+=u_wl
            ''' Verlet-Method '''
            t_dum[p]=dt_dum[p]*n
            ux_dum[p]=vx_dum[p]+dt_dum[p]*fx_dum[p]/(2*m_dum[p])
            uy_dum[p]=vy_dum[p]+dt_dum[p]*fy_dum[p]/(2*m_dum[p])
            x_dum[p]=x_dum[p]+ux_dum[p]*dt_dum[p]
            y_dum[p]=y_dum[p]+uy_dum[p]*dt_dum[p]
            fx_dum[p],fy_dum[p]=wall_force(x_dum[p],y_dum[p])
            vx_dum[p]=ux_dum[p]+dt_dum[p]*fx_dum[p]/(2*m_dum[p])
            vy_dum[p]=uy_dum[p]+dt_dum[p]*fy_dum[p]/(2*m_dum[p])

            #triggering condition
            if n<1000: #and ncycle[p]>0: # triggering condition
                t[0]=t_dum[p]
                ux[0]=ux_dum[p]
                uy[0]=uy_dum[p]
                x[0]=x_dum[p]
                y[0]=y_dum[p]
                fx[0]=fx_dum[p]
                fy[0]=fy_dum[p]
                vx[0]=vx_dum[p]
                vy[0]=vy_dum[p]
                id[0]=p
                dt[0]=dt_dum[p]
                ke[0]=ke_dum[p]
                pe_wall[0]=pe_wall_dum[p]
                pe_jone[0]=pe_jone_dum[p]
                tree.Fill()

        print("time: {}".format(t_dum[p]))


#   STEP3:
#       * write to disk
tree.Write()
file.Close()
