"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
from numpy import random
from datetime import datetime
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

random.seed(int(datetime.now().strftime("%s")))
l=20
natoms=10
nsteps=1000
ndt_steps=1

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

''' Tree 1 '''
tree=TTree("time_sequence", "particle object storage")
t_dum   =   array('f',[0]*natoms)
ux_dum  =   array('f',[0]*natoms)
uy_dum  =   array('f',[0]*natoms)
x_dum   =   array('f',[0]*natoms)
y_dum   =   array('f',[0]*natoms)
fx_dum  =   array('f',[0]*natoms)
fy_dum  =   array('f',[0]*natoms)
vx_dum  =   array('f',[2]*natoms)
vy_dum  =   array('f',[3]*natoms)
dt_dum  =   array('f',[0.01]*natoms)
m_dum   =   array('f',[1]*natoms)

t       =   array('f',[0])
ux      =   array('f',[0])
uy      =   array('f',[0])
x       =   array('f',[0])
y       =   array('f',[0])
fx      =   array('f',[0])
fy      =   array('f',[0])
vx      =   array('f',[2])
vy      =   array('f',[3])
n       =   array('i',[0])

id      =   array('i',[0])
dt      =   array('f',[0])
m       =   array('f',[1])

u_wall  =   array('f',[0])
ke      =   array('f',[0])
u_lj    =   array('f',[0])

tree.Branch('t',t,'t/F')
tree.Branch('ux',ux,'ux/F')
tree.Branch('uy',uy,'uy/F')
tree.Branch('x',x,'x/F')
tree.Branch('y',y,'y/F')
tree.Branch('fx',fx,'fx/F')
tree.Branch('fy',fy,'fy/F')
tree.Branch('vx',vx,'vx/F')
tree.Branch('vy',vy,'vy/F')
tree.Branch('n',n,'n/I')

tree.Branch('id',id,'id/I')
tree.Branch('dt',dt,'dt/F')
tree.Branch('m',m,'m/F')

tree.Branch('u_wall',u_wall,'u_wall/F')
tree.Branch('ke',ke,'ke/F')
tree.Branch('u_lj',u_lj,'u_lj/F')


''' Tree 2 '''
tree2=TTree("system_energy", "energy vs time")

u_wall_avg      =   array('f',[0])
ke_avg          =   array('f',[0])
u_lj_avg        =   array('f',[0])

tree2.Branch('t',t,'t/F')
tree2.Branch('n',n,'n/I')
tree2.Branch('u_wall_avg',u_wall_avg,'u_wall_avg/F')    #potential energy
tree2.Branch('ke_avg',ke_avg,'ke_avg/F')                #kinetic energy
tree2.Branch('u_lj_avg',u_lj_avg,'u_lj_avg/F')
tree2.Branch('dt',dt,'dt/F')

''' Tree 3 '''
tree3=TTree("avg_system_energy", "energy vs dt")

u_lj_avg_list   =   array('f')
ke_avg_list     =   array('f')
u_wall_avg_list =   array('f')

E_p_n           =   array('f',[0])
KE_p_n          =   array('f',[0])

tree3.Branch('E_p_n',E_p_n,'E_p_n/F')
tree3.Branch('KE_p_n',KE_p_n,'KE_p_n/F')
tree3.Branch('dt',dt,'dt/F')

''' Random Positions '''
for i in range(0,len(x_dum)):
    for j in range(i,len(x_dum)):
        while (np.sqrt((x_dum[i]-x_dum[j])**2+(y_dum[i]-y_dum[j])**2)<1) or (np.sqrt((x_dum[i]-l)**2+(y_dum[i]-l)**2)<1) or (np.sqrt(x_dum[i]**2+y_dum[i]**2)<1):
            print("too close together i {} j {}".format(i,j))
            x_dum[i]=random.random()*(l-1)
            y_dum[i]=random.random()*(l-1)
            if(i==j):
                break

for i in range(0,len(x_dum)):
    for j in range(i,len(x_dum)):
        while (np.sqrt((x_dum[i]-x_dum[j])**2+(y_dum[i]-y_dum[j])**2)<1) or (np.sqrt((x_dum[i]-l)**2+(y_dum[i]-l)**2)<1) or (np.sqrt(x_dum[i]**2+y_dum[i]**2)<1):
            print("too close together i {} j {}".format(i,j))
            x_dum[i]=random.random()*(l-1)
            y_dum[i]=random.random()*(l-1)
            if(i==j):
                break

r=np.zeros((natoms,natoms))
u_lj_mat=np.zeros((natoms,natoms))


#   STEP3:
#       * Verlet Algorithm + Energy
for inc in range(0,ndt_steps):
    dt_dum=array('f',[0.01*(inc+1)]*natoms)
    for time_step in range(0,nsteps):
        u_wall_list=array('f',[0]*natoms)
        ke_list=array('f',[0]*natoms)
        u_lj_list=array('f',[0]*natoms)
        for i in range(0,natoms):
            ''' Lennard-Jones '''
            for j in range(i,natoms):
                r[i][j]=((x_dum[i]-x_dum[j])**2+(y_dum[i]-y_dum[j])**2)**1/2.
                r[j][i]=r[i][j]
                if i==j:
                    r[i][j]=0
                    u_lj_mat[i][j]=0
                else:
                    u_lj_mat[i][j]=lennard_jones_potential(r[i][j])
                    u_lj_mat[j][i]=u_lj_mat[i][j]
        
            ''' Each Atom's Instantaneous Energy '''
            u_wall[0]=np.sum(u_lj_mat[i])
            ke[0]=0.5*m_dum[i]*(vx_dum[i]**2+vy_dum[i]**2)
            u_lj[0]=wall_potential(x_dum[i],y_dum[i])
            
            ''' Add Energy of each Atom to List '''
            #print("u_wall: {}\t ke: {}\t u_lj: {}\t n: {} dt: {}".format(u_wall[0],ke[0],u_lj[0],time_step,dt_dum[0]))
            u_lj_list[i]=np.sum(u_lj_mat[i])
            ke_list[i]=0.5*m_dum[i]*(vx_dum[i]**2+vy_dum[i]**2)
            u_wall_list[i]=wall_potential(x_dum[i],y_dum[i])
            
            ''' Verlet-Method '''
            t_dum[i]=dt_dum[i]*time_step
            ux_dum[i]=vx_dum[i]+dt_dum[i]*fx_dum[i]/(2*m_dum[i])
            uy_dum[i]=vy_dum[i]+dt_dum[i]*fy_dum[i]/(2*m_dum[i])
            x_dum[i]=x_dum[i]+ux_dum[i]*dt_dum[i]
            y_dum[i]=y_dum[i]+uy_dum[i]*dt_dum[i]
            fx_dum[i],fy_dum[i]=wall_force(x_dum[i],y_dum[i])
            vx_dum[i]=ux_dum[i]+dt_dum[i]*fx_dum[i]/(2*m_dum[i])
            vy_dum[i]=uy_dum[i]+dt_dum[i]*fy_dum[i]/(2*m_dum[i])

            #triggering condition
            if time_step<1000:
                t[0]=t_dum[i]
                ux[0]=ux_dum[i]
                uy[0]=uy_dum[i]
                x[0]=x_dum[i]
                y[0]=y_dum[i]
                fx[0]=fx_dum[i]
                fy[0]=fy_dum[i]
                vx[0]=vx_dum[i]
                vy[0]=vy_dum[i]
                n[0]=time_step
                
                id[0]=i
                dt[0]=dt_dum[i]
                m[0]=m_dum[i]
                
                tree.Fill()
    
        if time_step<1000:
            #print("time: {}".format(t_dum[p]))
            u_lj_avg[0]=np.sum(u_lj_list[0:natoms-3])/(natoms-3)
            ke_avg[0]=np.sum(ke_list[0:natoms-3])/(natoms-3)
            u_wall_avg[0]=np.sum(u_wall_list[0:natoms-3])/(natoms-3)
            
            u_lj_avg_list.append(u_lj_avg[0])
            ke_avg_list.append(ke_avg[0])
            u_wall_avg_list.append(u_wall_avg[0])
            
            tree2.Fill()

    E_p_n[0]=np.sum(u_lj_avg_list+u_wall_avg_list)/len(u_lj_avg_list) #total energy of system/n-steps
    KE_p_n[0]=np.sum(ke_avg_list)/len(ke_avg_list)
    dt[0]=dt_dum[0]
    tree3.Fill()


#   STEP3:
#       * write to disk
tree.Write()
tree2.Write()
tree3.Write()

file.Close()
