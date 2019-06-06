""" ============================================================
*  author: carlos p
*  purpose: particle trees
*
============================================================ """

from array import array
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

alpha_list  = [0,0.1,1,10]
temp_list   = [0.1,1,10]
l=20
natom=5
nstep=int(10E0)
ndt=20
nrecord=int(1E3)
nal=int(len(alpha_list))
ntemp=int(len(temp_list))

#====================================
# tree1: TURFV
#====================================
t               = array('f',[0])
ux              = array('f',[0]*natom)
uy              = array('f',[0]*natom)
x               = array('f',[0]*natom)
y               = array('f',[0]*natom)
fx              = array('f',[0]*natom)
fy              = array('f',[0]*natom)
vx              = array('f',[2]*natom) #set initial velocity
vy              = array('f',[3]*natom) #set initial velocity
# system
n               = array('i',[0])
id              = array('i',[0]*natom)
dt              = array('f',[0])
dt_id           = array('i',[0])
m               = array('f',[1])
# energy per particle
ke_v            = array('f',[0]*natom)
ke_u            = array('f',[0]*natom)
u_wall          = array('f',[0]*natom)
u_jones         = array('f',[0]*natom)
# gaussian number
alpha           = array('f',[0])
alpha_id        = array('i',[0])
xi_1            = array('f',[0])
xi_2            = array('f',[0])
eta_1           = array('f',[0])
eta_2           = array('f',[0])
beta_1          = array('f',[0])
beta_2          = array('f',[0])
temp            = array('f',[0])
temp_id         = array('i',[0])
Kb              = array('f',[1])
#====================================
# tree2: Energy as a Function of time
#====================================
ke_v_sum_i      = array('f',[0])
ke_u_sum_i      = array('f',[0])
u_jones_sum_i   = array('f',[0])
u_wall_sum_i    = array('f',[0])
#
ke_v_avg_i      = array('f',[0])
ke_u_avg_i      = array('f',[0])
u_jones_avg_i   = array('f',[0])
u_wall_avg_i    = array('f',[0])
#
total_v_i       = array('f',[0]) #total average energy, velocity=v @t=n
total_u_i       = array('f',[0]) #total average energy, velocity=u @t=n
#====================================
# tree3: Average Energy as a Function of dt
#====================================
total_ke_v_sum = array('f',[0])
total_ke_u_sum = array('f',[0])
total_v_sum    = array('f',[0])
total_u_sum    = array('f',[0])
total_ke_v_avg  = array('f',[0])
total_ke_u_avg  = array('f',[0])
total_e_v_avg   = array('f',[0]) #average energy, velocity=v
total_e_u_avg   = array('f',[0]) #average energy, velocity=u
dt              = array('f',[0])
#====================================
# tree4: System Parameters
#====================================
vx0             =   array('i',[int(vx[0]) ])
vy0             =   array('i',[int(vy[0]) ])
atoms           =   array('i',[natom])
n_time_steps    =   array('i',[0])
n_dt_steps      =   array('i',[0])

file=TFile("data/Langevin_OpenBounds_AtomBox_{}atoms_{}tSteps_{}dts_vx0_{:.0f}_vy0_{:.0f}.root".format(natom,nstep,ndt,vx[0],vy[0]), "recreate" );

'''========================
    Tree 1
    ==========================='''
tree=TTree("timeSequence", "data with respect to time")
tree.Branch('atoms',atoms,'atoms/I')
tree.Branch('t',t,'t/F')
tree.Branch('ux',ux,'ux[atoms]/F')
tree.Branch('uy',uy,'uy[atoms]/F')
tree.Branch('x',x,'x[atoms]/F')
tree.Branch('y',y,'y[atoms]/F')
tree.Branch('fx',fx,'fx[atoms]/F')
tree.Branch('fy',fy,'fy[atoms]/F')
tree.Branch('vx',vx,'vx[atoms]/F')
tree.Branch('vy',vy,'vy[atoms]/F')
tree.Branch('id',id,'id[atoms]/I')
tree.Branch('ke_v',ke_v,'ke_v[atoms]/F')
tree.Branch('ke_u',ke_u,'ke_u[atoms]/F')
tree.Branch('u_wall',u_wall,'u_wall[atoms]/F')
tree.Branch('u_jones',u_jones,'u_jones[atoms]/F')
tree.Branch('nstep',n,'n/I')
tree.Branch('dt',dt,'dt/F')
tree.Branch('dt_id',dt_id,'dt_id/I')
tree.Branch('m',m,'m/F')
tree.Branch('alpha',alpha,'alpha/F')
tree.Branch('alpha_id',alpha_id,'alpha_id/I')
tree.Branch('temp',temp,'temp/F')
tree.Branch('temp_id',temp_id,'temp_id/I')

'''========================
    Tree 2
    ==========================='''
tree2=TTree("averagedEnergy", "average particle energy vs time")
tree2.Branch('ke_v',ke_v_avg_i,'ke_v_avg_i/F')    #ith kinetic energy
tree2.Branch('ke_u',ke_u_avg_i,'ke_u_avg_i/F')    #ith kinetic energy
tree2.Branch('u_jones',u_jones_avg_i,'u_jones_avg_i/F')
tree2.Branch('u_wall',u_wall_avg_i,'u_wall_avg_i/F')
tree2.Branch('total_v',total_v_i,'total_v_i/F')
tree2.Branch('total_u',total_u_i,'total_u_i/F')
tree2.Branch('t',t,'t/F')
tree2.Branch('dt',dt,'dt/F')
tree2.Branch('dt_id',dt_id,'dt_id/I')
tree2.Branch('alpha',alpha,'alpha/F')
tree2.Branch('alpha_id',alpha_id,'alpha_id/I')
tree2.Branch('temp',temp,'temp/F')
tree2.Branch('temp_id',temp_id,'temp_id/I')

'''========================
    Tree 3
    ==========================='''
tree3=TTree("avgSystemEnergy", "energy vs dt")
tree3.Branch('total_ke_v_avg',total_ke_v_avg,'total_ke_v_avg/F')
tree3.Branch('total_ke_u_avg',total_ke_u_avg,'total_ke_u_avg/F')
tree3.Branch('total_e_v_avg',total_e_v_avg,'total_e_v_avg/F')
tree3.Branch('total_e_u_avg',total_e_u_avg,'total_e_u_avg/F')
tree3.Branch('dt',dt,'dt/F')
tree3.Branch('dt_id',dt_id,'dt_id/I')
tree3.Branch('alpha',alpha,'alpha/F')
tree3.Branch('alpha_id',alpha_id,'alpha_id/I')
tree3.Branch('temp',temp,'temp/F')
tree3.Branch('temp_id',temp_id,'temp_id/I')

'''========================
    Tree 4
    ==========================='''
tree4=TTree("systemParameters", "system parameters")
tree4.Branch('vx0',vx0,'vx0/I')
tree4.Branch('vy0',vy0,'vy0/I')
tree4.Branch('natom',atoms,'atoms/I')
tree4.Branch('ntime',n_time_steps,'n_time_steps/I')
tree4.Branch('ndt',n_dt_steps,'n_dt_steps/I')
