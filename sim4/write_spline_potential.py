""" ============================================================
*  author: carlos p
*  purpose: write data for visualizing the lj-spline potential
*
============================================================ """

from array import array
import numpy as np
from numpy import random
from datetime import datetime
random.seed(int(datetime.now().strftime("%s")))
import time
start_time=time.time()

import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
file=TFile("data/potential_e_r.root", "recreate" )
tree=TTree("timeSequence", "data with respect to time")
f       = array('f',[0])
fx       = array('f',[0])
fy       = array('f',[0])
f_old   = array('f',[0])
r       = array('f',[0])
u       = array('f',[0])
u_old   = array('f',[0])
lj      = array('i',[0])
sp      = array('i',[0])
tree.Branch('f',    f,      'f/F')
tree.Branch('fx',    fx,      'fx/F')
tree.Branch('fy',    fy,      'fy/F')
tree.Branch('f_old',f_old,  'f_old/F')
tree.Branch('r',    r,      'r/F')
tree.Branch('u',    u,      'u/F')
tree.Branch('u_old',u_old,  'u/F')
tree.Branch('lj',   lj,     'lj/I')
tree.Branch('sp',   sp,     'sp/I')
nstep=2000


def old_lennard_jones_potential(r):
    u_lj=(1/r)**12-2*(1/r)**6
    return u_lj

def lj_potential(r):
    #with spline compensation
    ri=1.1080
    rc=1.5475
    if r<ri:
        u_lj    = (1/r)**12-2*(1/r)**6
        lj      = 1
        sp      = 0
    if ri<=r<rc:
        u_lj    = -6.129*(r-rc)**2 + -4.655*(r-rc)**3
        lj      = 0
        sp      = 1
    elif r>=rc:
        u_lj    = 0
        lj      = 0
        sp      = 0

    return u_lj,lj,sp

'''
def lj_force(r):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    #r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    if r<ri:
        f       = (12/r**13 - 12/r**7) #neg sign added
        lj      = 1
        sp      = 0
    if ri<=r<rc:
        f       = 12.258*(r-rc)**1 + 13.966*(r-rc)**2
        lj      = 0
        sp      = 1
    elif r>=rc:
        f       = 0
        lj      = 0
        sp      = 0

    return f'''

def lj_force(xij,yij):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    ir = 1/np.sqrt(xij**2 + yij**2)
    
    if r<ri:
        fx = (12/r**13 - 12/r**7)*(xij*ir)
        fy = (12/r**13 - 12/r**7)*(xij*ir)
        lj = 1
        sp = 0

    if ri<=r<rc:
        fx = 12.258*((r-rc)**1)*(xij*ir) + 13.966*((r-rc)**2)*(xij*ir)
        fy = 12.258*((r-rc)**1)*(yij*ir) + 13.966*((r-rc)**2)*(yij*ir)
        lj = 0
        sp = 1
    elif r>=rc:
        fx = 0
        fy = 0
        lj = 0
        sp = 0
    
    return fx, fy

'''
def old_lj_force(r):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    #r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    f  = (12/r**13 - 12/r**7)

    return f'''

def old_lj_force(xij,yij):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    ir = 1/np.sqrt(xij**2 + yij**2)
    fx = (12/r**13 - 12/r**7)*(xij*ir)
    fy = (12/r**13 - 12/r**7)*(yij*ir)
    
    return fx,fy

def main():
    for n in range(0,nstep):
        #===============
        # R distance Loop
        #===============
        r[0]=0.01*(n+1)
        u[0],lj[0],sp[0] = lj_potential(r[0])
        u_old[0] = old_lennard_jones_potential(r[0])
        
        xij = 0.01*(n+1)
        yij = 0.00*(n+1)
        rij = np.sqrt(xij**2 + yij**2)
        
        fx1, fy1 = lj_force(xij,yij)
        f1 = np.sqrt(fx1**2 + fy1**2)
        f[0] = f1
        fx[0] = fx1
        fy[0] = fy1
        
        fx2,fy2 = old_lj_force(xij,yij)
        f2 = np.sqrt(fx2**2 + fy2**2)
        f_old[0] = fx2
        
        tree.Fill()
    
    #===============
    # Write2Disk
    #===============
    tree.Write()  # everything respect to time
    file.Close()

if __name__ == "__main__":
    main()
