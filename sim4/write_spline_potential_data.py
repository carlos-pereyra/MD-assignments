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
r       = array('f',[0])
u       = array('f',[0])
u_old   = array('f',[0])
lj      = array('i',[0])
sp      = array('i',[0])
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

def lj_force(r):
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

def main():
    for n in range(0,nstep):
        #===============
        # R distance Loop
        #===============
        r[0]=0.001*(n+1)
        u[0],lj[0],sp[0]=lennard_jones_potential(r[0])
        u_old[0]=old_lennard_jones_potential(r[0])
        #print("r={:.3f} u={:.3f}".format(r[0],u[0]))
        tree.Fill()
    
    #===============
    # Write2Disk
    #===============
    tree.Write()  # everything respect to time
    file.Close()

if __name__ == "__main__":
    main()
