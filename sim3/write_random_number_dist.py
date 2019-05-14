"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import random
import os
from datetime import datetime
from array import array

import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

random.seed(datetime.now())
ndata=int(10E7)
file=TFile("data/g_dist_n{}.root".format(ndata),"recreate");
tree=TTree("distribution", "data-storage")
xi1 = array('f',[0])
xi2 = array('f',[0])
nu1 = array('f',[0])
nu2 = array('f',[0])
g   = array('f',[0])
tree.Branch('xi1',xi1,'xi1/F')
tree.Branch('xi2',xi2,'xi2/F')
tree.Branch('nu1',nu1,'nu1/F')
tree.Branch('nu2',nu2,'nu2/F')
tree.Branch('g',g,'g/F')

#step0: function definitions
def g_distribution(nu1,nu2):
    g=np.exp(-(nu1**2+nu2**2)/2)/(2*np.pi)
    return g

#step1: generate set of random numbers
while ndata>0:
    ndata-=1
    if ndata%int(10E6)==0:
        print("{}".format(ndata))
    
    
    xi1=random.uniform(0, 1)
    if xi1==0:
        print("change xi1: {}".format(xi1))
        box.xi_1=1-box.xi_1
    
    xi2=random.uniform(0, 1)
    if xi2==0:
        print("change xi2: {}".format(xi2))
        xi2=1-xi2
    
    nu1=np.sqrt(-2*np.log(xi1))*np.cos(2*np.pi*xi2)
    nu2=np.sqrt(-2*np.log(xi2))*np.sin(2*np.pi*xi2)

    g=g_distribution(nu1,nu2)
    tree.Fill()

print("n: {}, DONE".format(ndata))

#step2: write to disk
tree.Write()
file.Close()
