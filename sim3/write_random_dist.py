"""
*  author: carlos p
*  purpose: generate random distribution from uniform dist.
*
"""

import numpy as np
import random
import os
from datetime import datetime
from array import array
import time
start_time=time.time()

import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

random.seed(datetime.now())
DBG=0
ndata=int(10E3)
nbins=int(1E3)
file=TFile("data/data_n{}.root".format(ndata),"recreate");
tree=TTree("distribution", "data-storage")
tree_xi_1=TTree("xi_1_distribution", "data-storage")
tree_xi_2=TTree("xi_2_distribution", "data-storage")
tree_nu_1=TTree("nu_1_distribution", "data-storage")
tree_nu_2=TTree("nu_2_distribution", "data-storage")

file_dat=open("data/data_n{}.dat".format(ndata), "w")

n           = array('i',[0])
xi1         = array('f',[0])
xi1_list    = array('f',[0]*nbins)
xi2         = array('f',[0])
xi2_list    = array('f',[0]*nbins)
nu1         = array('f',[0])
nu1_list    = array('f',[0]*nbins)
nu2         = array('f',[0])
nu2_list    = array('f',[0]*nbins)
g           = array('f',[0])
g_list      = array('f',[0]*nbins)

tree.Branch('n',n,'n/I')
tree.Branch('xi1',xi1,'xi1/F')
tree.Branch('xi2',xi2,'xi2/F')
tree.Branch('nu1',nu1,'nu1/F')
tree.Branch('nu2',nu2,'nu2/F')
tree.Branch('g',g,'g/F')
tree_xi_1.Branch('xi1',xi1,'xi1/F')
tree_xi_2.Branch('xi2',xi2,'xi2/F')
tree_nu_1.Branch('nu1',nu1,'nu1/F')
tree_nu_2.Branch('nu2',nu2,'nu2/F')

#step0: function definitions
def g_distribution(nu1,nu2):
    gdist=np.exp(-(nu1**2+nu2**2)/2)/(2*np.pi)
    return gdist

#step1: generate set of random numbers
while ndata>0:
    ndata-=1
    if ndata%int(10E6)==0:
        print("--- %s seconds --- %s" % (time.time() - start_time,ndata))

    ''' xi1 calculation '''
    xi1[0]=random.uniform(0, 1) #xi1
    if xi1[0]==0:
        print("change xi-1: {}".format(xi1[0]))
        xi1[0]=1-xi1[0]

    bin_xi_1=int(xi1[0]*nbins/1.)
    if 0<bin_xi_1<nbins:
        xi1_list[bin_xi_1]+=1
        tree_xi_1.Fill()

    ''' xi2 calculation '''
    xi2[0]=random.uniform(0, 1)
    if xi2[0]==0:
        print("change xi-2: {}".format(xi2[0]))
        xi2[0]=1-xi2[0]

    bin_xi_2=int(xi2[0]*nbins/1.)
    if 0<bin_xi_2<nbins:
        xi2_list[bin_xi_2]+=1
        tree_xi_2.Fill()

    ''' nu1 calculation '''
    nu1[0]=np.sqrt(-2*np.log(xi1[0]))*np.cos(2*np.pi*xi2[0])
    bin_nu_1=int(abs(-4-nu1[0])*nbins/8.)
    if 0<bin_nu_1<nbins:
        if DBG: print("nu1: {} bin1: {}".format(nu1,bin))
        nu1_list[bin_nu_1]+=1
        tree_nu_1.Fill()

    ''' nu2 calculation '''
    nu2[0]=np.sqrt(-2*np.log(xi1[0]))*np.sin(2*np.pi*xi2[0])
    bin_nu_2=int(abs(-4-nu2[0])*nbins/8.)
    if 0<bin_nu_2<nbins:
        if DBG: print("nu2: {} bin2: {}".format(nu2[0],bin))
        nu2_list[bin_nu_2]+=1
        tree_nu_2.Fill()

    #g=g_distribution(nu1,nu2)
    #g_list[bin]+=1
    #   print("xi1: {} xi2: {} nu1: {} nu2: {}".format(xi1,xi2,nu1,nu2))


print("DONE")
for i in range(nbins):
    if DBG: print("bin: {} xi1: {} nu1: {}".format(i,xi1_list[i], nu1_list[i]))
    n[0]=i
    xi1[0]=xi1_list[i]
    xi2[0]=xi2_list[i]
    nu1[0]=nu1_list[i]
    nu2[0]=nu2_list[i]
    tree.Fill()
    file_dat.write("{}\t{}\t{}\t{}\t{}\n".format(n[0],xi1[0],xi2[0],nu1[0],nu2[0]))


#step2: write to disk
tree.Write()
tree_xi_1.Write()
tree_xi_2.Write()
tree_nu_1.Write()
tree_nu_2.Write()

file.Close()
file_dat.close()
