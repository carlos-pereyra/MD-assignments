"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import *
import numpy as np
from array import array

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
    position r;\
    velocity v;\
    half u;\
    force f;\
    Float_t k;\
    Float_t m;\
    Float_t dt;\
    Int_t id;\
    Float_t KE_V;\
    Float_t KE_U;\
    Float_t PE_R;\
    Float_t KE_PE_V_SUM;\
    Float_t KE_PE_U_SUM;\
    Float_t KE_PE_V_AVG;\
    Float_t KE_PE_U_AVG;\
    Int_t ncycle;\
    Float_t freq;\
    Float_t omega_dt;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

file=TFile("dt-stormer-verlet.root", "read");
tree=file.Get("time-sequence")

for i, event in enumerate(tree):
    print("{}".format(event))

#tree.Scan("box.atoms.KE_PE_V_AVG:box.atoms.dt:box.atoms.freq:box.atoms.omega_dt","","colsize=20")

c=TCanvas("c", "canvas", 500, 500)
c.GetFrame().SetBorderSize(12)
c.cd()

tree.SetLineStyle(1)
tree.SetLineColor(2)
tree.SetLineWidth(4)
tree.SetMarkerColor(1)
tree.SetMarkerStyle(2)
tree.SetTitle('Omega_dt')
tree.Draw("box.atoms.KE_PE_V_AVG:box.atoms.dt","dt<2","L")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetTitle('')
g1.GetXaxis().SetTitle('Omega_dt')
g1.GetYaxis().SetTitle("Average Total Energy")
g1.Draw("A*L")

c.SetGrid()
c.Modified()
c.Update()
text=raw_input()

#print options
#print("Tree content:\n{}\n".format(np.asarray([tree.box.atoms[i].omega_dt for i, event in enumerate(tree.atoms)])))

#if __name__ == "__main__":
