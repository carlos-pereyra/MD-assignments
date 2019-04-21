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
#tree.Scan("box.atoms.KE_PE_V_AVG:box.atoms.KE_PE_U_AVG:box.atoms.KE_PE_V_SUM:box.atoms.KE_PE_U_SUM:box.atoms.omega_dt","","colsize=18")

for i, event in enumerate(tree):
    print("{}".format(event))

c=TCanvas("c", "canvas", 500, 500)
c.GetFrame().SetBorderSize(12)
c.cd()

#style=TStyle("Plain","Default Style")

tree.SetTitle('Omega_dt')
tree.Draw("box.atoms.KE_PE_V_AVG:box.atoms.dt","dt<2.1","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetTitle('')
g1.GetXaxis().SetTitle('#Omega_{0} dt')
g1.GetYaxis().SetTitle("Average Total Energy")
g1.GetHistogram().SetMaximum(15);
g1.Draw("A*L")

tree.SetLineStyle(1)
tree.SetLineColor(2)
tree.SetLineWidth(4)
tree.SetMarkerColor(1)
tree.SetMarkerStyle(2)
tree.Draw("box.atoms.KE_PE_U_AVG:box.atoms.dt","dt<2.1","L same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.GetHistogram().SetMaximum(15);
g2.Draw("same")

#gROOT.SetStyle("ATLAS");
c.UseCurrentStyle()
c.SetGrid()
c.Modified()
c.Update()
text=raw_input()
