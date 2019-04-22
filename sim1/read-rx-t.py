"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import ROOT
from ROOT import TTree,TFile,gROOT,gStyle
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend

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

file=TFile("stormer-verlet-dat.root", "read");
tree=file.Get("time-sequence")
#tree.Scan("box.atoms.r.x:box.atoms.t","id==1","colsize=18")

for i, event in enumerate(tree):
    print("{}".format(event))

c=TCanvas("c", "canvas", 500, 500)
c.GetFrame().SetBorderSize(12)
c.cd()

legend = TLegend()

tree.Draw("box.atoms.r.x:box.atoms.t","id==1 && t<40","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
g1.SetMarkerColor(1)
g1.SetMarkerStyle(2)
g1.SetLineStyle(1)
g1.SetLineWidth(2)
#
g1.SetTitle('')
g1.GetXaxis().SetTitle('t')
g1.GetYaxis().SetTitle('r(t)')
g1.Draw("AL")

tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==10 && t<40","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetMarkerColor(4)
g2.SetMarkerStyle(2)
g2.SetLineStyle(4)
g2.SetLineWidth(2)
g2.SetLineColor(4)
g2.Draw("same")

tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==20 && t<40","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g3.SetName("g3")
g3.SetMarkerColor(4)
g3.SetMarkerStyle(2)
g3.SetLineStyle(6)
g3.SetLineWidth(2)
g3.SetLineColor(6)
g3.Draw("same")

gStyle.SetOptStat(0)
legend.AddEntry("g1","id=1","l")
legend.AddEntry("g2","id=10","l")
legend.AddEntry("g3","id=20","l")
legend.Draw()

#c.SetGrid()
c.Modified()
c.Update()
text=raw_input()
