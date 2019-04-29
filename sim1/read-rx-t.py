"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall

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
    std::vector<float> KE_V;\
    std::vector<float> KE_U;\
    std::vector<float> PE_R;\
    Float_t KE_PE_V_SUM;\
    Float_t KE_PE_U_SUM;\
    Float_t KE_PE_V_AVG;\
    Float_t KE_PE_U_AVG;\
    Int_t ncycle;\
    Float_t freq;\
    Float_t freq_intpl;\
    Float_t omega_dt;\
    Int_t r0;\
    Int_t v0;\
    Float_t t1;\
    Float_t t2;\
    Int_t size;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

file=TFile("data/stormer-verlet-dat-r0-v1.root", "read")
tree=file.Get("time-sequence")
#tree.Scan("box.atoms.r.x:box.atoms.t","id==1","colsize=18")

for i, event in enumerate(tree):
    print("{}".format(event))

c=TCanvas("c", "canvas", 600, 500)
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderSize(0);
c.SetLeftMargin(0.15)
c.SetRightMargin(0.15)
c.SetBottomMargin(0.15)
c.SetTopMargin(0.06)

legend = TLegend(0.64, 0.20, .97, 0.37)

#graph1
tree.Draw("box.atoms.r.x:box.atoms.t","id==0 && t<t2+5","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
g1.SetTitle('r0={} v0={}'.format(tree.box.atoms[0].r0,tree.box.atoms[0].v0))
g1.SetMarkerColor(1)
g1.SetMarkerStyle(2)
g1.SetLineWidth(4)
g1.GetXaxis().SetTitle('t')
g1.GetYaxis().SetTitle('r(t)')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.Draw("LA PLC PFC")

#graph2
tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==10 && t<t2+5","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetMarkerColor(4)
g2.SetMarkerStyle(2)
g2.SetLineWidth(4)
g2.SetLineStyle(2)
#g2.SetLineColor(4)
g2.Draw("L PLC PFC")

#graph3
tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==20 && t<t2+5","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g3.SetName("g3")
g3.SetMarkerColor(4)
g3.SetMarkerStyle(2)
g3.SetLineStyle(2)
g3.SetLineWidth(4)
#g3.SetLineColor(6)
g3.Draw("L PLC PFC")

#graph4
tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==30 && t<t2+5","same")
g4=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g4.SetName("g4")
g4.SetMarkerColor(4)
g4.SetMarkerStyle(2)
g4.SetLineStyle(2)
g4.SetLineWidth(4)
#g4.SetLineColor(6)
g4.Draw("L PLC PFC")

legend.SetTextSize(0.04)
legend.AddEntry("g1","#Omega_0 dt={:.2e}".format(tree.box.atoms[1].omega_dt),"l")
legend.AddEntry("g2","#Omega_0 dt={:.2e}".format(tree.box.atoms[10].omega_dt),"l")
legend.AddEntry("g3","#Omega_0 dt={:.2e}".format(tree.box.atoms[20].omega_dt),"l")
legend.AddEntry("g4","#Omega_0 dt={:.2e}".format(tree.box.atoms[30].omega_dt),"l")

gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0)
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

legend.Draw()
#c.SetGrid()
c.Modified()
c.Update()
gPad.Print("images/r(t)_r{}_v{}.pdf".format(tree.box.atoms[1].r0,tree.box.atoms[1].v0))

text=raw_input()
