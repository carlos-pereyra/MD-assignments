"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import ROOT
from ROOT import TTree,TFile,gStyle,gROOT,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend,kRainBow,kRust,kFall

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
    Int_t size;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

file=TFile("data/stormer-verlet-dat-r1-v0.root", "read");
tree=file.Get("time-sequence")
#tree.Scan("box.atoms.dt:box.atoms.omega_dt","","colsize=18")
#tree.Scan("(Sum$(box.atoms.KE_V)+Sum$(box.atoms.PE_R))/Length$(box.atoms.KE_V)","id==0","colsize=20")
#tree.Scan("(Sum$(box.atoms.KE_V)+Sum$(box.atoms.PE_R))/Length$(box.atoms.KE_V):box.atoms.id","","colsize=20")
#tree.Scan("box.atoms.KE_V+box.atoms.PE_R","id==0","colsize=25")
#tree.Scan("box.atoms.KE_PE_V_SUM:box.atoms.n","id==0","colsize=25")
#tree.Scan("box.atoms.KE_PE_V_SUM:box.atoms.size:box.atoms.KE_PE_V_SUM/size:omega_dt","id==0","colsize=25")
for i, event in enumerate(tree):
    print("{}".format(event))

legend = TLegend(0.35, 0.17, .97, 0.32)
c=TCanvas("c", "canvas", 500, 500)
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderSize(0);
c.SetLeftMargin(0.15)
c.SetRightMargin(0.15)
c.SetBottomMargin(0.15)
c.SetTopMargin(0.06)

tree.Draw("box.atoms.KE_PE_V_SUM/size:omega_dt","omega_dt<2.1","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
g1.SetTitle('r0={} v0={}'.format(tree.box.atoms[0].r0,tree.box.atoms[0].v0))
#g2.SetLineStyle(1)
#g2.SetLineColor(2)
g1.SetLineWidth(4)
#g2.SetMarkerColor(1)
#g2.SetMarkerStyle(2)
g1.GetXaxis().SetTitle('#Omega_{0}dt')
g1.GetYaxis().SetTitle('Total Energy_{avg}')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.GetHistogram().SetMaximum(1);
g1.GetHistogram().SetMinimum(0);
g1.Draw("LA PLC")


tree.Draw("box.atoms.KE_PE_U_SUM/size:omega_dt","omega_dt<2.1","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
#g2.SetLineStyle(1)
#g2.SetLineColor(2)
g2.SetLineWidth(4)
#g2.SetMarkerColor(4)
#g2.SetMarkerStyle(2)
g2.Draw("L PLC PFC")

legend.SetTextSize(0.04)
legend.AddEntry("g1","#sum(PE+KE(v)), r0={:.1f}, v0={:.1f}".format(tree.box.atoms[0].r0,tree.box.atoms[0].v0),"l")
legend.AddEntry("g2","#sum(PE+KE(u)), r0={:.1f}, v0={:.1f}".format(tree.box.atoms[0].r0,tree.box.atoms[0].v0),"l")
legend.Draw()

#gROOT.SetStyle("ATLAS");
gStyle.SetOptStat(0)
gStyle.SetPalette(kFall)
#gStyle.SetOptTitle(0)

#c.UseCurrentStyle()
c.SetGrid()
c.Modified()
c.Update()
gPad.Print("images/total-energy_r{}_v{}.pdf".format(tree.box.atoms[0].r0,tree.box.atoms[0].v0))
text=raw_input()
