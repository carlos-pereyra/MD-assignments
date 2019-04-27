"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle
from ROOT import TCanvas,TGraph,TPad
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall

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
    Float_t rx0;\
    Float_t ry0;\
    Float_t vx0;\
    Float_t vy0;\
    position r;\
    velocity v;\
    half u;\
    force f;\
    Float_t k;\
    Float_t m;\
    Float_t dt;\
    std::vector<float> pe;\
    Float_t KE_V;\
    Float_t KE_U;\
    Float_t KE_PE_V_SUM;\
    Float_t KE_PE_U_SUM;\
    Float_t KE_PE_V_AVG;\
    Float_t KE_PE_U_AVG;\
    Int_t ncycle;\
    Float_t freq;\
    Float_t omega_dt;\
    Int_t id;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

file=TFile("data/stormer-verlet-dat-rx3-vx-1.root", "read")
tree=file.Get("time-sequence")
#tree.Scan("box.atoms.pe_x","id==1","colsize=18")

for i, event in enumerate(tree):
    print("{}".format(event))

c1=TCanvas("c", "canvas 1", 600, 500)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.04)

'''
c2=TCanvas("c2", "canvas 2", 600, 500)
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderSize(0);
c2.SetLeftMargin(0.15)
c2.SetRightMargin(0.15)
c2.SetBottomMargin(0.15)
c2.SetTopMargin(0.04)

c3=TCanvas("c3", "canvas 3", 600, 500)
c3.GetFrame().SetBorderMode(0);
c3.GetFrame().SetBorderMode(0);
c3.GetFrame().SetBorderSize(0);
c3.SetLeftMargin(0.15)
c3.SetRightMargin(0.15)
c3.SetBottomMargin(0.15)
c3.SetTopMargin(0.04)
'''

#c.cd()
legend = TLegend(0.60, 0.80, .98, 0.98)
'''
#graph1
'''
c1.cd()
tree.Draw("box.atoms.r.y:box.atoms.r.x","id==0 && x>0 && x<40 && y>0 && y<40","L PLC PFC")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
#g1.SetMarkerColor(1)
#g1.SetMarkerStyle(2)
#g1.SetLineWidth(2)
#g1.GetXaxis().SetTitle('t')
#g1.GetYaxis().SetTitle('r_x(t)')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.Draw()

#c.SetGrid()

tree.Draw("box.atoms.r.y:box.atoms.r.x","id==1 && x>0 && x<40 && y>0 && y<40","same L PLC PFC")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.Draw()

gStyle.SetPalette(kRainBow);
c1.Modified()
c1.Update()

'''
#graph2
'''
'''
c2.cd()
tree.Draw("box.atoms.f.x:box.atoms.r.x","id==0 && x>0 && x<40 && y>0 && y<40","")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetMarkerColor(4)
g2.SetMarkerStyle(2)
#g2.SetLineStyle(4)
g2.SetLineWidth(2)
#g2.SetLineColor(4)
g2.Draw("L PLC PFC")

#c.SetGrid()
c2.Modified()
c2.Update()


'''
#graph3
'''
c3.cd()
#graph2
tree.Draw("box.atoms.f.y:box.atoms.r.x","id==0 && x>0 && x<40 && y>0 && y<40","")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetMarkerColor(4)
g2.SetMarkerStyle(2)
#g2.SetLineStyle(4)
g2.SetLineWidth(2)
#g2.SetLineColor(4)
g2.Draw("L PLC PFC")

#c.SetGrid()
c3.Modified()
c3.Update()
'''

'''
#graph3
tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==20 && t<20","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g3.SetName("g3")
g3.SetMarkerColor(4)
g3.SetMarkerStyle(2)
#g3.SetLineStyle(6)
g3.SetLineWidth(2)
#g3.SetLineColor(6)
g3.Draw("C PLC PFC")

#graph4
tree.Draw("box.atoms.r.x:box.atoms.t>>g2","id==30 && t<20","same")
g4=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g4.SetName("g4")
g4.SetMarkerColor(4)
g4.SetMarkerStyle(2)
#g4.SetLineStyle(6)
g4.SetLineWidth(2)
#g4.SetLineColor(6)
g4.Draw("C PLC PFC")

legend.SetTextSize(0.04)
legend.AddEntry("g1","#Omega_0 dt={:.2e}".format(tree.box.atoms[1].omega_dt),"l")
legend.AddEntry("g2","#Omega_0 dt={:.2e}".format(tree.box.atoms[10].omega_dt),"l")
legend.AddEntry("g3","#Omega_0 dt={:.2e}".format(tree.box.atoms[20].omega_dt),"l")
legend.AddEntry("g4","#Omega_0 dt={:.2e}".format(tree.box.atoms[30].omega_dt),"l")

#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);
#gStyle.SetPalette(kRust);
gROOT.SetStyle("ATLAS");

legend.Draw()
'''

#gPad.Print("images/r(t)_r{}_v{}.pdf".format(tree.box.atoms[1].r0,tree.box.atoms[1].v0))

text=raw_input()
