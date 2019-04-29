"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
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
    std::vector<float> ke_avg;\
    std::vector<float> u_wall_avg;\
    std::vector<float> u_lj_avg;\
    std::vector<float> total_energy_avg;\
    Int_t id;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
};");

file=TFile("test/particle_box_rx5_ry5_vx2_vy3.root", "read")
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
c1.SetTopMargin(0.06)

#c.cd()
legend = TLegend(0.60, 0.80, .98, 0.98)
'''
#graph1
'''
c1.cd()
tree.Draw("box.atoms.r.y:box.atoms.r.x","id==0 && x<40 && y<40 && y>0","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
g1.SetTitle('(x0,y0)=({:.0f},{:.0f}), (vx0,vy0)=({:.0f},{:.0f})'.format(tree.box.atoms[0].rx0,tree.box.atoms[0].ry0,tree.box.atoms[0].vx0,tree.box.atoms[0].vy0))
g1.SetMarkerStyle(2)
g1.SetLineWidth(4)
g1.GetXaxis().SetTitle('#xi')
g1.GetYaxis().SetTitle('#nu')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.Draw("LA PLC PFC")

'''
tree.Draw("box.atoms.r.y:box.atoms.r.x","id==1","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetLineWidth(4)
g2.Draw("L PLC PFC")


tree.Draw("box.atoms.r.y:box.atoms.r.x","id==2","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g3.SetLineWidth(4)
g3.Draw("L PLC PFC")

tree.Draw("box.atoms.r.y:box.atoms.r.x","id==3","same")
g4=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g4.SetLineWidth(4)
g4.Draw("L PLC PFC")

tree.Draw("box.atoms.r.y:box.atoms.r.x","id==4","same")
g5=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g5.SetLineWidth(4)
g5.Draw("L PLC PFC")
'''

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0)
#c1.SetGrid()
c1.Modified()
c1.Update()

#legend.SetTextSize(0.04)
#legend.AddEntry("g1","#Omega_0 dt={:.2e}".format(tree.box.atoms[1].omega_dt),"l")
#legend.AddEntry("g2","#Omega_0 dt={:.2e}".format(tree.box.atoms[10].omega_dt),"l")
#legend.AddEntry("g3","#Omega_0 dt={:.2e}".format(tree.box.atoms[20].omega_dt),"l")
#legend.AddEntry("g4","#Omega_0 dt={:.2e}".format(tree.box.atoms[30].omega_dt),"l")
#legend.Draw()

#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");


gPad.Print("images/xy_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0,tree.box.atoms[0].ry0,tree.box.atoms[0].vx0,tree.box.atoms[0].vy0))

text=raw_input()
