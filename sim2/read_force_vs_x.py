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
    std::vector<float> u_wall;\
    std::vector<float> u_lj;\
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

file=TFile("test/particle_box_rx5_ry5_vx2_vy3.root", "read")
tree=file.Get("time-sequence")

for i, event in enumerate(tree):
    print("{}".format(event))

c=TCanvas("c", "canvas 1", 600, 500)
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderMode(0);
c.GetFrame().SetBorderSize(0);
c.SetLeftMargin(0.15)
c.SetRightMargin(0.15)
c.SetBottomMargin(0.15)
c.SetTopMargin(0.06)

#c.cd()
legend = TLegend(0.70, 0.20, .98, 0.38)
'''
#graph1
'''
c.cd()
tree.Draw("box.atoms.f.x:box.atoms.r.t","id==0 && t<20 && box.atoms.f.x>-300","")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g1.SetName("g1")
g1.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
                                                        ,tree.box.atoms[0].ry0
                                                        ,tree.box.atoms[0].vx0
                                                        ,tree.box.atoms[0].vy0))
#g1.SetMarkerStyle(2)
g1.SetLineWidth(4)
g1.GetXaxis().SetTitle('Time')
#g1.GetYaxis().SetTitle('wall force_x')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.Draw("LA PLC PFC")
gStyle.SetOptStat(0)
gStyle.SetPalette(kAurora)
gStyle.SetOptTitle(0)
#gROOT.SetStyle("ATLAS");
c.Modified()
c.Update()
'''
gPad.Print("images/fx_x_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0
                                                                    ,tree.box.atoms[0].ry0
                                                                    ,tree.box.atoms[0].vx0
                                                                    ,tree.box.atoms[0].vy0))
'''

'''
#canvas 2
'''
#d=TCanvas("d", "canvas 2", 600, 500)
#d.GetFrame().SetBorderMode(0);
#d.GetFrame().SetBorderMode(0);
#d.GetFrame().SetBorderSize(0);
#d.SetLeftMargin(0.15)
#d.SetRightMargin(0.15)
#d.SetBottomMargin(0.15)
#d.SetTopMargin(0.06)
#d.cd()

tree.Draw("box.atoms.r.x:box.atoms.r.t","id==0 && t<20","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetLineWidth(4)
g2.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
                                                        ,tree.box.atoms[0].ry0
                                                        ,tree.box.atoms[0].vx0
                                                        ,tree.box.atoms[0].vy0))
g2.Draw("L PLC PFC")
#gROOT.SetStyle("ATLAS");

legend.SetTextSize(0.04)
legend.AddEntry("g1","Wall Force-X","l")
legend.AddEntry("g2","#xi","l")
legend.Draw()
c.Modified()
c.Update()

gPad.Print("images/fx_x_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0
,tree.box.atoms[0].ry0
,tree.box.atoms[0].vx0
,tree.box.atoms[0].vy0))


text=raw_input()
