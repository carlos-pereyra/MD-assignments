"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

nlines=5
import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from array import array

gROOT.ProcessLine("struct Particle {\
std::vector<float> t;\
std::vector<float> ux;\
std::vector<float> uy;\
std::vector<float> x;\
std::vector<float> y;\
std::vector<float> fx;\
std::vector<float> fy;\
std::vector<float> vx;\
std::vector<float> vy;\
Float_t k;\
Float_t m;\
Float_t dt;\
Int_t id;\
std::vector<float> l_mv;\
std::vector<float> u_wl;\
std::vector<float> u_lj;\
std::vector<float> total_u;\
std::vector<float> u_per_atom;\
Float_t            u_per_time;\
};\
class Box{\
public:\
std::vector<Particle> atoms;\
std::vector<float> u_per_atom;\
};");

file=TFile("data/test.root", "read")
tree=file.Get("time_sequence")
tree_energy_t=file.Get("system_energy")         #new tree
tree_energy_dt=file.Get("avg_system_energy")    #new tree

#tree.Scan("x:y:vx:vy:t","id==1 && dt>0.009 && dt<0.011","colsize=14")
#for i, event in enumerate(tree):
    #print("{}".format(event))

''' Canvas '''
c1=TCanvas("c", "canvas 1", 600, 500)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.06)
legend = TLegend(0.60, 0.80, .98, 0.98)

if nlines>=1:
    tree.Draw("y:x","id==0 && 0.009<dt<0.011","")
    g0=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g0.SetName("g0")
    #g1.SetTitle('(x0,y0)=({:.0f},{:.0f}), (vx0,vy0)=({:.0f},{:.0f})'.format(tree.box.atoms[0].rx0,tree.box.atoms[0].ry0,tree.box.atoms[0].vx0,tree.box.atoms[0].vy0))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(4)
    g0.GetXaxis().SetTitle('#xi')
    g0.GetYaxis().SetTitle('#nu')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree.Draw("y:x","id==1 && 0.009<dt<0.011","same") # && 0.019<dt<0.021
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g1.SetLineWidth(4)
    g1.Draw("L PLC PFC")

if nlines>=3:
    tree.Draw("y:x","id==2 && 0.009<dt<0.011","same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g2.SetLineWidth(4)
    g2.Draw("L PLC PFC")

if nlines>=4:
    tree.Draw("y:x","id==3 && 0.009<dt<0.011","same")
    g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g3.SetLineWidth(4)
    g3.Draw("L PLC PFC")

if nlines>=5:
    tree.Draw("y:x","id==4 && 0.009<dt<0.011","same")
    g4=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g4.SetLineWidth(4)
    g4.Draw("L PLC PFC")

if nlines>=6:
    tree.Draw("y:x","id==5 && 0.009<dt<0.011","same")
    g5=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g5.SetLineWidth(4)
    g5.Draw("L PLC PFC")

if nlines>=7:
    tree.Draw("y:x","id==6 && 0.009<dt<0.011","same")
    g6=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g6.SetLineWidth(4)
    g6.Draw("L PLC PFC")

if nlines>=8:
    tree.Draw("y:x","id==7 && 0.009<dt<0.011","same")
    g7=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g7.SetLineWidth(4)
    g7.Draw("L PLC PFC")

if nlines>=9:
    tree.Draw("y:x","id==8 && 0.009<dt<0.011","same")
    g8=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
    g8.SetLineWidth(4)
    g8.Draw("L PLC PFC")


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


#gPad.Print("images/xy_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0,tree.box.atoms[0].ry0,tree.box.atoms[0].vx0,tree.box.atoms[0].vy0))

text=raw_input()
