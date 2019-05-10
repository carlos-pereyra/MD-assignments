"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

nlines=1
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall

file=TFile("data/particle_in_box_data.root", "read")
tree=file.Get("time_sequence")
tree2=file.Get("system_energy")
tree3=file.Get("avg_system_energy")
tree4=file.Get("system")

#tree.Scan("x:y:ke:u_lj:u_wall","id==9")    #local issues?
#tree2.Scan("u_wall_avg:u_lj_avg:dt:t","")  #energy vs time
#tree3.Scan("","")                          #energy vs dt
#tree4.Scan()

tree.GetEntry()
tree2.GetEntry()
tree3.GetEntry()
tree4.GetEntry()

''' Canvas '''
c1=TCanvas("c", "canvas 1", 600, 500)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.10)
legend = TLegend(0.60, 0.80, .98, 0.98)

if nlines>=1:
    tree.Draw("y:x","id==0 && x<20 && x>0 && y<20 && y>0 && 0.1<dt<0.02","")
    g0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g0.SetName("g0")
    g0.SetTitle("n-atoms=%d, vx_0=%d, vy_0=%d" % (tree4.natom,tree.vx,tree.vy))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(3)
    g0.GetXaxis().SetTitle('#xi')
    g0.GetYaxis().SetTitle('#nu')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree.Draw("y:x","id==1 && 0.009<dt<0.011","same") # && 0.019<dt<0.021
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g1.SetLineWidth(3)
    g1.Draw("L PLC PFC")

if nlines>=3:
    tree.Draw("y:x","id==2 && 0.009<dt<0.011","same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g2.SetLineWidth(3)
    g2.Draw("L PLC PFC")

if nlines>=4:
    tree.Draw("y:x","id==3 && 0.009<dt<0.011","same")
    g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g3.SetLineWidth(3)
    g3.Draw("L PLC PFC")

if nlines>=5:
    tree.Draw("y:x","id==4 && 0.009<dt<0.011","same")
    g4=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g4.SetLineWidth(3)
    g4.Draw("L PLC PFC")

if nlines>=6:
    tree.Draw("y:x","id==5 && 0.009<dt<0.011","same")
    g5=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g5.SetLineWidth(3)
    g5.Draw("L PLC PFC")

if nlines>=7:
    tree.Draw("y:x","id==6 && 0.009<dt<0.011","same")
    g6=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g6.SetLineWidth(3)
    g6.Draw("L PLC PFC")

if nlines>=8:
    tree.Draw("y:x","id==7 && 0.009<dt<0.011","same")
    g7=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g7.SetLineWidth(3)
    g7.Draw("L PLC PFC")

if nlines>=9:
    tree.Draw("y:x","id==8 && 0.009<dt<0.011","same")
    g8=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g8.SetLineWidth(3)
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

text=raw_input()

gPad.Print("images/particle_ntrace{}_nAtoms{}_tSteps{}_dtSteps{}_TEST3.pdf"
           .format(nlines,tree4.natom,tree4.ntime,tree4.ndt))
