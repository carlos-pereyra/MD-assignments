"""
*  author: carlos p
*  purpose: plot HMO r(t) vs t
*
"""

nlines=4
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from ROOT import TH1F,TStyle,kYellow

file=TFile("data/HMO_dat_r0_v1.root", "read")
tree=file.Get("time")
tree_freq=file.Get("freq")
tree_energy=file.Get("energy")

#tree.Scan("r:t:dt","","")
tree.GetEntry()
tree_freq.GetEntry()
tree_energy.GetEntry()


''' Canvas-1 '''
c1=TCanvas("c1", "canvas_1", 850, 400)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.1)
c1.SetRightMargin(0.1)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.1)
legend = TLegend(0.64, 0.20, .97, 0.37)

#pad1=TPad("pad1","xi_1_histogram",0,0.0,0.51,0.98);
#pad1.SetLeftMargin(0.1)
#pad1.SetRightMargin(0.1)
#pad1.SetBottomMargin(0.15)
#pad1.SetTopMargin(0.1)
#pad1.Draw()

#pad2=TPad("pad2","xi_2_histogram",0.51,0.0,0.98,0.98)
#pad2.SetLeftMargin(0.1)
#pad2.SetRightMargin(0.1)
#pad2.SetBottomMargin(0.15)
#pad2.SetTopMargin(0.1)
#pad2.Draw()

#pad1.cd()
if nlines>=1:
    tree.Draw("r:t","id==0 && t<20","")
    g0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g0.SetName("g0")
    g0.SetTitle("r0=%.0f, v0=%.0f" % (tree.r0,tree.v0))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(3)
    g0.GetXaxis().SetTitle('time')
    g0.GetYaxis().SetTitle('r')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree.Draw("r:t","id==10 && t<20","same")
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(3)
    g1.Draw("L PLC PFC")

if nlines>=3:
    tree.Draw("r:t","id==15 && t<20","same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g2.SetName("g2")
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(3)
    g2.Draw("L PLC PFC")

if nlines>=4:
    tree.Draw("r:t","id==20 && t<20","same")
    g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g3.SetName("g3")
    g3.SetMarkerStyle(2)
    g3.SetLineWidth(3)
    g3.Draw("L PLC PFC")

if nlines>=5:
    tree.Draw("r:t","id==25 && t<20","same")
    g4=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g4.SetName("g4")
    g4.SetMarkerStyle(2)
    g4.SetLineWidth(3)
    g4.Draw("L PLC PFC")

#pad2.cd()
legend.SetTextSize(0.04)
#legend.AddEntry("g0","#Omega_0 dt={:.2e}".format(w0),"l")
#legend.AddEntry("g1","#Omega_0 dt={:.2e}".format(w1),"l")
#legend.AddEntry("g2","#Omega_0 dt={:.2e}".format(w2),"l")
#legend.AddEntry("g3","#Omega_0 dt={:.2e}".format(w3),"l")

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

#pad1.SetGrid()
#pad2.SetGrid()
#legend.Draw()
c1.Modified()
c1.Update()

text=raw_input()

gPad.Print("images/HMO_r%.0f_v%.0f.pdf" % (tree.r0,tree.v0))
