"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

nlines=5
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from ROOT import TH1F,TStyle,kYellow

file=TFile("data/HMO_dat.root", "read")
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
    tree.Draw("r:t","id==20 && t<20","")
    g0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g0.SetName("g0")
    g0.SetTitle("Wdt=%f" % (tree_energy.Wdt))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(3)
    g0.GetXaxis().SetTitle('#xi')
    g0.GetYaxis().SetTitle('#nu')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree.Draw("r:t","id==1 && t<20","same")
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(3)
    g1.Draw("L PLC PFC")

if nlines>=3:
    tree.Draw("r:t","id==2 && t<20","same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g2.SetName("g2")
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(3)
    g2.Draw("L PLC PFC")

if nlines>=4:
    tree.Draw("r:t","id==3 && t<20","same")
    g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g3.SetName("g3")
    g3.SetMarkerStyle(2)
    g3.SetLineWidth(3)
    g3.Draw("L PLC PFC")

if nlines>=5:
    tree.Draw("r:t","id==36 && t<20","same")
    g5=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g5.SetName("g1")
    g5.SetMarkerStyle(2)
    g5.SetLineWidth(3)
    g5.Draw("L PLC PFC")

#pad2.cd()

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

#pad1.SetGrid()
#pad2.SetGrid()
c1.Modified()
c1.Update()

text=raw_input()

#gPad.Print("images/n1_n2_hist_n.pdf")
