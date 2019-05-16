"""
*  author: carlos p
*  purpose: plot HMO frequency vs Wdt
*
"""

nlines=2
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
c1=TCanvas("c1", "canvas_1", 600, 400)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.1)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.1)
legend = TLegend(0.17, 0.70, .79, 0.87)

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
    tree_freq.Draw("f_avg:Wdt","","")
    g0=TGraph(tree_freq.GetSelectedRows(),tree_freq.GetV2(),tree_freq.GetV1())
    g0.SetName("g0")
    g0.SetTitle("r0=%.0f, v0=%.0f" % (tree.r0,tree.v0))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(3)
    g0.GetXaxis().SetTitle('#Omega_0 dt')
    g0.GetYaxis().SetTitle('Freq [Hz]')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree_freq.Draw("f_i_avg:Wdt","","same")
    g1=TGraph(tree_freq.GetSelectedRows(),tree_freq.GetV2(),tree_freq.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(3)
    g1.Draw("L PLC PFC")

#pad2.cd()

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

legend.SetTextSize(0.045)
legend.AddEntry("g0","frequency average with full step velocity","l")
legend.AddEntry("g1","frequency average with half step velocity","l")
legend.Draw()

#pad1.SetGrid()
#pad2.SetGrid()
c1.SetGrid()
c1.Modified()
c1.Update()

text=raw_input()

gPad.Print("images/HMO_freq_r%.0f_v%.0f.pdf" % (tree.r0,tree.v0))
