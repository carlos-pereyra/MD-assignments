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
from ROOT import TH1F,TStyle,kYellow

file=TFile("data/data_n100000.root", "read")
tree=file.Get("distribution")
tree_xi_1=file.Get("xi_1_distribution")
tree_xi_2=file.Get("xi_2_distribution")
tree_nu_1=file.Get("nu_1_distribution")
tree_nu_2=file.Get("nu_2_distribution")

#tree.Scan("xi1:xi2:nu1:nu2:g","")
tree.GetEntry()
tree_xi_1.GetEntry()
tree_xi_2.GetEntry()
tree_nu_1.GetEntry()
tree_nu_2.GetEntry()

''' Canvas-1 '''
c1=TCanvas("c1", "canvas_1", 850, 400)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.01)
c1.SetRightMargin(0.01)
c1.SetBottomMargin(0.01)
c1.SetTopMargin(0.01)

pad1=TPad("pad1","xi_1_histogram",0,0.0,0.51,0.98);
pad1.SetLeftMargin(0.1)
pad1.SetRightMargin(0.1)
pad1.SetBottomMargin(0.15)
pad1.SetTopMargin(0.1)
pad1.Draw()

pad2=TPad("pad2","xi_2_histogram",0.51,0.0,0.98,0.98)
pad2.SetLeftMargin(0.1)
pad2.SetRightMargin(0.1)
pad2.SetBottomMargin(0.15)
pad2.SetTopMargin(0.1)
pad2.Draw()

pad1.cd()
h1=TH1F("h1","xi_1",1000,0,1.1)
tree_xi_1.Draw("xi1>>h1","","HIST PFC")
h1.GetXaxis().SetTitle('#xi 1')
#h1.GetYaxis().SetTitle('Counts')
h1.GetXaxis().SetTitleSize(0.06)
h1.GetYaxis().SetTitleSize(0.06)
h1.SetFillColor(40)
h1.SetFillStyle(3008)
h1.Draw()

pad2.cd()
h2=TH1F("h2","xi_2",1000,0,1.1)
tree_xi_2.Draw("xi2>>h2","","HIST PFC")
h2.GetXaxis().SetTitle('#xi 2')
#h2.GetYaxis().SetTitle('Counts')
h2.GetXaxis().SetTitleSize(0.06)
h2.GetYaxis().SetTitleSize(0.06)
h2.SetFillColor(40)
h2.SetFillStyle(3008)
h2.Draw()

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0)
pad1.SetGrid()
pad2.SetGrid()
c1.Modified()
c1.Update()

#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");
#gPad.Print("images/n1_n2_hist_d.pdf")

''' Canvas-2 '''
c2=TCanvas("c2", "canvas_2", 850, 400)
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderSize(0);
c2.SetLeftMargin(0.01)
c2.SetRightMargin(0.01)
c2.SetBottomMargin(0.01)
c2.SetTopMargin(0.01)

pad3=TPad("pad3","nu_1_histogram",0,0.0,0.51,0.98);
pad3.SetLeftMargin(0.1)
pad3.SetRightMargin(0.1)
pad3.SetBottomMargin(0.15)
pad3.SetTopMargin(0.1)
pad3.Draw()

pad4=TPad("pad4","nu_2_histogram",0.51,0.0,0.98,0.98)
pad4.SetLeftMargin(0.1)
pad4.SetRightMargin(0.1)
pad4.SetBottomMargin(0.15)
pad4.SetTopMargin(0.1)
pad4.Draw()

pad3.cd()
h3=TH1F("h3","nu_1",1000,-4,4)
tree_nu_1.Draw("nu1>>h3","","HIST PFC")
h3.GetXaxis().SetTitle('#eta 1')
#h3.GetYaxis().SetTitle('Counts')
h3.GetXaxis().SetTitleSize(0.06)
h3.GetYaxis().SetTitleSize(0.06)
h3.SetFillColor(40)
h3.SetFillStyle(3008)
h3.Draw()

pad4.cd()
h4=TH1F("h4","nu_2",1000,-4,4)
tree_nu_2.Draw("nu2>>h4","","HIST PFC")
h4.GetXaxis().SetTitle('#eta 2')
#h4.GetYaxis().SetTitle('Counts')
h4.GetXaxis().SetTitleSize(0.06)
h4.GetYaxis().SetTitleSize(0.06)
h4.SetFillColor(40)
h4.SetFillStyle(3008)
h4.Draw()

c2.Modified()
c2.Update()

text=raw_input()

#gPad.Print("images/n1_n2_hist_n.pdf")
